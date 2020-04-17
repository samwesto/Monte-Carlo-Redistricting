import numpy as np
import pandas as pd
import fiona
import random
import geopandas as gpd
import csv
import shapely
import os
import itertools
import time
import MCMC_functions as MC
import scipy.stats as sct

class TimeError(Exception):
	pass





###Rodden Algorithm
###################

class RoddenChain:
	Method = 'constructive'

	def __init__(self,df,condf,pop_constraint,comp_constraint,comp=None,relabel=None,method=None,OAframe=None): #df = main ward data, OAframe = OA area for MMI compactness, comp = area or MMI depending on compactness measure, pop and comp constraints are numeric bounds for sampling. 
		self.df = df
		self.condf = condf
		if comp == 'MMI':
			self.OAframe = OAframe
			self.comp = MC.MMIcompact
		else:
			self.comp = MC.Areacompact
		if relabel == 'bynumber':
			self.relabel = MC.relabel_bynumber
		else: 
			self.relabel = MC.relabel_bylocation
		if method == 'RoddenWang':
			self.method = 'RoddenWang'
		else:
			self.method = 'Normal'
		self.pop_constraint = pop_constraint
		self.comp_constraint = comp_constraint


	def Run(self,maps,con_num):
		self.start_time = time.time()
		self.pop_constraint
		condf = self.condf
		df = self.df
		OAframe = self.OAframe

		dist = []
		no_maps = 0
		while no_maps <= maps:

			samp = False
			cycle_start = time.time()
			while True:
				int_cons,samp = RunRalg(self,con_num)

				if time.time() - cycle_start > 1800 or time.time() - self.start_time > 7200:
					raise TimeError("It looks like the chain can't find any more states")


				int_cons = self.comp(self,df,int_cons,range(len(int_cons)))

				if samp == True and any(x > self.comp_constraint for x in int_cons['comp_score']) == False:
					break

			newdist = int_cons['wards']
			no_maps +=1
			vec = []
			for x in range(len(newdist.tolist())):
				xx = newdist[x]
				for y in xx:
					vec.append((y,x))
			vec = sorted(vec, key=lambda k: k[0])
			vec = [x[1]+1 for x in vec]
			dist.append([self.relabel(dist,vec),0])
		return dist



#Creates a new independent map each time
def RunRalg(self,con_number):
	int_cons = self.condf.copy()
	int_df = self.df
	OAframe = self.OAframe

	self.avpop = sum(int_df['All Ages'])/con_number

	#Part 1 of Rodden Algorithm - using method of combing until desired number of int_constituencies
	while True: 
		con = random.choice(range(len(int_cons['con']))) #Pick rando constituency
		concent = int_cons['geometry'][con].centroid 
		neighcent = []
		for neigh in int_cons['neighbours'][con]: #Find center of each con, then compute distance from these to center of con
			neighcent.append(concent.distance(int_cons['geometry'][neigh].centroid)) 
		mergecon = int_cons['neighbours'][con][neighcent.index(min(neighcent))] #Finds the neighbour with the closest geographic distance
		int_cons.at[con,'wards'] = int_cons.at[con,'wards'] + int_cons.at[mergecon,'wards'] 
		int_cons.at[[con],'geometry'] = gpd.GeoSeries([int_cons.at[con,'geometry'],int_cons.at[mergecon,'geometry']]).unary_union #Assign new shape to merged int_constituency
		int_cons = int_cons.drop(mergecon).reset_index(drop=True) #Drop old int_constituency
		int_cons = Rneighbours(self,int_cons)
		if int_cons['con'].nunique() == con_number: #Break loop when int_constituency number is the desired number (con_number for our reference dist)
			break

	if self.method == 'RoddenWang':
			#Second part of Rodden (All Ages evening) - creating maps that meet population int_constraint
		ind = list(range(len(int_cons['pop_score']))) #Do the pop and comp scores for the new set up, only at this point does it matter.
		int_cons = MC.population(self,int_cons,self.avpop,ind)
		
		Time = time.time()
		while any(abs(x) > self.pop_constraint for x in int_cons['pop_score']):
			constit = [0]*len(df['constituency'])
			for i,wards in zip(int_cons['con'],int_cons['wards']):
				for j in wards:
					constit[j] = i 
			int_df['constituency'] = constit

			Gibbs = 0.2
			edgeframe = updateedgeframe(int_df)
			no_cons = 3
			avpop = sum(int_df['All Ages'])/no_cons
			b = 10


			q = float(0.04)
			b2 = 1
			b1 = 1
			int_df,int_cons,Gibbs = UpdateS(q,b,self.avpop,int_df,edgeframe,int_cons,Gibbs)
			
			if (time.time() - Time) > 30: #Ends search and starts a new one. Necessary due to the deterministic nature of the chain. 
				return int_cons,False

		return int_cons,True	

	else:
		#Second part of Rodden (population evening) - creating maps that meet population int_constraint
		ind = list(range(len(int_cons['pop_score'])))
		int_cons = MC.population(self,int_cons,self.avpop,ind)

		Time = time.time()
		while (sum([1 if abs(x) > self.pop_constraint else 0 for x in int_cons['pop_score']]) != 0):
			neighdata = []
			for h in range(len(int_cons['neighbours'])):
				for i in int_cons['neighbours'][h]:
					if [i,h] not in neighdata:
						neighdata.append([h,i])
			popdiff = []
			for k in neighdata:
				popdiff.append(int_cons['pop_score'][k[0]] - int_cons['pop_score'][k[1]])
			list3 = [abs(x) for x in popdiff]
			ind = list3.index(max(list3)) #Finds the max All Ages difference
			if popdiff[ind] < 0: #x is the ward with the bigger All Ages
				[i,j] = neighdata[ind]
			else:
				[j,i] = neighdata[ind]

			swaps = []
			wards = int_cons.at[j,'wards']
			for ward in wards:
				for neigh in int_df['neighbours'][ward]: #Find wards that have neighbours in other int_constituency
					if neigh in int_cons.at[i,'wards'] and type(gpd.GeoSeries(int_df['geometry'][[x for x in int_cons.at[j,'wards'] if x != ward]]).unary_union) != shapely.geometry.multipolygon.MultiPolygon: #Either got ward not going into ward or some touching error with shapes. 
						swaps.append(ward)
						break
			if len(swaps) == 0: #If can't make a move here, restart
				return int_cons,False

			jcent = int_cons['geometry'][j].centroid
			icent = int_cons['geometry'][i].centroid
			distance = []
			for ward in swaps:
				cent = int_df['geometry'][ward].centroid
				dist1 = icent.distance(cent)
				dist2 = jcent.distance(cent)
				distance.append(dist1-dist2)
			chosenward = swaps[distance.index(min(distance))]
			int_df.at[chosenward,'constituency'] == int_df['constituency'][int_cons.at[j,'wards'][0]]
			int_cons.at[i,'wards'] = int_cons.at[i,'wards'] + [chosenward]


			int_cons.at[j,'wards'] = [x for x in int_cons.at[j,'wards'] if x != chosenward]
			int_cons.at[i,'geometry'] = gpd.GeoSeries([int_cons.at[i,'geometry'],int_df.at[chosenward,'geometry']]).unary_union
			wards = int_cons.at[j,'wards']
			int_cons.at[[j],'geometry'] = gpd.GeoSeries(int_df['geometry'][wards].tolist()).unary_union
			ind = [i,j]
			int_cons = MC.population(self,int_cons,self.avpop,ind)
			int_cons = Rneighbours(self,int_cons)

			#print(int_cons['wards'])
			#print(int_cons['pop_score'])

			if (time.time() - Time) > 30: #Ends search and starts a new one. Necessary due to the deterministic nature of the chain. 
				return int_cons,False

		return int_cons,True	



#Finds the neighbours of a constituency by adding ward shapes
def Rneighbours(self,condf):
	df = self.df
	Neighbours = []
	for i in range(len(condf['con'])):
		wardsi = condf['wards'][i]
		neigh_i = []
		for j in range(len(condf['con'])):
			if i != j:
				wardsj = condf['wards'][j]
				for wardi,wardj in itertools.product(wardsi,wardsj):
						if wardj in df['neighbours'][wardi]:
							neigh_i.append(j)
							break
		Neighbours.append(neigh_i)
	condf['neighbours'] = Neighbours
	return condf	









#### FlipSwap algorithm
######################

class FlipSwap:

	Method = 'MCMC'

	def __init__(self,df,condf,pop_constraint,comp_constraint,alg_type,comp=None,relabel=None,OAframe=None): #df = main ward data, OAframe = OA area for MMI compactness, comp = area or MMI depending on compactness measure, pop and comp constraints are numeric bounds for sampling, alg_type = F, FS or R for flip, flip and swap or random sampler.
		self.df = df
		self.condf = condf
		if comp == 'MMI':
			self.OAframe = OAframe
			self.comp = MC.MMIcompact
		else:
			self.comp = MC.Areacompact
		if relabel == 'bynumber':
			self.relabel = MC.relabel_bynumber
		else: 
			self.relabel = MC.relabel_bylocation
		self.pop_constraint = pop_constraint
		self.comp_constraint = comp_constraint
		self.alg_type = alg_type

		number_cons = len(set(df['constituency']))
		self.avpop = sum(df['All Ages'])/number_cons

	def update(self,df,condf,edgepairs,mlep,b1,b2,E2=None,Gibbs=None): #Tune b1 for map movement (tempering) and b2 for sampling (only compactness)
	#Updates the chain each time. df: ward data, OAframe: ordiance survey data for compactness, edgepairs: list of proposals, mlep: updated weighting of maximum options from map (needed to normalise sampler)
	#self.avpop: average population over all constituencies, samp: indicator if population constraint is met, b: tempering coefficient, popconstraint:percentage population constraint, Gibbs:value of Gibbs for the previous map (used in case map is not changed)
	#Outout = df,edgepairs: general info. samp: whether to keep map in list, Gibbs: Gibbs for new map. 

		if Gibbs == None: #For the initial map only
			condf = self.comp(self,df,condf,range(len(condf.index)))
			E1 = np.mean(condf['comp_score'])
			condf = MC.population(self,condf,self.avpop,range(len(condf.index)))
			E2 = np.mean(condf['pop_score'])

			Gibbs = np.exp(-b1*E1-b2*E2)

		##Flip Swap or Random
		if self.alg_type == 'F':
			op1 = 0
			op2 = Gibbs
		elif self.alg_type == 'FS':
			op1 = np.exp(-4.2*E2)
			op2 = Gibbs
		elif self.alg_type == 'R':
			op1 = 0
			op2 = 1

		if random.random() < 0.5:
			return(df,edgepairs,Gibbs,E2,condf)

		if random.random() > len(edgepairs)/mlep: #Weight it so that a uniform dist is produced (maps with less exit points don't create skew)
			return(df,edgepairs,Gibbs,E2,condf)


		#print("hello pussy",E2,np.exp(-4.2*E2)) #4 is the best so far my experiment. 10 is really bad ~ 0 acceptance. 6 >2% acc, 5 4% acc, 4 10% acc, 3 3%
		if random.random() < op1: #make this 0 for OF alg, np.exp(-4.2*E2) for OFS 			#swap if population deviance is low
			while True:
				y = random.sample(edgepairs,1)
				ward,con = y[0]
				othercon = df.at[ward,'constituency']
				edgepairscop,dfcop,condfcop,E1cop,E2cop = flip(self,edgepairs,df,condf,ward,con,b1,b2)
				edgepairs2 = [edge for edge in edgepairscop if edge[0] in condf['wards'][con-1] and edge[1] == othercon and edge[0] != ward]
				if len(edgepairs2) != 0: #Run again if swap doesn't exist. 
					break
			y2 = random.sample(edgepairs2,1)
			ward,con = y2[0]
			edgepairscop,dfcop,condfcop,E1cop,E2cop = flip(self,edgepairscop,dfcop,condfcop,ward,con,b1,b2)
			Gibbscop = np.exp(-b1*E1cop-b2*E2cop)
			if random.random() > Gibbs:
				return df,edgepairs,Gibbs,E2,condf
		else:
			y = random.sample(edgepairs,1)
			ward,con = y[0]
			edgepairscop,dfcop,condfcop,E1cop,E2cop = flip(self,edgepairs,df,condf,ward,con,b1,b2) #Got an internal acceptance
			Gibbscop = np.exp(-b1*E1cop-b2*E2cop)
			if random.random() > op2: #1 if you want a random sampler - also need 0 at the top (line 242)
				return df,edgepairs,Gibbs,E2,condf

		return dfcop,edgepairscop,Gibbscop,E2cop,condfcop


	def Run(self,iterations,b1,b2):
		#Runs the whole aglorithm, keeping appropriate samples and recording their weightings for resampling. df,OAframe: information (see update), iterations: time to run, b:tempering coeff, popconstraint: percentage population constraint.
		#Output = all maps and their Gibbs values
		df = self.df.copy()
		condf = self.condf.copy()

		#Initialise edge wards
		edgepairs = set()
		for i in range(len(df['wd18cd'])):
			for n in df['neighbours'][i]:
				proposed_con = [x for x in df.loc[df['constituency'] == df['constituency'][i]].index.tolist() if x != i]
				if df['constituency'][n] != df['constituency'][i] and MC.connected(df,proposed_con) and len(df.loc[df['constituency'] == df['constituency'][i]].index.values) > 1: #Do MC.connected at this stage
					edgepairs.add((i,df['constituency'][n]))


		mlep = len(edgepairs) #Find maximum number of exit points for a map. - updates this number to approximate theoretical maximum across whole space

		t = 0
		accept = 0
		initial_con = df['constituency']
		df,edgepairs,Gibbs,E2,condf = self.update(df,condf,edgepairs,mlep,b1,b2) #First time is different because there is no Gibbs
		dist = []
		while t <= iterations:

			if len(edgepairs) > mlep: #Update mlep so it is always the biggest recorded
				mlep = len(edgepairs)
			df,edgepairs,Gibbs,E2,condf = self.update(df,condf,edgepairs,mlep,b1,b2,E2,Gibbs)
			vec = [0]*len(df['constituency'])
			for i in range(3):
				for x in condf.at[i,'wards']: #Not all wards present
					vec[x] = i+1
			vec = self.relabel(dist,vec) #Do a reassignment of constituencies


			if (any(abs(x)>self.pop_constraint for x in condf['pop_score']) == False) and (any(x<self.comp_constraint for x in condf['comp_score']) == False): #Only when desired condition (population) is satisfied
				accept +=1 
				#put into format, each index is a ward that corresponds to its index in df and the value is a number of a constituency
				vec = [0]*len(df['constituency'])
				for i in range(3):
					for x in condf.at[i,'wards']: #Not all wards present
						vec[x] = i+1
				vec = self.relabel(dist,vec) #Do a reassignment of constituencies
				dist.append([vec,Gibbs])
			t += 1

		return dist



def flip(self,edgepairs,df,condf,ward,con,b1,b2):
	condfcop = condf.copy()
	dfcop = df.copy()
	edgepairscop = edgepairs.copy()

	#Add ward to constituency
	ind_c1 = con-1
	condfcop.at[ind_c1,'wards'] = [x for x in condfcop.at[ind_c1,'wards']] +[ward]

	condfcop.at[ind_c1,'geometry'] = gpd.GeoSeries(dfcop['geometry'][condfcop.at[ind_c1,'wards']]).unary_union

	#Remove ward from constituency
	othercon = dfcop.at[ward,'constituency']
	ind_c = othercon - 1
	condfcop.at[ind_c,'wards'] = [x for x in condfcop.at[ind_c,'wards'] if x != ward]
	condfcop.at[ind_c,'geometry'] = gpd.GeoSeries(dfcop['geometry'][condfcop.at[ind_c,'wards']]).unary_union


	#Find measures of new swap
	ind = [ind_c1,ind_c]
	condfcop = self.comp(self,dfcop,condfcop,ind)
	E1cop = np.mean(condfcop['comp_score'])
	condfcop = MC.population(self,condfcop,self.avpop,ind)
	E2cop = np.mean(condfcop['pop_score'])


	#Update edgepairs
	dfcop.at[ward,'constituency'] = con
	for edge in edgepairscop.copy():
		if edge[0] in set([x for x in dfcop.loc[dfcop['constituency'].isin([othercon,con])].index.values] + [x for x in dfcop.loc[ward,'neighbours']]) or edge[1] in [othercon,con]:
			edgepairscop.remove(edge)


	both_con_wards = [x for x in dfcop.loc[dfcop['constituency'].isin([othercon,con])].index.values]
	for n in both_con_wards + list(set([x for xs in dfcop.loc[both_con_wards,'neighbours'].values for x in xs])):
		for m in dfcop['neighbours'][n]:
			proposed_con = [x for x in dfcop.loc[dfcop['constituency'] == dfcop['constituency'][n]].index.tolist() if x != n]
			if dfcop['constituency'][n] != dfcop['constituency'][m] and MC.connected(dfcop,proposed_con) and len(dfcop.loc[dfcop['constituency'] == dfcop['constituency'][n]].index.values) > 1:
				edgepairscop.add((n,dfcop['constituency'][m]))

	for edge in edgepairscop.copy():
		if len(dfcop.loc[dfcop['constituency'] == dfcop['constituency'][edge[0]]].index.values) <= 1 or edge[0] in dfcop.loc[dfcop['constituency'] == edge[1]].index.values:
			edgepairscop.remove(edge)

	return edgepairscop,dfcop,condfcop,E1cop,E2cop




##Swendsen-Wang Algorithm
#########################


class SwendsenWang:
	method = 'MCMC'

	def __init__(self,df,condf,pop_constraint,comp_constraint,comp=None,relabel=None,OAframe=None): #df = main ward data, OAframe = OA area for MMI compactness, comp = area or MMI depending on compactness measure, pop and comp constraints are numeric bounds for sampling, alg_type = F, FS or R for flip, flip and swap or random sampler.
		self.df = df
		self.condf = condf
		if comp == 'MMI':
			self.comp = MC.MMIcompact
			self.OAframe = OAframe
		else:
			self.comp = MC.Areacompact
		if relabel == 'bynumber':
			self.relabel = MC.relabel_bynumber
		else: 
			self.relabel = MC.relabel_bylocation
		self.pop_constraint = pop_constraint
		self.comp_constraint = comp_constraint

		number_cons = len(set(df['constituency']))
		self.avpop = sum(df['All Ages'])/number_cons
	


	def Run(self,iterations,q=None,b=None,b1=None,b2=None):
		df = self.df.copy()
		condf = self.condf.copy()
		condf['con'] = [x+1 for x in condf['con']]
		default_params = [0.04,20,1,5]
		params = [q,b,b1,b2]
		self.params = [x if x != None else default_params[i] for x,i in enumerate(params)]

		edgeframe = updateedgeframe(df)

		t = 0 
		dist = []
		while t <= iterations:
			Gibbs = 0.2
			df,condf,Gibbs = UpdateS(self,df,edgeframe,condf,Gibbs)
			t +=1

			if (any(abs(x)>self.pop_constraint for x in condf['pop_score']) == False) and (any(x<self.comp_constraint for x in condf['comp_score']) == False):
				vec = []
				cons = list(set(df['constituency'].values.tolist()))
				for x in range(len(df['constituency'])):
					vec.append(cons.index(df['constituency'][x]))
				vec = [x+1 for x in vec]
				dist.append([vec,Gibbs])
		return dist 




def UpdateS(self,df,edgeframe,condf,Gibbs):
	q,b,b1,b2 = self.params
	start_cycle = time.time()
	tryy = 0
	while True:
		edgeframe = updateedgeframe(df,q,edgeframe)
		Allblocks = blocks(df,edgeframe)
		bblocks = boundaryblocks(Allblocks,df)

		lim = len(bblocks)
		while True:
			while True: #Sample from trunc poisson 
				R = sct.poisson.rvs(1,size=1)[0]
				if R != 0 and R <= lim:
					break

			ans = "Calm" #For connectivity at the bottom
			r = 0
			bblockspick = bblocks.copy()
			selected_blocks = []
			while r <= R and len(bblockspick) >= 1:
				ans = "Calm"
				blocks_wards = [x for bblockspick in selected_blocks for x in bblockspick]
				chosenblock = random.sample(bblockspick,1)
				chosenblock = chosenblock[0]
				bblockspick.remove(chosenblock)
				for ward in chosenblock:
					if ans == "Error":
						break
					for neigh in df['neighbours'][ward]:
						if neigh in blocks_wards:
							ans = "Error"
							break

				con = condf.at[condf.loc[condf['con'] == df['constituency'][chosenblock[0]]].index[0],'wards'] #Constituency that block is part of
				if ans != "Error" and MC.connected(df,[x for x in con if x not in chosenblock]): #Initial check. 
					interiors = []
					for interior in gpd.GeoSeries(df['geometry'][[x for x in con if x not in chosenblock]]).unary_union.interiors:
						interiors.append(list(interior.coords))
					if interiors == []:
						selected_blocks.append(chosenblock)
						r +=1
			if r >= R:
				cont = True
				break
			if tryy >= 7:
				cont = False
				break
			tryy += 1
		if cont:
			break


	#Proposal
	dfcop = df.copy()
	condfcop = condf.copy()
	changed_cons = []
	for block in selected_blocks:
		options = set()
		for ward in block:
			for neigh in dfcop['neighbours'][ward]:
				if dfcop['constituency'][ward] != dfcop['constituency'][neigh]:
					options.add(dfcop['constituency'][neigh])
		choice = random.sample(options,1)
		changed_cons = changed_cons + dfcop['constituency'][block].tolist() + choice
		dfcop.loc[block,'constituency'] = choice

	for con in set(changed_cons):
		ind_w = dfcop.loc[dfcop['constituency'] == con].index.values.tolist()
		condfcop.at[condfcop.loc[condfcop['con'] == con].index[0],'wards'] = ind_w
	ind = condfcop.loc[condfcop['con'].isin(changed_cons)].index.values

	#Accepting
	oldprop = len(bblocks)
	bblocks2 = boundaryblocks(Allblocks,dfcop)
	newprop = len(bblocks2)


	f_oldprop = sct.poisson.cdf(oldprop,1)
	f_newprop = sct.poisson.cdf(newprop,1)


	boundaryedge = [0,0]
	for i in range(2): #Counts Gibbs of edges that cross boundary for each map
		mydf = [df,dfcop][i]
		for edge in range(len(edgeframe['node1'])):
			if mydf['constituency'][edgeframe['node1'][edge]] != mydf['constituency'][edgeframe['node2'][edge]]:
				boundaryedge[i] +=1


	ind = list(ind)
	condfcop = self.comp(self,dfcop,condfcop,ind)
	condfcop = MC.population(self,condfcop,self.avpop,ind)

	E1new = np.mean(condfcop['comp_score'])
	E2new = np.mean([abs(x) for x in condfcop['pop_score']])
	E1old = np.mean(condf['comp_score'])
	E2old = np.mean([abs(x) for x in condf['pop_score']])


	#Accept with this prob
	Gibbscop = min(1,(1-q)**(boundaryedge[1]-boundaryedge[0])*(oldprop/newprop)**(R)*(f_oldprop/f_newprop)*np.exp(-b*(b2*E2new+b1*E1new)+b*(b2*E2old+b1*E1old)))

	if random.random() < Gibbscop:
		return dfcop,condfcop,Gibbscop
	

	return df,condf,Gibbs


#Running
def updateedgeframe(df,q=None,edgeframe=None):
	if edgeframe is None:
		edgelist1 = []
		edgelist2 = []
		edgelist = []
		for ward in range(len(df['neighbours'])):
			for neigh in df['neighbours'][ward]:
				if [neigh,ward] not in edgelist:
					edgelist1.append(ward)
					edgelist2.append(neigh)
					edgelist.append([ward,neigh])
		edgeframe = pd.DataFrame({'node1':edgelist1,'node2':edgelist2,'on':range(len(edgelist))})

	else:
		for edge in range(len(edgeframe['node1'])):
			ward1,ward2 = [edgeframe['node1'][edge],edgeframe['node2'][edge]]
			if random.random() <= q and df['constituency'][ward1] == df['constituency'][ward2]:
				edgeframe['on'][edge] = 1
			else:
				edgeframe['on'][edge] = 0
	return(edgeframe)


def boundary(df):
	boundary = set()
	for i in range(len(df['wd18cd'])):
		for n in df['neighbours'][i]:
			if df['constituency'][n] != df['constituency'][i]:
				boundary.add(i)
	return(boundary)




def blocksold(edgeframe):
	block = []
	for ward in range(len(df['constituency'])):
		wardsdone = []
		block_ward = set()

		while len(wardsdone) != len(block_ward):
			current_ward = wardsleft[0]
			block_ward.add(current_ward)
			for edge in range(len(edgeframe['node2'])):
				if edgeframe['node1'][edge] == current_ward and edgeframe['on'][edge] == 1:
					block_ward.add(edgeframe['node2'][edge])
			wardsleft = [x for x in wardsleft if x != current_ward]

		block.append(list(block_ward))
	return(block)


def blocks(df,edgeframe):
	block = []
	uniwardsdone = []
	for ward in range(len(df['constituency'])):
		if ward not in uniwardsdone:
			wardsdone = []
			block_ward = set()
			block_ward.add(ward)
			while len(wardsdone) != len(block_ward):
				current_ward = [x for x in block_ward if x not in wardsdone][0]
				for edge in range(len(edgeframe['node2'])):
					if edgeframe['node1'][edge] == current_ward and edgeframe['on'][edge] == 1:
						block_ward.add(edgeframe['node2'][edge])
					elif edgeframe['node2'][edge] == current_ward and edgeframe['on'][edge] == 1:
						block_ward.add(edgeframe['node1'][edge])
				wardsdone.append(current_ward)
			block.append(list(block_ward))
			uniwardsdone = uniwardsdone + wardsdone
	return(block)



def boundaryblocks(Allblocks,df):
	bblocks = []
	boundaryy = boundary(df)
	for block in Allblocks:
		for ward in block:
			if ward in boundaryy:
				bblocks.append(block)
				break
	return(bblocks)


