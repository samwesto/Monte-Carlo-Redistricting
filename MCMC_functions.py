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
import cvxpy as cp
import scipy.optimize as opt
from itertools import product


def Initialise(Dataframe,Shapefile,Adjacencylist,OAframedata,starting_map,OAcol = False): #Builds the two data frames from ready made files.
	
	df = pd.read_csv(Dataframe)
	shape = gpd.read_file(Shapefile)

	#Check for correct data
	if any(type(x) != shapely.geometry.polygon.Polygon for x in shape['geometry']):
		raise TypeError('Shape was not a shapely polygon')

	if any(type(x) not in [str,int] for x in df.iloc[:,1]):
		raise TypeError('Population was not a string or integer')
		if any(type(x) != int for x in df.iloc[:,1]):
			df.iloc[:,1] = df.iloc[:,1].astype(int)

	if any(type(x) not in [str,int,] for x in df.iloc[:,0]):
		raise TypeError('ward codes was not a string or int')

	starting_map = relabel_bynumber([],starting_map) #Make numeric


	#Add adjlist to df
	with open(Adjacencylist) as inputfile:
	    reader = csv.reader(inputfile)
	    inputm = list(reader)
	df['neighbours'] = inputm

	for i in range(len(df['neighbours'])):
		for j in range(len(df['neighbours'][i])):
			df['neighbours'][i][j] = int(df['neighbours'][i][j])
	


	##FOR MY USE ONLY
	#OAframe
	if OAcol != False:
		ref = [336,339,340,341,342,344,349,350,357,360,364,365,371,372,373,375,377,378,380,384,385,387]
		#ref = [707,702, 717, 706, 710, 714, 716] #For micro dist
		OAs = []
		with open(OAcol) as inputfile:
		    reader = csv.reader(inputfile)
		    inputm = list(reader)
		for i in ref:
			OAs.append(inputm[i])



	OAframe = gpd.read_file(OAframedata)

	#Add geometry to df
	if OAs == False:
		df = gpd.GeoDataFrame({'wd18cd':df.iloc[:,0],'All Ages':df.iloc[:,1],'constituency':starting_map,'neighbours':df['neighbours'],'geometry':shape['geometry']},geometry='geometry')
	else:
		df = gpd.GeoDataFrame({'wd18cd':df.iloc[:,0],'All Ages':df.iloc[:,1],'constituency':starting_map,'neighbours':df['neighbours'],'geometry':shape['geometry'],'OAs':OAs},geometry='geometry') #Change to df['O']


	congeom = []
	for x in range(len(set(starting_map))):
		x = x+1
		congeom.append(gpd.GeoSeries(df['geometry'][df.loc[df['constituency'] == x].index.values]).unary_union)
	if len(set(starting_map)) == len(starting_map):
		neighs = list(df['neighbours'])
		wards = [[y] for y in range(len(df['wd18cd']))]
	else: 
		neighs = range(len(set(starting_map)))
		wards = []
		for x in range(len(set(starting_map))):
			x = x+1
			wards.append(df.loc[df['constituency'] == x].index.values)

	cons = gpd.GeoDataFrame({'con':range(len(set(starting_map))),'geometry':congeom,'neighbours':neighs, 'wards':wards,'pop_score':range(len(set(starting_map))),'comp_score':range(len(set(starting_map)))},geometry='geometry')

	return df,cons,OAframe




####Sorting and distributing
############################

def relabel_bynumber(dist,mapi): #assgins constituency name to the constituency that best matches the original layout. Done via linear assignment
	order = sorted(set(mapi), key=mapi.index)
	vec = [order.index(x)+1 for x in mapi]
	return vec

def relabel_bylocation(dist,mapi): #assgins constituency name to the constituency that best matches the original layout. Done via linear assignment
	if dist == []:
		return mapi
	map1 = dist[0][0]
	matrix = np.full(shape = [len(set(map1)),len(set(map1))], fill_value = +1000)

	for ward in range(len(map1)):
		matrix[map1[ward]-1][mapi[ward]-1] -= 1 

	row, col = opt.linear_sum_assignment(matrix)

	
	return [list(col).index(x-1)+1 for x in mapi]


def impsample(dist,size):
	values = [1/float(x[1]) for x in dist] #Requires dist to be in [vec,Gibbs] format
	a = [x[0] for x in dist]
	denom = sum(values)
	return random.choices(a, k=size , weights= [1/(float(x[1])*denom) for x in dist]) #sample from maps with weight 1/Gibbs to make it uniform



###Metrics
##########

def MMIcompact(self,condf,ind = None): 
	#Uses moments of inertia, weighted by population
	#Out = metric of comparison to circle with same area and same centroid. 
	OAframe = self.OAframe
	df = self.df
	if ind == None:
		return condf
	i = 0
	for con,ws in list(zip(condf['geometry'][ind],condf['wards'][ind])):
		center = con.centroid
		dist = []
		circdist = []
		OAs = []
		distances = []
		for ward in ws:
			OAs = OAs + df['OAs'][ward]
		for OA in OAs:
			index = OAframe.loc[OAframe['OA Code'] == OA].index.values[0]
			pop = OAframe['pop'][index]
			point = OAframe['geometry'][index]
			distances.append(point.distance(center))
			dist.append((point.distance(center))**2*float(pop)) 
		circle = center.buffer(max(distances)) #Creates smallest encompassing circle
		OAs = []
		#Search for OAs in encompassing circle algorithm
		wardsdone = []
		for ward in ws:
			OAs = OAs + df['OAs'][ward]
			wardsdone.append(ward)
		while True:
			stillgoing = []
			for ward in ws:
				for neigh in df['neighbours'][ward]:
					if neigh not in wardsdone:
						if df['geometry'][neigh].intersection(circle).is_empty == False:
							OAs = OAs + df['OAs'][neigh]
							stillgoing.append(neigh)
						wardsdone.append(neigh)
			ws = stillgoing
			if stillgoing == []:
				break
		#for ward in range(len(df['constituency'])):
		#	if df['geometry'][ward].intersection(circle).is_empty == False:
		#		OAs = OAs + df['OAs'][ward]
		for OA in OAs:
			index = OAframe.loc[OAframe['OA Code'] == OA].index.values[0]
			pop = OAframe['pop'][index]
			point = OAframe['geometry'][index]
			if OAframe['geometry'][index].within(circle):
				circdist.append((point.distance(center))**2*float(pop)) #Add in the population weighting here...
		realarea = sum(dist)
		circarea = sum(circdist)
		condf.loc[ind[i],'comp_score'] = realarea/circarea
		i += 1
	return condf


def Areacompact(self,df,condf,ind=None):
	if ind == None:
		return condf
	for i in ind:
		con = df.loc[condf.loc[i,'wards'],['geometry']].unary_union
		condf.loc[i,'comp_score'] = (4*np.pi*con.area)/(con.length**2)
	return condf



def population(self,condf,avpop,ind = None): #returns non abs value
	#measures the population deviance of all constituencies under a given map. avpop: average pop over all constituencies, con: the constituency gaining a ward, ward: the ward moving to new constituency
	#out = list of populations by constituency (averaged later to make E2)
	df = self.df
	if len(ind) == 0:
		return condf
	for i in ind:
		condf.loc[i,'pop_score'] = (sum(df['All Ages'][condf.at[i,'wards']]) - avpop)/avpop
	return condf




####Wasserstein distance
########################
def Vmatrix(refdist,i,j):
	if len(refdist) < max(i,j)+1:
		raise ValueError('Requires ' +str(max(i,j)+1) + ' or more maps to calculate this distance')
	Allparts = []
	refdist = [x[0] for x in refdist]
	for redist in refdist:
		V = np.array([range(len(redist))]).transpose() #Creates dummy first column to append to
		for y in range(1,len(set(redist))+1): #For each district, shift all numbers by +1 due to R indexing start at 1
			v = np.array([[1 if ele == y else 0 for ele in redist]]) #Create district vector where value is an indicator if ward is in that discrict: vi, i,..,k
			V = np.concatenate([V,np.divide(v,sum(v[0])).transpose()],axis=1) #Normalise and append vector to V
		V = V[:,1:] #Remove first column from matrix V
		Allparts.append(V) 
	return Allparts


#Indicdent matrix
def Imatrix(df):
	#Generate edge list
	edgelist = []
	for ward in range(len(df['neighbours'])):
		for neigh in df['neighbours'][ward]:
			if [neigh,ward] not in edgelist:
					edgelist.append([ward,int(neigh)])

	in_matrix = np.array([range(len(edgelist))])
	for x in range(len(df.index)):
		in_row = np.array([[1 if x in edge else 0 for edge in edgelist]])
		in_matrix = np.concatenate([in_matrix,in_row],axis=0)
	in_matrix = in_matrix[1:,:]

	return in_matrix

def distance(a_indicator, b_indicator,edge_incidence): #elements taken from Github - https://github.com/vrdi/optimal-transport/blob/master/wasserplan.py
	n_edges = edge_incidence.shape[1]
	edge_weights = cp.Variable(n_edges)
	diff = b_indicator - a_indicator
	objective = cp.Minimize(cp.sum(cp.abs(edge_weights)))
	conservation = (edge_incidence @ edge_weights) == diff
	prob = cp.Problem(objective, [conservation])
	prob.solve(solver='ECOS') 
	return np.sum(np.abs(edge_weights.value))


def pairwise_d(no_districts,Allparts,P,district1,district2):
	dist = np.zeros((no_districts,no_districts))
	for ina in range(no_districts):
		for inb in range(no_districts):
			dis = distance(Allparts[district1][:,ina],Allparts[district2][:,inb],P)
			dist[ina][inb] = dis
	return dist


def partdistance(no_districts,Allparts,P,district1,district2):
	distances = pairwise_d(no_districts,Allparts,P,district1,district2)
	inda, indb = opt.linear_sum_assignment(distances)

	chosenones  = {a_in: b_in 
	for a_in, b_in in zip(inda,indb)}

	total_dist = 0
	for a_in, b_in in chosenones.items():
		total_dist += distances[a_in][b_in]

	return total_dist

def LWasserstein(adjlist,refdist,i,j): #Recommend using sklearn.manifold.MDS for embedding into euclidean space. 
	Allparts = Vmatrix(refdist,i,j)
	P = Imatrix(adjlist)

	no_districts = Allparts[0].shape[1]
	return partdistance(no_districts,Allparts,P,i,j)





#Ward relationships
###################


def neighbours(Shapefile): #Takes in shapefile and creates adjacency list for wards, saving it in the Data folder. 
	shape = gpd.read_file(Shapefile)

	neighbours = []
	for i in range(len(shape['geometry'])):
		if type(shape['geometry'][i]) != shapely.geometry.polygon.Polygon:
			print(type(shape['geometry'][i]))
		neighbours_i = []
		for j in range(len(shape['geometry'])):
			if type(gpd.GeoSeries(shape['geometry'][[i,j]]).unary_union) == shapely.geometry.collection.GeometryCollection:
				print(type(gpd.GeoSeries(shape['geometry'][[i,j]]).unary_union))
			if shape['geometry'][i].touches(shape['geometry'][j]) and type(gpd.GeoSeries(shape['geometry'][[i,j]]).unary_union) != shapely.geometry.multipolygon.MultiPolygon: #Connect at more than one point
				neighbours_i.append(j)
		neighbours.append(neighbours_i)

	with open("Data/WMadjlist.csv", "w", newline="") as f:
	    writer = csv.writer(f)
	    for row in neighbours:
	        writer.writerow(row)

def connected(df,con,start_vertex=None, wardsmet= None): 
	if len(con) == 0:
		return False
	if wardsmet is None: 
		wardsmet=set()
	if start_vertex is None:
		start_vertex = random.choice(con) #Random start vertex if not provided
	wardsmet.add(start_vertex)
	if len(wardsmet) != len(con):
		for h in df['neighbours'][start_vertex]:
			if h in con and h not in wardsmet:
				if connected(df,wardsmet=wardsmet,start_vertex=h,con=con):
					return(True)
	else:
		return True
	return(False)



