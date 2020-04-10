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

class TimeError(Exception):
	pass

class RoddenChain:
	Method = 'constructive'

	def __init__(self,df,condf,OAframe,comp,pop_constraint,comp_constraint): #df = main ward data, OAframe = OA area for MMI compactness, comp = area or MMI depending on compactness measure, pop and comp constraints are numeric bounds for sampling. 
		self.df = df
		self.OAframe = OAframe
		self.condf = condf
		self.comp = comp
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

				if self.comp == 'area': #Use compactness measure determined by user
					 int_cons = MC.Areacompact(self,int_cons,range(len(int_cons)))
				else:
					int_cons = MC.MMIcompact(self,int_cons,range(len(int_cons)))

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
			dist.append(MC.relabel_bynumber(vec))
		return dist



#Creates a new independent map each time
def RunRalg(self,con_number):
	int_cons = self.condf.copy()
	int_df = self.df
	OAframe = self.OAframe

	avpop = sum(int_df['All Ages'])/con_number

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


	#Second part of Rodden (population evening) - creating maps that meet population int_constraint
	ind = list(range(len(int_cons['pop_score'])))
	int_cons = MC.population(self,int_cons,avpop,ind)

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
		int_cons = MC.population(self,int_cons,avpop,ind)
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





