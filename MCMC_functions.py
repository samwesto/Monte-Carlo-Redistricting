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

	starting_map = relabel_bynumber(starting_map) #Make numeric


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


	cons = gpd.GeoDataFrame({'con':starting_map,'geometry':list(df['geometry']),'neighbours':list(df['neighbours']), 'wards':[[x] for x in range(len(df['wd18cd']))],'pop_score':range(len(df['wd18cd'])),'comp_score':range(len(df['wd18cd']))},geometry='geometry')

	return df,cons,OAframe




def relabel_bynumber(mapi): #assgins constituency name to the constituency that best matches the original layout. Done via linear assignment
	order = sorted(set(mapi), key=mapi.index)
	vec = [order.index(x)+1 for x in mapi]
	return vec




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


def Areacompact(self,condf,ind=None):
	df = self.df
	if ind == None:
		return condf
	for con,i in list(zip(condf['geometry'][ind],ind)):
		condf.at[i,'comp_score'] = (4*np.pi*con.area)/(con.length**2)
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


