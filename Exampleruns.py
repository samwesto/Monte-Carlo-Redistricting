import MCMC_functions as MC
import MC_algorithms as MCA


#Input files
Dataframe = 'Data/Reference/WMdfREF.csv'
Shapefile = 'Data/Reference/WMshapefileREF.shp'
Adjacencylist = 'Data/Reference/WMadjlistREF.csv'
OAframedata = 'Data/Compactness/OAframe.shp'
OAcol = 'Data/Compactness/postward.csv'



#Swendsen Wang test
starting_map = [1,2,2,1,1,2,3,2,2,2,3,3,3,3,1,3,1,3,3,2,2,1]
df,condf,OAframe = MC.Initialise(Dataframe,Shapefile,Adjacencylist,OAframedata,starting_map,OAcol)


SWang = MCA.SwendsenWang(df,condf,0.4,0.05,OAframe=OAframe)

distribution = SWang.Run(100)
print(distribution)

distance = MC.LWasserstein(df,distribution,0,1)
print(distance)


#Flip Swap test
MainFlipChain = MCA.FlipSwap(df,condf,0.4,0.15,'FS','area',OAframe=OAframe)

distribution = MainFlipChain.Run(100,1,4)
print(distribution)

distance = MC.LWasserstein(df,distribution,0,1)
print(distance)

#Rodden test
starting_map = range(22)
df,condf,OAframe = MC.Initialise(Dataframe,Shapefile,Adjacencylist,OAframedata,starting_map,OAcol)


Rodden = MCA.RoddenChain(df,condf,0.2,0.15,'area',OAframe=OAframe)

distR = Rodden.Run(3,3)
print(distR)
distanceR = MC.LWasserstein(df,distR,0,1)
print(distanceR)





