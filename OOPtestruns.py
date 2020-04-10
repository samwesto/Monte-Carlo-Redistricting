import MCMC_functions as MC
import MC_algorithms as MCA


#Input files
Dataframe = 'Data/Reference/WMdfREF.csv'
Shapefile = 'Data/Reference/WMshapefileREF.shp'
Adjacencylist = 'Data/Reference/WMadjlistREF.csv'
OAframedata = 'Metrics/Compactness/OAframe.shp'
OAcol = 'Metrics/Compactness/postward.csv'

starting_map = [1,2,2,1,1,2,3,2,2,2,3,3,3,3,1,3,1,3,3,2,2,1]
df,condf,OAframe = MC.Initialise(Dataframe,Shapefile,Adjacencylist,OAframedata,starting_map,OAcol)


MainFlipChain = MCA.FlipSwap(df,condf,OAframe,'area',0.2,0.15,'FS')

print(MainFlipChain.Run(10,0.2,0.4))


starting_map = range(22)
df,condf,OAframe = MC.Initialise(Dataframe,Shapefile,Adjacencylist,OAframedata,starting_map,OAcol)


Rodden = MCA.RoddenChain(df,condf,OAframe,'area',0.2,0.15)

print(Rodden.Run(10,3))