import MCMC_functions as MC
import MC_algorithms as MCA


#Input files
Dataframe = 'Data/Reference/WMdfREF.csv'
Shapefile = 'Data/Reference/WMshapefileREF.shp'
Adjacencylist = 'Data/Reference/WMadjlistREF.csv'
OAframedata = 'Metrics/Compactness/OAframe.shp'
OAcol = 'Metrics/Compactness/postward.csv'

starting_map = range(22)
df,condf,OAframe = MC.Initialise(Dataframe,Shapefile,Adjacencylist,OAframedata,starting_map,OAcol)

MainRchain = MCA.RoddenChain(df,condf,OAframe,'area',0.2,0.15)

print(MainRchain.Run(10,3))