'''2018-09-12
       Adapt EFAProdModel to implement Graham Shaeffer model.
'''

from pymc import  *
from pylab import *
from tables import *
from numpy import log
import tables
import warnings
from datetime import datetime as dt

from EFAProdModel import EFAProdModel as oldEFAProdModel




class EFAProdModel(oldEFAProdModel):
    def __init__(self,SiteNumber,Sites,VirginSites,CoastLength,oldB=1):
        '''EFAProdModel(SiteNumber,Sites,VirginSites,CoastLength,oldB=1)
        * Sites and VirginSites are lists of instances of the Site and VirginSite classes
        * CoastLength is a dictionary giving lengths in metres for each site. e.g.{0:10000,2:10002,4:10004,8:1008,16:10016}            
        * oldB is the original relative-abundance.   Value is probably one
        
         '''
        self.SiteNumber=SiteNumber
        self.nsite=len(self.SiteNumber)
        self.Sites=Sites
        self.VirginSites=VirginSites
        self.CoastLength=CoastLength
        self.oldB=oldB


        self.ProductivityParameters()  #Generate nodes for productivity function       
        self.GrandMeans()              #Generate nodes for Grand Mean and virgin Grand Mean of biomass density
        self.EffectStandardDeviations()#Generate nodes for standard deviations of effects
        self.SiteEffects()             #                   site effects
        self.TransectEffects()         #                   transect effects
        self.YearEffects()             #                   Year Effects        
        self.SiteAreas()               #                   Site Areas
        self.VirginBiomass()           #                   Virgin biomass
        self.CalcRelAbund()            #                   Relative Abundance
                              
        # Mean Weight nodes will be created as they are needed
        self.MeanWeight =[]
             
        self.VirginTransects()#Nodes to represent expected and observed number of animals for initial survey (untrimmed transects)
        self.SurveyTransects()#Nodes to represent expected and observed number of animals for trimmed transects
        self.SiteAreasFromData()#Nodes to represent site-areas as sampled directly from statistics
        

        print('Finished Constructing Model')
    


        
    def ProductivityParameters(self):
        #parameters of Graham-Shaeffer productivity function  
        self.fmax=Uniform('fmax',.01,0.25)        
        self.xmax=Lambda('xmax',lambda dummy=self.fmax:0.5)        
        self.a=Lambda('a',lambda dummy=self.fmax:1.0) 
        self.b=Lambda('b',lambda dummy=self.fmax:1.0)
            

if __name__ == "__main__":
  import pickle
  pfile="../MCMC/Jervis/Jervis.pickle"
  SiteNumber,Sites,VirginSites,CoastLength=pickle.load(open(pfile,"rb"))  
  
  ConstantValues.MinYear=1999  #Earliest year to be considered in calculations
  
  test=EFAProdModel(SiteNumber,Sites,VirginSites,CoastLength,oldB=1)
  #print(test.GetMeanWeight(16,[2005,2,14]))  

