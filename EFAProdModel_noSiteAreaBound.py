'''
2019-0619
    Special version where Site-area is not restricted to +/- 2 standard Areas
'''

from pymc import  *
from pylab import *
from tables import *
from numpy import log
import tables
import warnings
from datetime import datetime as dt

from ADO import adoBaseClass as daoBaseClass
from constants import ConstantValues
from EFAsite import  Site,VirginSite
from Harvest import  Harvest
from LastEffect import LastEffect
from MakeName import MakeName
from MeanWeight import MeanWeight
from NextRelBmass import CurRelBmass
from NodeNameToDecimalYear import MatchNodeRefTime
from ToDecimalYear import ToDecimalYear,FromDecimalYear
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
        

        print('Finished Constructing Model')
    

    def SiteAreas(self): 
        #Site areas
        nsite=len(self.SiteNumber)
        w=[ self.CoastLength[s]  for s in self.SiteNumber]
        l=[ t.AvgTranLen  for t in self.VirginSites]
        se=[ t.sterrTranLen  for t in self.VirginSites]
        tau=[1/t/t for t in se]
        
        
        
        self.AvgTran=[ Normal(MakeName('AvgTran',self.Sites[i].SiteNumber),l[i],tau[i]  )   for i in range(nsite)]
        self.sArea=[Lambda(MakeName('sArea',self.SiteNumber[i]),lambda \
            a=self.AvgTran[i],w=w[i]: a*w) for i in range(nsite)]
        
        
 
      
      
    def SiteAreasFromData(self): 
        #Site areas
        nsite=len(self.SiteNumber)
        w=[ self.CoastLength[s]  for s in self.SiteNumber]
        l=[ t.AvgTranLen  for t in self.VirginSites]
        se=[ t.sterrTranLen  for t in self.VirginSites]
        tau=[1/s/s  for s in se]
        
        self.AvgTranLenFromData=self.AvgTranLen
        self.sAreaFromData=self.Area
        
        

if __name__ == "__main__":
  import pickle
  pfile="../MCMC/Jervis/Jervis.pickle"
  SiteNumber,Sites,VirginSites,CoastLength=pickle.load(open(pfile,"rb"))  
  
  ConstantValues.MinYear=1999  #Earliest year to be considered in calculations
  
  test=EFAProdModel(SiteNumber,Sites,VirginSites,CoastLength,oldB=1)
  #print(test.GetMeanWeight(16,[2005,2,14]))  

