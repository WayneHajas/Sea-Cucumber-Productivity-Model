'''
20190405
    Special version where site, year and transect effects are not required to sum-to-zero

'''

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
    

        
    def SiteEffects(self):
        #Generate Site Effects
        self.SiteNode=[Normal(MakeName('site',s),0,self.tauSite) for s in self.SiteNumber    ]
        
       
    def YearEffects(self):
        #Generate Year Effects
        self.Year=self.GetYearOfSurvey()
        self.YearNode=[Normal(MakeName('year',s),0,self.tauYear) for s in self.Year    ]    



    def TransectEffects(self):
        self.TranEffect=[] #Transect Effects   
        self.DummyTran=[]  #Used to impose prior distributions on last transect of site 
        for s in self.SiteNumber:
            SE=self.GetSiteNode(s) #Site Effect
            curSite=self.GetSite(s)
            
            #Generate transect-effects for the site 
            CurTran=[ Normal(MakeName('TranEffect',[s,t]),0,self.tauTransect)  for t in curSite.TranNum] #All but the last transect in the site
            self.TranEffect+=  CurTran       


            

if __name__ == "__main__":
  import pickle
  pfile="../MCMC/Jervis/Jervis.pickle"
  SiteNumber,Sites,VirginSites,CoastLength=pickle.load(open(pfile,"rb"))  
  
  ConstantValues.MinYear=1999  #Earliest year to be considered in calculations
  
  test=EFAProdModel(SiteNumber,Sites,VirginSites,CoastLength,oldB=1)
  #print(test.GetMeanWeight(16,[2005,2,14]))  

