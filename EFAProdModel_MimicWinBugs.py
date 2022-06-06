'''2016-01-06
       Force 
         *Site Effects to sum to zero
         *Year Effects to sum to zero
         *Transect Effects to sum to zero for each site
2016-01-12
    Use pymc.Deterministic to explicitly incorporate the relative-biomasses as nodes.
2016-01-13
    Get all the nodes defined in __init__    
    I would like to break up __init__ but some nodes do not seem to get sampled properly.
2016-02-15
   Make the prior distributions more like the old WinBugs version
   F:\Archive\s-drive\analyses\CukeExpHarvRecov\Zeballos.beta\ZeballosSpatial.odc for example
   Corrected serious error in calculation of virgin biomass.  Now use the zeroth year-effect
2016-02-15
    Start from the new implementatino of the WCh model and create teh Gram-schaeffer model
2016-02-29
	Apply a prior specifying that virgin grand-mean must be approximate to the nonvirgin grand-mean.
	Apply a prior to sdYear to try and make it small
 2016-03-04
     Implemented virgin transect-counts as random variables
2016-03-29
    VPredVal:  make conversion between population and biomass 
2017-01-17
    Change the name of the model to EFAProdModel and start a new implementation
    
2017-02-06
    Change constructor-variables so be SiteNumber,Sites,VirginSites instead of mdbfile,Survey,CoastLength
    
    Previous versions read from a .mdb file - and that required a 32-bit version of python.
    This version assumes recieves data in the form of instances of classes.  No need to read .mdb files.  Can use the 64-bit version of python.
2017-03-13
    Do not create year-effects for years where there was harvest but no survey    
2017-03-24
    Apply prior-distributions to relative-biomass values.  
    dummy~normal(log(RBM),0.754)   
    
2017-04-07
    Replaced in error in cacluating b from binit and xmax
2017-04-11
    Error Correction.  Relative biomass included in calculation of predicted number of animals on transect

2017-06-30
    New Version: EFAProdModel_MimicWinBugs
    
    The intent here is to be able to use the same data as went into the WinBugs model.
    Data is taken from a library
    Code mimics the WinBugs model as closely as possible
2017-07-07
    Try changing the prior on xmax to encourage small values
    
2017-07-13
    Trace average weights
2017-09-14
    Refactored the code

'''

from pymc import  *
from pylab import *
from tables import *
from numpy import log,average,exp
import tables
import warnings
from datetime import datetime as dt
from ADO import adoBaseClass as daoBaseClass
from ToDecimalYear import ToDecimalYear,FromDecimalYear
from constants import ConstantValues
from EFAsite import  Site,VirginSite
from  Harvest import  Harvest
from MeanWeight import MeanWeight
from NodeNameToDecimalYear import MatchNodeRefTime
from NextRelBmass import NextRelBmass
from MakeName import MakeName
from EFAProdModel_fixSD import EFAProdModel as oldEFAProdModel
from EFAProdModel_fixSD import NormalizeEffect



class EFAProdModel(oldEFAProdModel):
    def __init__(self,datalib,oldB=1):
        '''EFAProdModel(datalib)
        * datlib is a library with data corresponding to old WinBUGS analysis
        
         '''
        self.datalib=datalib
        self.nsite=self.datalib.nsite
        self.MaxnSurvey=self.datalib.MaxnSurvey
        self.Nmonth=self.datalib.Nmonth
        self.CoastLength=self.datalib.CoastLength
        self.oldB=oldB
        
        self.ProductivityParameters()  #Generate nodes for productivity function       
        self.GrandMeans()              #Generate nodes for Grand Mean and virgin Grand Mean of biomass density
        self.GrandMean=Lambda('GrandMean',lambda t=self.lnGrandMean:exp(t))
        self.EffectStandardDeviations()#Generate nodes for standard deviations of effects
        self.SiteEffects()             #                   site effects
        self.TransectEffects()         #                   transect effects
        self.YearEffects()             #                   Year Effects        
        self.SiteAreas()               #                   Site Areas
        self.VirginBiomass()           #                   Virgin biomass
        self.RelativeHarvest()         #                   Harvest as a fraction of virgin biomass
        self.CalcRelAbund()            #                   Relative Abundance
        
        self.MeanWeight()
        self.VTransects()              #Virgin Transects
        self.STransects()              #Survey Transects
        
          
  

      
 
    def SiteEffects(self):
        #Generate Site Effects        
        self.oriSite=[Normal(MakeName('orisite',s),0,1,trace=False) for s in range(self.nsite)    ]        
        self.OriSiteMean=Lambda('OrisiteMean', lambda x=self.oriSite:average(x),trace=False)
        self.OriSiteStdv=Lambda('OrisiteStdv', lambda x=self.oriSite:std(x),trace=False)
        self.SiteNode=[Lambda(OSE.__name__[3:], lambda t=OSE,mu=self.OriSiteMean,sigma1=self.OriSiteStdv,sigma2=self.sdSite:NormalizeEffect(t,mu,sigma1,sigma2),trace=True)  for OSE in self.oriSite]

                          
    def YearEffects(self):
        #Year Effects have mean of zero and force the estimated standard deviation to be sdYear
        self.oriYear=[Normal(MakeName('oriYear',i),0,1,trace=False)   for i in range(self.datalib.MaxnSurvey)]           
        self.OriYearMean=Lambda('OriYearMean', lambda x=self.oriYear:average(x),trace=False)
        self.OriYearStdv=Lambda('OriYearStdv', lambda x=self.oriYear:std(x),trace=False)
        self.Year=[Lambda(MakeName('Year',i), lambda t=self.oriYear[i],mu=self.OriYearMean,sigma1=self.OriYearStdv,sigma2=self.sdYear:(t-mu)/sigma1*sigma2,trace=True)  for i in range(self.datalib.MaxnSurvey)]


           

    def SiteAreas(self): 
        self.Tlbeta=[]      
        self.AvgNumQuad=[]
        self.SiteArea=[]
        self.TLbeta=[]
        for i in range(self.datalib.nsite):
            self.TLbeta+=[Beta(MakeName('TLbeta',i),1.5,1.55)]
            self.AvgNumQuad+=[Lambda(MakeName('AvgNumQuad',i),lambda TLbeta=self.TLbeta[i],PRTranTau=self.datalib.PRTranTau[i],PRTranAvg=self.datalib.PRTranAvg[i]: (TLbeta-.5) *4/sqrt( PRTranTau)+PRTranAvg)]
            self.SiteArea+=[Lambda(MakeName('SiteArea',i),lambda CoastLength=self.datalib.CoastLength[i],AvgNumQuad=self.AvgNumQuad[i]: CoastLength*5*AvgNumQuad )]
            
 
    def VirginBiomass(self):
        #Virgin Biomass in kilograms
        self.VBmass=[Lambda(MakeName('VBmass',s),\
            lambda VGrandMean=self.VGrandMean,SiteEff=self.SiteNode[s],YearEff=self.Year[0],SiteArea=self.SiteArea[s]:VGrandMean*exp(SiteEff+YearEff)*SiteArea)\
            for s in range(self.datalib.nsite)]
                                    
 

    def TransectEffects(self):
        #For individual sites, force the mean of transect effects to be zero and the estimated standard deviation to be sdTran
        self.oriTran00=[Normal(MakeName('oriTran00',i),0,1,trace=False)   for i in range(self.datalib.TranPerSite[0])]           
        self.oriTran02=[Normal(MakeName('oriTran02',i),0,1,trace=False)   for i in range(self.datalib.TranPerSite[1])]           
        self.oriTran04=[Normal(MakeName('oriTran04',i),0,1,trace=False)   for i in range(self.datalib.TranPerSite[2])]           
        self.oriTran08=[Normal(MakeName('oriTran08',i),0,1,trace=False)   for i in range(self.datalib.TranPerSite[3])]           
        self.oriTran16=[Normal(MakeName('oriTran16',i),0,1,trace=False)   for i in range(self.datalib.TranPerSite[4])]       
        
        self.OriTran00Mean=Lambda('OriTran00Mean', lambda x=self.oriTran00:average(x),trace=False)
        self.OriTran02Mean=Lambda('OriTran02Mean', lambda x=self.oriTran02:average(x),trace=False)
        self.OriTran04Mean=Lambda('OriTran04Mean', lambda x=self.oriTran04:average(x),trace=False)
        self.OriTran08Mean=Lambda('OriTran08Mean', lambda x=self.oriTran08:average(x),trace=False)
        self.OriTran16Mean=Lambda('OriTran16Mean', lambda x=self.oriTran16:average(x),trace=False)      
        
        self.OriTran00Stdv=Lambda('OriTran00Stdv', lambda x=self.oriTran00:std(x),trace=False)
        self.OriTran02Stdv=Lambda('OriTran02Stdv', lambda x=self.oriTran02:std(x),trace=False)
        self.OriTran04Stdv=Lambda('OriTran04Stdv', lambda x=self.oriTran04:std(x),trace=False)
        self.OriTran08Stdv=Lambda('OriTran08Stdv', lambda x=self.oriTran08:std(x),trace=False)
        self.OriTran16Stdv=Lambda('OriTran16Stdv', lambda x=self.oriTran16:std(x),trace=False)

        self.Tran00=[ Lambda(MakeName('Tran00',i), lambda t=self.oriTran00[i],mu=self.OriTran00Mean,sigma1=self.OriTran00Stdv,sigma2=self.sdTransect:(t-mu)/sigma1*sigma2)  for i in range(self.datalib.TranPerSite[0])]
        self.Tran02=[ Lambda(MakeName('Tran02',i), lambda t=self.oriTran02[i],mu=self.OriTran02Mean,sigma1=self.OriTran02Stdv,sigma2=self.sdTransect:(t-mu)/sigma1*sigma2)  for i in range(self.datalib.TranPerSite[1])]
        self.Tran04=[ Lambda(MakeName('Tran04',i), lambda t=self.oriTran04[i],mu=self.OriTran04Mean,sigma1=self.OriTran04Stdv,sigma2=self.sdTransect:(t-mu)/sigma1*sigma2)  for i in range(self.datalib.TranPerSite[2])]
        self.Tran08=[ Lambda(MakeName('Tran08',i), lambda t=self.oriTran08[i],mu=self.OriTran08Mean,sigma1=self.OriTran08Stdv,sigma2=self.sdTransect:(t-mu)/sigma1*sigma2)  for i in range(self.datalib.TranPerSite[3])]
        self.Tran16=[ Lambda(MakeName('Tran16',i), lambda t=self.oriTran16[i],mu=self.OriTran16Mean,sigma1=self.OriTran16Stdv,sigma2=self.sdTransect:(t-mu)/sigma1*sigma2)  for i in range(self.datalib.TranPerSite[4])]
                
    def CalcRelAbund(self):
        #Relative biomass.  Harvest is assumed to occur immediately after a survey
        self.RelBmass=[]
        for s in range(self.datalib.nsite):
            CurSite=[]
            m=0
            CurSite+=[Lambda(MakeName('RelBmass',[s,m]),lambda dummy=self.VBmass[0]:0.999)]#relative biomass starts out at near-one
            for m in range(1,self.datalib.Nmonth):
                CurSite+=[Lambda(MakeName('RelBmass',[s,m]),\
                lambda OldBiomass=CurSite[-1],OldTime=0,NewTime=1./12, a=self.a,b=self.b,xmax=self.xmax, fmax=self.fmax,VBiomass=1,CurHarvest=self.RelHarv[s][m-1]:\
                    NextRelBmass(OldBiomass=OldBiomass,OldTime=OldTime,NewTime=NewTime,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBiomass,CurHarvest=CurHarvest),trace=m in self.datalib.SurveyIndex)]
            self.RelBmass+=[CurSite]
 
    def MeanWeight(self):
        #In grams        
        self.AvgWeight=[  \
                        [TruncatedNormal(MakeName('AvgWeight',[s,j]),self.datalib.MnWgt[s][j],self.datalib.TauWgt[s][j],a=1,b=1000,trace=True)   for j in range(self.datalib.MaxnSurvey)]\
                        for s in range(self.datalib.nsite) ]
        

    def RelativeHarvest(self):
        #Harvest as a fraction of virgin biomass
        self.RelHarv=[
            [Lambda(MakeName('RelHarv',[s,m]),lambda Harvest=self.datalib.Harvest[s][m], VBmass=self.VBmass[s]:Harvest/VBmass,trace=(self.datalib.Harvest[s][m]>0))\
            for m in range(self.datalib.Nmonth)] \
            for s in range(self.datalib.nsite)]


            
    def VTransects(self):            
        #Site 0 virgin(untrimmed) predicted and observed counts
        s=0
        self.PredVPop00=[
                        Lambda(MakeName('PredVPop00',[s,t]),\
                           lambda GrandMean=self.VGrandMean,nQuad=self.datalib.QuadNoTrim[s][t], \
                               SiteEff=self.SiteNode[s],YearEff=self.Year[0],TranEff=self.Tran00[t],AvgWeight=self.AvgWeight[s][0]:
                           GrandMean*nQuad*20*\
									exp(SiteEff+YearEff+TranEff)/
												AvgWeight*1000 ,trace=False)                     
                       for t in range(self.datalib.TranPerSite[s])]
        self.ObsVCuke00=[]
        for t in range(self.datalib.TranPerSite[s]):
            try:
              self.ObsVCuke00+=[Poisson(MakeName('ObsVCuke00',[s,t]),self.PredVPop00[t],\
                           observed=True,value=self.datalib.NoTrim[s][t]) ]
            except:
              self.ObsVCuke00+=[Lambda(MakeName('ObsVCuke00',[s,t]),lambda dummy=0:-1,trace=False)]                     
        
        
                      
        #Site 2 virgin(untrimmed) predicted and observed counts
        s=1
        self.PredVPop02=[
                        Lambda(MakeName('PredVPop02',[s,t]),\
                           lambda GrandMean=self.VGrandMean,nQuad=self.datalib.QuadNoTrim[s][t], \
                               SiteEff=self.SiteNode[s],YearEff=self.Year[0],TranEff=self.Tran02[t],AvgWeight=self.AvgWeight[s][0]:
                           GrandMean*nQuad*20*\
									exp(SiteEff+YearEff+TranEff)/
												AvgWeight*1000 ,trace=False )                     
                       for t in range(self.datalib.TranPerSite[s])]
        self.ObsVCuke02=[]
        for t in range(self.datalib.TranPerSite[s]):
            try:
              self.ObsVCuke02+=[Poisson(MakeName('ObsVCuke02',[s,t]),self.PredVPop02[t],\
                           observed=True,value=self.datalib.NoTrim[s][t]) ]
            except:
              self.ObsVCuke02+=[Lambda(MakeName('ObsVCuke02',[s,t]),lambda dummy=0:-1,trace=False)]    
                      
        #Site 4 virgin(untrimmed) predicted and observed counts
        s=2
        self.PredVPop04=[
                        Lambda(MakeName('PredVPop04',[s,t]),\
                           lambda GrandMean=self.VGrandMean,nQuad=self.datalib.QuadNoTrim[s][t], \
                               SiteEff=self.SiteNode[s],YearEff=self.Year[0],TranEff=self.Tran04[t],AvgWeight=self.AvgWeight[s][0]:
                           GrandMean*nQuad*20*\
									exp(SiteEff+YearEff+TranEff)/
												AvgWeight*1000 ,trace=False )                     
                       for t in range(self.datalib.TranPerSite[s])]
        self.ObsVCuke04=[]
        for t in range(self.datalib.TranPerSite[s]):
            try:
              self.ObsVCuke04+=[Poisson(MakeName('ObsVCuke04',[s,t]),self.PredVPop04[t],\
                           observed=True,value=self.datalib.NoTrim[s][t]) ]
            except:
              self.ObsVCuke04+=[Lambda(MakeName('ObsVCuke04',[s,t]),lambda dummy=0:-1,trace=False)]    
                      
        #Site 8 virgin(untrimmed) predicted and observed counts
        s=3
        self.PredVPop08=[
                        Lambda(MakeName('PredVPop08',[s,t]),\
                           lambda GrandMean=self.VGrandMean,nQuad=self.datalib.QuadNoTrim[s][t], \
                               SiteEff=self.SiteNode[s],YearEff=self.Year[0],TranEff=self.Tran08[t],AvgWeight=self.AvgWeight[s][0]:
                           GrandMean*nQuad*20*\
									exp(SiteEff+YearEff+TranEff)/
												AvgWeight*1000  ,trace=False)                     
                       for t in range(self.datalib.TranPerSite[s])]
        self.ObsVCuke08=[]
        for t in range(self.datalib.TranPerSite[s]):
            try:
              self.ObsVCuke08+=[Poisson(MakeName('ObsVCuke08',[s,t]),self.PredVPop08[t],\
                           observed=True,value=self.datalib.NoTrim[s][t]) ]
            except:
              self.ObsVCuke08+=[Lambda(MakeName('ObsVCuke08',[s,t]),lambda dummy=0:-1,trace=False)]    
        
        
                      
        #Site 16 virgin(untrimmed) predicted and observed counts
        s=4
        self.PredVPop16=[
                        Lambda(MakeName('PredVPop16',[s,t]),\
                           lambda GrandMean=self.VGrandMean,nQuad=self.datalib.QuadNoTrim[s][t], \
                               SiteEff=self.SiteNode[s],YearEff=self.Year[0],TranEff=self.Tran16[t],AvgWeight=self.AvgWeight[s][0]:
                           GrandMean*nQuad*20*\
									exp(SiteEff+YearEff+TranEff)/
												AvgWeight*1000 ,trace=False )                     
                       for t in range(self.datalib.TranPerSite[s])]
        self.ObsVCuke16=[]
        for t in range(self.datalib.TranPerSite[s]):
            try:
              self.ObsVCuke16+=[Poisson(MakeName('ObsVCuke16',[s,t]),self.PredVPop16[t],\
                           observed=True,value=self.datalib.NoTrim[s][t]) ]
            except:
              self.ObsVCuke16+=[Lambda(MakeName('ObsVCuke16',[s,t]),lambda dummy=0:-1,trace=False)]    
                      
             
        
    def STransects(self):
        #Predicted and observed number of cukes:  Site 0
        s=0
        self.PredPop00=[]
        self.ObsCuke00=[]
        for t in range(self.datalib.TranPerSite[s]):
            CurPred=[]
            CurObs=[]
            for y in range(self.datalib.MaxnSurvey):
                GrandMean=self.GrandMean
                nQuad=self.datalib.Quad00[t]
                RelBmass=self.RelBmass[s][self.datalib.SurveyIndex[y]]
                SiteEff=self.SiteNode[s]
                YearEff=self.Year[y]
                TranEff=self.Tran00[t]
                AvgWeight=self.AvgWeight[s][y]
                CurPred+=[Lambda(MakeName('PredPop00',[t,self.datalib.SurveyIndex[y]]),\
                           lambda GrandMean=GrandMean,nQuad=nQuad, RelBmass=RelBmass,\
                               SiteEff=SiteEff,YearEff=YearEff,TranEff=TranEff,AvgWeight=AvgWeight:
                                   GrandMean*nQuad*20*RelBmass*exp(SiteEff+YearEff+TranEff)/AvgWeight*1000 ,trace=False )  ]    
                try:            
                   CurObs+=[Poisson(MakeName('ObsCuke00',[t,self.datalib.SurveyIndex[y]]),CurPred[-1],\
                           observed=True,value=self.datalib.Cuke00[t][y])  ]  
                except:
                   CurObs+=[Lambda(MakeName('PredPop00',[t,self.datalib.SurveyIndex[y]]),lambda dummy=0:-1,trace=False)  ]
            self.PredPop00+=[CurPred]
            self.ObsCuke00+=[CurObs]
            
        #Predicted and observed number of cukes:  Site 2
        s=1
        self.PredPop02=[]
        self.ObsCuke02=[]
        for t in range(self.datalib.TranPerSite[s]):
            CurPred=[]
            CurObs=[]
            for y in range(self.datalib.MaxnSurvey):
                GrandMean=self.GrandMean
                nQuad=self.datalib.Quad02[t]
                RelBmass=self.RelBmass[s][self.datalib.SurveyIndex[y]]
                SiteEff=self.SiteNode[s]
                YearEff=self.Year[y]
                TranEff=self.Tran02[t]
                AvgWeight=self.AvgWeight[s][y]
                CurPred+=[Lambda(MakeName('PredPop02',[t,self.datalib.SurveyIndex[y]]),\
                           lambda GrandMean=GrandMean,nQuad=nQuad, RelBmass=RelBmass,\
                               SiteEff=SiteEff,YearEff=YearEff,TranEff=TranEff,AvgWeight=AvgWeight:
                                   GrandMean*nQuad*20*RelBmass*exp(SiteEff+YearEff+TranEff)/AvgWeight*1000 ,trace=False )  ]    
                try:            
                   CurObs+=[Poisson(MakeName('ObsCuke02',[t,self.datalib.SurveyIndex[y]]),CurPred[-1],\
                           observed=True,value=self.datalib.Cuke02[t][y])  ]  
                except:
                   CurObs+=[Lambda(MakeName('PredPop02',[t,self.datalib.SurveyIndex[y]]),lambda dummy=0:-1,trace=False)  ]
            self.PredPop02+=[CurPred]
            self.ObsCuke02+=[CurObs]
            
        #Predicted and observed number of cukes:  Site 4
        s=2
        self.PredPop04=[]
        self.ObsCuke04=[]
        for t in range(self.datalib.TranPerSite[s]):
            CurPred=[]
            CurObs=[]
            for y in range(self.datalib.MaxnSurvey):
                GrandMean=self.GrandMean
                nQuad=self.datalib.Quad04[t]
                RelBmass=self.RelBmass[s][self.datalib.SurveyIndex[y]]
                SiteEff=self.SiteNode[s]
                YearEff=self.Year[y]
                TranEff=self.Tran04[t]
                AvgWeight=self.AvgWeight[s][y]
                CurPred+=[Lambda(MakeName('PredPop04',[t,self.datalib.SurveyIndex[y]]),\
                           lambda GrandMean=GrandMean,nQuad=nQuad, RelBmass=RelBmass,\
                               SiteEff=SiteEff,YearEff=YearEff,TranEff=TranEff,AvgWeight=AvgWeight:
                                   GrandMean*nQuad*20*RelBmass*exp(SiteEff+YearEff+TranEff)/AvgWeight*1000,trace=False  )  ]    
                try:            
                   CurObs+=[Poisson(MakeName('ObsCuke04',[t,self.datalib.SurveyIndex[y]]),CurPred[-1],\
                           observed=True,value=self.datalib.Cuke04[t][y])  ]  
                except:
                   CurObs+=[Lambda(MakeName('ObsCuke04',[t,self.datalib.SurveyIndex[y]]), lambda dummy=0:-1,trace=False)  ]  
            self.PredPop04+=[CurPred]
            self.ObsCuke04+=[CurObs]
            
        #Predicted and observed number of cukes:  Site 8
        s=3
        self.PredPop08=[]
        self.ObsCuke08=[]
        for t in range(self.datalib.TranPerSite[s]):
            CurPred=[]
            CurObs=[]
            for y in range(self.datalib.MaxnSurvey):
                GrandMean=self.GrandMean
                nQuad=self.datalib.Quad08[t]
                RelBmass=self.RelBmass[s][self.datalib.SurveyIndex[y]]
                SiteEff=self.SiteNode[s]
                YearEff=self.Year[y]
                TranEff=self.Tran08[t]
                AvgWeight=self.AvgWeight[s][y]
                CurPred+=[Lambda(MakeName('PredPop08',[t,self.datalib.SurveyIndex[y]]),\
                           lambda GrandMean=GrandMean,nQuad=nQuad, RelBmass=RelBmass,\
                               SiteEff=SiteEff,YearEff=YearEff,TranEff=TranEff,AvgWeight=AvgWeight:
                                   GrandMean*nQuad*20*RelBmass*exp(SiteEff+YearEff+TranEff)/AvgWeight*1000 ,trace=False )  ]       
                try:            
                   CurObs+=[Poisson(MakeName('ObsCuke08',[t,self.datalib.SurveyIndex[y]]),CurPred[-1],\
                           observed=True,value=self.datalib.Cuke08[t][y])  ]  
                except:
                   CurObs+=[Lambda(MakeName('PredPop08',[t,self.datalib.SurveyIndex[y]]),lambda dummy=0:-1,trace=False)]
            
        #Predicted and observed number of cukes:  Site 16
        s=4
        self.PredPop16=[]
        self.ObsCuke16=[]
        for t in range(self.datalib.TranPerSite[s]):
            CurPred=[]
            CurObs=[]
            for y in range(self.datalib.MaxnSurvey):
                GrandMean=self.GrandMean
                nQuad=self.datalib.Quad16[t]
                RelBmass=self.RelBmass[s][self.datalib.SurveyIndex[y]]
                SiteEff=self.SiteNode[s]
                YearEff=self.Year[y]
                TranEff=self.Tran16[t]
                AvgWeight=self.AvgWeight[s][y]
                CurPred+=[Lambda(MakeName('PredPop16',[t,self.datalib.SurveyIndex[y]]),\
                           lambda GrandMean=GrandMean,nQuad=nQuad, RelBmass=RelBmass,\
                               SiteEff=SiteEff,YearEff=YearEff,TranEff=TranEff,AvgWeight=AvgWeight:
                                   GrandMean*nQuad*20*RelBmass*exp(SiteEff+YearEff+TranEff)/AvgWeight*1000,trace=False  )  ]    
                try:            
                   CurObs+=[Poisson(MakeName('ObsCuke16',[t,self.datalib.SurveyIndex[y]]),CurPred[-1],\
                           observed=True,value=self.datalib.Cuke16[t][y])  ]  
                except:
                   CurObs+=[Lambda(MakeName('PredPop16',[t,self.datalib.SurveyIndex[y]]),lambda dummy=0:-1,trace=False)]
            self.PredPop16+=[CurPred]
            self.ObsCuke16+=[CurObs]
                
                             
           
            
if __name__ == "__main__":
  
  
  import sys
  sys.path.append('D:\\Analyses\\CukeEFAProd\\data\\')
  import  JervisWinBugsData as datalib

  test=EFAProdModel(datalib)
  #print(test.GetMeanWeight(16,[2005,2,14]))  

