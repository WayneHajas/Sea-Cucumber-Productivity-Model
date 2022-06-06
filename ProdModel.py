#run D:/Analyses/CukeEFAProd/pyfunctions/ProdModel
import pdb
import  numpy as np
from numpy import ndarray
from pymc import  *
from pylab import *
from tables import *
import tables
import warnings
warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)


sys.path.append('D:/Analyses/CukeEFAProd/pyfunctions')
from SiteNode import SiteNode
from TransectNode import TransectNode
from MonthTranNode import MonthTranNode
from location import location
from EFAsite import EFAsite
from transect import transect,MonthTran


class ProdModel():
  def __init__(self,mdbfile,Survey,SiteNum,CoastLength,Constants):
     '''ProdModel(mdbfile,Survey,SiteNum,CoastLength)
     Survey is the name of the survey
     SiteNum is a list of site-numbers for each site in the survey
     CoastLength is a list of lengths for each site in the survey'''
     
     self.mdbfile=mdbfile
     self.Survey=Survey
     self.CoastLength=CoastLength
     self.Constants=Constants
     
     self.L=location(mdbfile,Survey,SiteNum,CoastLength,self.Constants)
     self.DefineLocNodes()
     
     self.DefineSites(SiteNum)
     self.sites=[EFAsite(mdbfile,Survey,sn,self.CoastLength[i],self.Constants) for i ,sn in enumerate(SiteNum) ]
     self.SN=[SiteNode(self.LN,s)    for s in self.sites[:-1] ]
     LastSite=SiteNode(self.LN,self.sites[-1],OtherSiteNode=self.SN)
     self.SN=append(self.SN,LastSite)
     
     self.transect=[self.ReadTransects(s) for s in self.sites]
     self.TranNodes=self.GenTranNodes(self.SN,self.transect)
     
     self.MonthTranNode=self.GenMonthTranNode(self.LN,self.SN,self.TranNodes,self.transect)
     
     self.insertLnodes()
     del self.LN
     self.insertSnodes()
     del self.SN
     self.TranEff=  self.insertTnodes(self.TranNodes)
     del self.TranNodes
     
     self.PredCuke=self.insertPredCuke(self.MonthTranNode)
     self.ObsVal=self.insertObsVal(self.MonthTranNode)
     del self.MonthTranNode
     
  def DefineLocNodes(self):
    '''DefineLocNodes()'''
        
    self.xmax=Beta('xmax',1.1,1.1)
    self.binit=Beta('binit',1.1,1.1)
    self.b=Lambda('b',lambda bb=self.binit,xx=self.xmax:AdjB(bb,xx))
    self.a=Lambda('a',lambda bb=self.b,xx=self.xmax:bb*xx/(1-xx))
    self.fmax=Lognormal('fmax',-2.3,10)
    
    #Mean biomass in lg/metre-squared
    self.GrandMean=Lognormal('GrandMean',-1.,10.)
    
    	
    #Standard deviations of effects
    self.tauYear	=Lognormal('tauYear',		   3.,1.)
    self.tauSite	=Lognormal('tauSite',		   3.,1.)
    self.tauTransect	=Lognormal('tauTransect', 3.,1.)
    
    #Year Effects
    self.YearEff= [Normal('YearEff'+'_%i' %sy,  0.0, self.tauYear) for sy in self.L.SurveyYears[:-1]]
    LastYear=Lambda('YearEff'+'_%i' %self.L.SurveyYears[-1],lambda x=self.YearEff:-sum(x))
    self.YearEff=append(self.YearEff,LastYear)     


 	  

  def GetYearEff(self,Year):	
      result=None
      for i in range (self.L.nYear):
      	if Year==self.L.SurveyYears[i]:return(self.YearEff[i])
      return(None)
  def DefineSites(self,SiteNum):
     self.sites=[EFAsite(mdbfile,Survey,SiteNum[i],cl,self.Constants) for  i,cl in enumerate(self.CoastLength)     ]      
     self.SN=[ SiteNode(self,s)  for s in self.sites[:-1]]
     LastSite=SiteNode(self.LN,self.sites[-1],OtherSiteNode=self.SN)
     self.SN=append(self.SN,LastSite)     

  def insertSnodes(self):
     self.SiteEff=map(lambda s:s.SiteEff,self.SN)
     self.VBmass=map(lambda s:s.VBmass,self.SN)
     self.sArea=map(lambda s:s.sArea,self.SN)
     self.beta_area=map(lambda s:s.beta_area,self.SN)
     
     self.RelBmass=map(lambda s:s.RelBmass,self.SN)
     self.MeanWeight=map(lambda s:s.MeanWeight,self.SN)

  def insertTnodes(self,Tran):
    if isinstance(Tran,(list,ndarray)):return(map(lambda t: self.insertTnodes(t),Tran))
    return(Tran.TranEff)
  

  def GetYearEff(self,Year):	
          result=None
          for i in range (self.L.nYear):
          	if Year==self.L.SurveyYears[i]:return(self.YearEff[i])
          return(None)
     
  def ReadTransects(self,site):
    result=map(lambda tn:transect(self.mdbfile,self.Survey,site.Site,tn,self.Constants,WithTrim=True)     ,site.TranNum)
    return(result)
     
  def GenTranNodes(self,SN,tran):
    if isinstance(SN,(list,ndarray)):    return(map(lambda sn,t:self.GenTranNodes(sn,t)  ,SN,tran))
    if isinstance(tran,(list,ndarray)):  
      result=map(lambda    t:self.GenTranNodes(SN,t)  ,tran[:-1])
      last=TransectNode(self.LN,SN,tran[-1],OtherTran=result)
      result=append(result,last)
      return(result)
   
    return (TransectNode(self.LN,SN,tran))
    
  def GenMonthTranNode(self,LN,SN,TN,TRAN):
    if isinstance(SN,(list,ndarray)):    return(map(lambda sn,tn,t:self.GenMonthTranNode(LN,sn,tn,t)  ,SN,TN,TRAN))
    if isinstance(TN,(list,ndarray)):    return(map(lambda tn,t:   self.GenMonthTranNode(LN,SN,tn,t)     ,TN,TRAN))
    if TRAN.AllQuad!=[]:  
      return(map( lambda ym:MonthTranNode(LN,SN,TN,ym[1],self.Constants)   , TRAN.AllQuad) )
  
  def insertPredCuke(self,MTN):
    if isinstance(MTN,(list,ndarray)):    
      return(map(lambda mtn:self.insertPredCuke(mtn)  ,MTN))
    return (MTN.PredCuke)
 
  def insertObsVal(self,MTN):
    if isinstance(MTN,(list,ndarray)):    
      return(map(lambda mtn:self.insertObsVal(mtn)  ,MTN))
    return (MTN.ObsVal)
 
def AntiLogit(x):return(1./(1+exp(-x)))

def AdjB(initB,xmax):
  if xmax<.5:return(initB)
  return(initB*xmax)
     
if __name__ == "__main__":
  mdbfile='D:\Analyses\CukeEFAProd\data\SeaCuke_Bio_97.mdb'
  Survey='Tolmie Channel'
  SiteNum=[0,2,4,8,16]
  CoastLength=map(lambda x:1000,SiteNum)
  from MakeConstantsDict import MakeConstantsDict
  Constants=MakeConstantsDict()
  
  
  test=ProdModel(mdbfile,Survey,SiteNum,CoastLength,Constants)

