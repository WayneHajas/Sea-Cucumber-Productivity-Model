'''
20200113
    Build more routines so that minimum linear and spatial, biomass and population, densities can be taken
    I am just using dummy values for coastlength
'''

from numpy import exp

from GetNames import *
from GetParamStats import GetParamValues2 as GetParamValues
from NextRelBmass import NextRelBmass,ApplyBounds
from GetCoastLength import GetCoastLength
from ToDecimalYear import ToDecimalYear



def GetRel_Abund_Names(mdbfile,root='RelBmass'):
    pname=GetNames(mdbfile)
    pname=[t for t in pname if t[:len(root)]==root]
    return(pname)
    
def GetAllRelAbund(mdbfile,burn=0,nthin=1000,root='RelBmass'):
    
    #Parameter names associated with relative biomass    
    rname=GetRel_Abund_Names(mdbfile,root=root)
    result=GetParamValues(mdbfile, rname,burn=burn,nthin=nthin)
    return(result)    

def GetMinRelAbund(hdf5file,burn=0,nthin=1000,minRelBmass=1.e-3,maxRelBmass=1-1e-3,maxTime=9999,root='RelBmass'):
    RelAbund=GetAllRelAbund(hdf5file,burn=burn,nthin=nthin,root=root)
    result= [min([t[i]  for t in RelAbund ])   for i,dummy in enumerate(RelAbund[0])]
    return(result)    
    
def GetRelAbund_WithName(hdf5file,burn=0,nthin=1000,root='RelBmass'):
    Rel_Abund_Names=GetRel_Abund_Names(hdf5file,root=root)
    RelAbund=GetAllRelAbund(hdf5file,burn=burn,nthin=nthin,root=root)
    nname=len(Rel_Abund_Names)
    compo=[]
    AllSite=[]
    for i,curname in enumerate(Rel_Abund_Names):
        CurResult={}
        site,year,month,day=[int(t) for t in curname.split('_')[-4:]]
        AllSite+=[site]
        DecimalYear=ToDecimalYear([year,month,day])
        CurResult['values']=RelAbund[i]
        CurResult['DecimalYear']=DecimalYear
        CurResult['year']=year
        CurResult['month']=month
        CurResult['day']=day
        CurResult['site']=site
        compo+=[CurResult]
    uniSite=list(set(AllSite))
    result={}
    for s in uniSite:
        result[s]=[t for t in compo if t['site']==s   ]
    return(result)
        
def GetVirginBiomass_WithName(hdf5file,burn=0,nthin=1000,root='VBmass'):
      pname=GetRel_Abund_Names(hdf5file,root=root)
      site=[int(curname.split('_')[-1]) for curname in pname]
      result={}
      for s in site:
          result[s]=GetParamValues(hdf5file, root+'_'+str(s),burn=burn,nthin=nthin)
      return(result)
      
def IndexToMeanWeight(RelBM,MW):
    '''
    Generate indeces to connect releative biomass values to corresponding mean-weight values
    Both RelBM and MW are results of GetRelAbund_WithName'''
    site=list(RelBM.keys())
    result={}
    for s in site:
        curIndex=[]
        for t in RelBM[s]:
            ssq=[ (t['DecimalYear'] - r['DecimalYear'])**2   for r in MW[s] ]
            curIndex+=[ssq.index(min(ssq))]
        result[s]=curIndex
    return(result)
   
def LinearBiomassDensity(hdf5file,survey,burn=0,nthin=1000):
    RelBM=GetRelAbund_WithName(hdf5file,burn=burn,nthin=nthin,root='Rel_Abund')
    MW=GetRelAbund_WithName(hdf5file,burn=burn,nthin=nthin,root='MeanWeight')
    Rel_MW_index=IndexToMeanWeight(RelBM,MW)
    VBiomass=GetVirginBiomass_WithName(hdf5file,burn=burn,nthin=nthin,root='VBmass')
    SiteEffect=GetVirginBiomass_WithName(hdf5file,burn=burn,nthin=nthin,root='site')
    YearEffect=GetVirginBiomass_WithName(hdf5file,burn=burn,nthin=nthin,root='year')
    sigmaTransect=GetParamValues(hdf5file,'sdTransect',burn=burn,nthin=nthin)
    
    CoastLength=GetCoastLength(survey,mdbfile='S:\\analyses\\EFA_Productivity_Model.20180929\\data\\CoastLength.mdb')
    
    n=len(VBiomass[0])
    site=list(RelBM.keys())
    result={}
    for s in site:
        result[s]=[]
        VB=VBiomass[s]
        SE=SiteEffect[s]
        curCL=CoastLength[s]
        for i,r in enumerate(RelBM[s]):
            curMW=Rel_MW_index[s][i]
            try:
                curYE=YearEffect[r['year']]
                CurResult={}
                CurResult['DecimalYear']=r['DecimalYear']
                CurResult['year']=r['year']
                CurResult['month']=r['month']
                CurResult['day']=r['day']
                CurResult['site']=r['site']
                
                CurResult['values']=[r['values'][j]*VB[j]/curCL*exp(SE[j]+curYE[j]+(sigmaTransect[j]**2)/2)   for j in range(n)]
                result[s]+=[CurResult]
            except:
                #No survey-results(year effect) for significant event
                dummy=True
    return(result)
   
def LinearPopDensity(hdf5file,survey,burn=0,nthin=1000):
    RelBM=GetRelAbund_WithName(hdf5file,burn=burn,nthin=nthin,root='Rel_Abund')
    MW=GetRelAbund_WithName(hdf5file,burn=burn,nthin=nthin,root='MeanWeight')
    Rel_MW_index=IndexToMeanWeight(RelBM,MW)
    VBiomass=GetVirginBiomass_WithName(hdf5file,burn=burn,nthin=nthin,root='VBmass')
    SiteEffect=GetVirginBiomass_WithName(hdf5file,burn=burn,nthin=nthin,root='site')
    YearEffect=GetVirginBiomass_WithName(hdf5file,burn=burn,nthin=nthin,root='year')
    sigmaTransect=GetParamValues(hdf5file,'sdTransect',burn=burn,nthin=nthin)
    
    CoastLength=GetCoastLength(survey,mdbfile='S:\\analyses\\EFA_Productivity_Model.20180929\\data\\CoastLength.mdb')
    
    n=len(VBiomass[0])
    site=list(RelBM.keys())
    result={}
    for s in site:
        result[s]=[]
        VB=VBiomass[s]
        SE=SiteEffect[s]
        curCL=CoastLength[s]
        for i,r in enumerate(RelBM[s]):
            try:
                CurResult={}
                curYE=YearEffect[r['year']]
                CurResult['DecimalYear']=r['DecimalYear']
                CurResult['year']=r['year']
                CurResult['month']=r['month']
                CurResult['day']=r['day']
                CurResult['site']=r['site']
                
                CurWeight=MW[s] [Rel_MW_index[s] [i]]['values']
                
                CurResult['values']=[r['values'][j]*VB[j]/curCL/CurWeight[j]*exp(SE[j]+curYE[j]+(sigmaTransect[j]**2)/2)   for j in range(n)]
                result[s]+=[CurResult]
            except:
                #No survey-results(year effect) for significant event
                dummy=True
    return(result)
   
def SpatialBiomassDensity(hdf5file,survey,burn=0,nthin=1000):
    RelBM=GetRelAbund_WithName(hdf5file,burn=burn,nthin=nthin,root='Rel_Abund')
    SiteEffect=GetVirginBiomass_WithName(hdf5file,burn=burn,nthin=nthin,root='site')
    YearEffect=GetVirginBiomass_WithName(hdf5file,burn=burn,nthin=nthin,root='year')
    sigmaTransect=GetParamValues(hdf5file,'sdTransect',burn=burn,nthin=nthin)
    
    lnGrandMean=GetParamValues(hdf5file, 'lnGrandMean',burn=burn,nthin=nthin)
    GrandMean=[exp(t)  for t in lnGrandMean]
    
    n=len(GrandMean)
    site=list(RelBM.keys())
    result={}
    for s in site:
        result[s]=[]
        SE=SiteEffect[s]
        for i,r in enumerate(RelBM[s]):
            try:
                CurResult={}
                curYE=YearEffect[r['year']]
                CurResult['DecimalYear']=r['DecimalYear']
                CurResult['year']=r['year']
                CurResult['month']=r['month']
                CurResult['day']=r['day']
                CurResult['site']=r['site']
                
                CurResult['values']=[r['values'][j]*GrandMean[j]*exp(SE[j]+curYE[j]+(sigmaTransect[j]**2)/2)   for j in range(n)]
                result[s]+=[CurResult]
            except:
                #No survey-results(year effect) for significant event
                dummy=True
    return(result)
   
def SpatialPopDensity(hdf5file,survey,burn=0,nthin=1000):
    RelBM=GetRelAbund_WithName(hdf5file,burn=burn,nthin=nthin,root='Rel_Abund')
    MW=GetRelAbund_WithName(hdf5file,burn=burn,nthin=nthin,root='MeanWeight')
    Rel_MW_index=IndexToMeanWeight(RelBM,MW)
    SiteEffect=GetVirginBiomass_WithName(hdf5file,burn=burn,nthin=nthin,root='site')
    YearEffect=GetVirginBiomass_WithName(hdf5file,burn=burn,nthin=nthin,root='year')
    sigmaTransect=GetParamValues(hdf5file,'sdTransect',burn=burn,nthin=nthin)
    
    lnGrandMean=GetParamValues(hdf5file, 'lnGrandMean',burn=burn,nthin=nthin)
    GrandMean=[exp(t)  for t in lnGrandMean]
    
    n=len(GrandMean)
    site=list(RelBM.keys())
    result={}
    for s in site:
        result[s]=[]
        SE=SiteEffect[s]
        for i,r in enumerate(RelBM[s]):
            try:
                CurResult={}
                curYE=YearEffect[r['year']]
                CurResult['DecimalYear']=r['DecimalYear']
                CurResult['year']=r['year']
                CurResult['month']=r['month']
                CurResult['day']=r['day']
                CurResult['site']=r['site']
                                
                CurWeight=MW[s] [Rel_MW_index[s] [i]]['values']

                CurResult['values']=[r['values'][j]*GrandMean[j]/CurWeight[j]*exp(SE[j]+curYE[j]+(sigmaTransect[j]**2)/2)   for j in range(n)]
                result[s]+=[CurResult]
            except:
                #No survey-results(year effect) for significant event
                dummy=True
    return(result)

def FindMinima(abundance):
    n=len (abundance[0][0]['values'])
    result=[]
    for i in range(n):
        CurVal=[]
        for s in abundance.keys():
            CurVal+=[ t['values'][i]   for t in abundance[s]]
        result+=[min(CurVal)]
    return(result)

def GetProdParam(mdbfile,burn=0,nthin=1000):
    '''Get values of a,b,xmax,fmax'''
    pname=['a','b','xmax','fmax']
    ProdParam=GetParamValues(mdbfile, pname,burn=burn,nthin=nthin)
    return(ProdParam)
 
def CalcProd(B,a,b,xmax,fmax,minB=0.001,maxB=0.999):
    #Relative biomass is too small to have any productivity 
    if B<=minB:return(0.)     
    #Relative biomass is too large to have any productivity 
    if B>=maxB:return(0.)    
    
    result=fmax*\
            ((    B /   xmax ) **a) *\
            (( (1-B)/(1-xmax))**b) 
    return(result)
    
def GetMaxProd(hdf5file,burn=0,nthin=1000,minB=0.001,maxB=0.999,maxTime=9999,root='RelBmass'):
     
    MinRelAbund=GetMinRelAbund(hdf5file,burn=burn,nthin=nthin,minRelBmass=minB,maxRelBmass=maxB,maxTime=maxTime,root=root)
    ModelParam=GetParamValues(hdf5file, ['a','b','fmax','xmax','VBmass_16'],burn=burn,nthin=nthin)
    niter= len(ModelParam[0])
    
    result=[]
    
    for i in range(niter):
        a,b,fmax,xmax,VBiomass=[t[i] for t in ModelParam]
        CurMin=MinRelAbund[i]
                
        #Biomass has been depleted past the point of maximum productivity
        if CurMin<=xmax:
            result+=[fmax]            
        
        #Maximum observed productivity will coincide with minimum relative biomass
        else:
           maxProd=CalcProd(CurMin,a,b,xmax,fmax,minB=minB,maxB=maxB)
           result+=[maxProd]
    return(result)
   
    
            
   
if __name__ == "__main__":
   
   hdf5file='D:\Analyses\CukeEFAProd\MCMC\Jervis\MCMC.hdf5' 
   quantile=[.025,.5,.975]
   minRelBmass=1.e-3
   maxRelBmass=1-1e-3
   from EFAsite import Site
   nthin=10
   
   import pickle
   PickleFile='D:\Analyses\CukeEFAProd\MCMC\Jervis\Jervis.pickle'
   SiteNumber,Sites,VirginSites,CoastLength=pickle.load(open(PickleFile,"rb"))

   s=Sites[-1]
   hdata=s.SiteHarv.hdata

   x=GetMinRelAbund(hdf5file,burn=0,nthin=nthin)
   for t in x:
       print(t)
   print()
   y=GetAllRelAbund(hdf5file,burn=0,nthin=nthin)
   niter=len(y[0])
   z=[ min([t[i]  for t in y ])   for i,dummy in enumerate(y[0])]
   
   for i in range(nthin):
       print(x[i],z[i])
  
   