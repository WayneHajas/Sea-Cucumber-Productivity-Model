# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 08:58:35 2017

@author: hajasw
Last modified so I could run ProdAtBiomass
"""
import pickle

import sys
from numpy import ndarray
from scipy.stats.mstats import mquantiles

sys.path.append( '..\..\..\pyfunctions'  ) 
from EFAsite import Site
from NextRelBmass import CurRelBmass,NextRelBmass
from GetParamStats import GetParamValues2 as GetParamValues
from ReadWinBugs import  ReadWinBugs,GetParamName
from GetMinRelAbund import GetMinRelAbund


def CalcMinRelBMass(a,b,fmax,xmax,VBMass,CurSite,MaxYear=9999):
    if isinstance(CurSite,(list,ndarray)):
        SiteResult=[CalcMinRelBMass(a,b,fmax,xmax,VBMass[i],s,MaxYear=MaxYear)    for i,s in enumerate(CurSite)]
        result=min(SiteResult)
        return(result)
    
    Harvest=CurSite.SiteHarv.hdata
    Harvest=[t for t in Harvest if t[0]<=MaxYear ]
	
    #Add a zero-harvest one day after every actual-harvest
    #Results in an estimate of biomass at full impact of harvest
    Harvest2=[]
    for h in Harvest:
      Harvest2+=[h]
      #Harvest2+=[[h[0]+1/365.25,0.]]
	
	
    nharvest=len(Harvest2)
    RelBMass=[1.0]
    for i in range(1,nharvest):
      try:
        RelBMass+=[\
            NextRelBmass(OldBiomass=RelBMass[-1],OldTime=Harvest2[i-1][0],NewTime=Harvest2[i][0],a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBMass,CurHarvest=Harvest2[i-1][1])            ]
      except:
        RelBMass+=[\
            NextRelBmass(OldBiomass=RelBMass[-1],OldTime=Harvest2[i-1][0],NewTime=Harvest2[i][0],a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBMass,CurHarvest=Harvest2[i-1][1])            ]
    result=min(RelBMass)   
    return(result)

def CalcMinRelBMass2(ModelValues,Sites,MaxYear=9999):
    a,b,fmax,xmax=ModelValues[:4]
    VBMass=ModelValues[4:]
    result=CalcMinRelBMass(a,b,fmax,xmax,VBMass,Sites,MaxYear=MaxYear)
    return(result)
    
def  CalcMinRelBMass3(hdf5file,burn=1000,nthin=None,minRelBmass=1.e-3,maxRelBmass=1-1e-3,maxTime=9999,root='Rel_Abund'):
      result=GetMinRelAbund(hdf5file,burn=burn,nthin=nthin,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass,maxTime=maxTime,root=root)
      return(result)
    
def CalcProductivity(a,b,fmax,xmax,MinRelBMass,B=[(i+.5)/100  for i in range(100)]):
    try:
        result=[(bmass>=MinRelBMass) * fmax * (bmass/xmax)**a * ((1-bmass)/(1-xmax))**b   for bmass in B]
    except:
        print('calcminrelbmass 45')
        result=[(bmass>=MinRelBMass) * fmax * (bmass/xmax)**a * ((1-bmass)/(1-xmax))**b   for bmass in B]
    return(result)

def CalcProductivity2(ModelValues,Sites,B=[(i+.5)/100  for i in range(100)],TruncMinBmass=True,MaxYear=9999):
    if isinstance(ModelValues[0],(list,ndarray)):
        niter=len(ModelValues[0])
        result=[ CalcProductivity2( [t[i] for t in ModelValues]  ,Sites,B=B,MaxYear=MaxYear,TruncMinBmass=TruncMinBmass)  for i in range(len(ModelValues[0]))]
        #result=[ CalcProductivity2( t  ,Sites,B=B,MaxYear=MaxYear,TruncMinBmass=TruncMinBmass)  for t in ModelValues]
        return(result)
    
    MinRelBMass=0.0
    if TruncMinBmass:
        MinRelBMass=CalcMinRelBMass2(ModelValues,Sites,MaxYear=MaxYear)
    a,b,fmax,xmax=ModelValues[:4]
    result=CalcProductivity(a,b,fmax,xmax,MinRelBMass,B=B)
    return(result)

def QuantProductivity(ModelValues,sites,B=[(i+.5)/100  for i in range(100)],q=[.005,.025,.05,.25,.5,.75,.95,.975,.995],TruncMinBmass=True,MaxYear=9999):
    EstProd=CalcProductivity2(ModelValues,sites,B=B,TruncMinBmass=TruncMinBmass,MaxYear=MaxYear)
    nB=len(B)
    
    ByBmass=[ mquantiles([t[i] for t in EstProd ] ,q)   for i in range(nB)]
    nq=len(q)
    result={}
    result['RelBMass']=B
    for i in range(nq):
        result[q[i]]=[ BBM[i]  for BBM in ByBmass]
    return(result)

def QuantProductivity_hdf5(hdf5file,sites,B=[(i+.5)/100  for i in range(100)],q=[.005,.025,.05,.25,.5,.75,.95,.975,.995],NewName=True,TruncMinBmass=True,MaxYear=9999,burn=0,nthin=10000):
    ParamName=['a','b','fmax','xmax','VBmass_0','VBmass_2','VBmass_4','VBmass_8','VBmass_16']
    if not(NewName):
        ParamName=['a','b','fmax','xmax','VBmass_0','VBmass_1','VBmass_2','VBmass_3','VBmass_4']
    ModelValues=GetParamValues(hdf5file,ParamName,burn=burn,nthin=nthin)
    
    
    #Transform the parameter values.
    niter=len(ModelValues[0])
    #ModelValues=[ [t[i] for t in ModelValues]    for i in range(niter)]
    result=QuantProductivity(ModelValues,sites,B=B,q=q,TruncMinBmass=TruncMinBmass,MaxYear=MaxYear)
    return(result)

def QuantProductivity_hdf5_fromRel_Abund(hdf5file,sites,B=[(i+.5)/100  for i in range(100)],q=[.005,.025,.05,.25,.5,.75,.95,.975,.995],NewName=True,TruncMinBmass=True,MaxYear=9999,burn=0,nthin=10000,nb=1000):
    ParamName=['a','b','fmax','xmax']
    a,b,fmax,xmax=GetParamValues(hdf5file,ParamName,burn=burn,nthin=nthin)
    MinRelBMass3=CalcMinRelBMass3(hdf5file,burn=burn,nthin=nthin,minRelBmass=1.e-3,maxRelBmass=1-1e-3,maxTime=9999,root='Rel_Abund')
    
    if not(TruncMinBmass):
       MinRelBMass3=[0 for t in MinRelBMass3 ]
    
    niter=len(a)
    Productivity_byIter=[ CalcProductivity(a[i],b[i],fmax[i],xmax[i],MinRelBMass3[i],B=B) for i in range(niter) ]
    qprod=[ mquantiles(t,prob=q)     for t in array(Productivity_byIter).transpose()]
    
    result={}
    result['RelBMass']=B
    for i in range(len(q)):
        result[q[i]]=[t[i]  for t in  qprod   ]
    
    return(result)
def QuantProductivity_WinBugs(IndFile,OutFile,sites,B=[(i+.5)/100  for i in range(100)],q=[.005,.025,.05,.25,.5,.75,.95,.975,.995],TruncMinBmass=True,MaxYear=9999):
    ParamName=['a','b','fmax','xmax','VBmass[1]','VBmass[2]','VBmass[3]','VBmass[4]','VBmass[5]']
    ModelValuesDict=ReadWinBugs(IndFile,OutFile,ParamName)
    
    niter=len(ModelValuesDict['a'])
    ModelValues=[ [ModelValuesDict['a'][i],\
                   ModelValuesDict['b'][i],\
                   ModelValuesDict['fmax'][i],\
                   ModelValuesDict['xmax'][i],\
                   ModelValuesDict['VBmass[1]'][i],\
                   ModelValuesDict['VBmass[2]'][i],\
                   ModelValuesDict['VBmass[3]'][i],\
                   ModelValuesDict['VBmass[4]'][i],\
                   ModelValuesDict['VBmass[5]'][i]]  for i in range(niter)]

    result=QuantProductivity(ModelValues,sites,B=B,q=q,TruncMinBmass=TruncMinBmass,MaxYear=MaxYear)
    return(result)


def ProdBySim(ModelValues,sites,B=[(i+.5)/100  for i in range(100)],q=[.05+i*.1 for i in range(10)],IndexRef=-1,MaxYear=9999):
    EstProdTrunc=CalcProductivity2(ModelValues,sites,B=B,TruncMinBmass=True, MaxYear=MaxYear)
    EstProd     =CalcProductivity2(ModelValues,sites,B=B,TruncMinBmass=False,MaxYear=MaxYear)
    niter=len(EstProd)
    
    Combo=[ [EstProd[i],EstProdTrunc[i]]  for i in range(niter) ]
    if IndexRef:
       Combo.sort(key=lambda t: t[0][IndexRef])
        
    nq=len(q)
    smaller=[Combo[int(qq*niter)] for qq in q ]
    
    NoTrunc={}
    Trunc  ={}
    for i in range(nq):
        NoTrunc[q[i]]=smaller[i][0]
        Trunc[q[i]]  =smaller[i][1]
    
    
    result={}
    result['RelBMass']=B
    result['Trunc'  ]=Trunc
    result['NoTrunc']=NoTrunc
    return(result)

def ProdBySim_hdf5(hdf5file,sites,B=[(i+.5)/100  for i in range(100)],q=[.05+i*.1 for i in range(10)],IndexRef=-1,NewName=True,TruncMinBmass=True,MaxYear=9999,burn=0,nthin=10000,minb=-1):
    ParamName=['a','b','fmax','xmax','VBmass_0','VBmass_2','VBmass_4','VBmass_8','VBmass_16']
    if not(NewName):
        ParamName=['a','b','fmax','xmax','VBmass_0','VBmass_1','VBmass_2','VBmass_3','VBmass_4']
    ModelValues=GetParamValues(hdf5file,ParamName,burn=burn,nthin=nthin)
    
    #Transform the parameter values.
    niter=len(ModelValues[0])
    #ModelValues=[ [t[i] for t in ModelValues]    for i in range(niter)]
    #if minb>0:
    #    ModelValues=[t for t in ModelValues if t[1]>=minb  ]
    result=ProdBySim(ModelValues,sites,B=B,q=q,IndexRef=IndexRef,MaxYear=MaxYear)
    return(result)

def ProdBySim_WinBUGS(IndFile,OutFile,sites,B=[(i+.5)/100  for i in range(100)],q=[.05+i*.1 for i in range(10)],IndexRef=-1,NewName=True,TruncMinBmass=True,MaxYear=9999,burn=0,nthin=10000,minb=-1):
    ParamName=['a','b','fmax','xmax','VBmass[1]','VBmass[2]','VBmass[3]','VBmass[4]','VBmass[5]']
    ModelValuesDict=ReadWinBugs(IndFile,OutFile,ParamName)
    
    niter=len(ModelValuesDict['a'])
    ModelValues=[ [ModelValuesDict['a'][i],\
                   ModelValuesDict['b'][i],\
                   ModelValuesDict['fmax'][i],\
                   ModelValuesDict['xmax'][i],\
                   ModelValuesDict['VBmass[1]'][i],\
                   ModelValuesDict['VBmass[2]'][i],\
                   ModelValuesDict['VBmass[3]'][i],\
                   ModelValuesDict['VBmass[4]'][i],\
                   ModelValuesDict['VBmass[5]'][i]]  for i in range(niter)]
    if minb>0:
        ModelValues=[t for t in ModelValues if t[1]>=minb  ]
    
    result=ProdBySim(ModelValues,sites,B=B,q=q,IndexRef=IndexRef,MaxYear=MaxYear)
    return(result)


if __name__ == "__main__":  

    import pickle  
    
    
    
    hdf5file=['..\MCMC\JervisWinBugsData\MCMC.hdf5','..\MCMC\JervisWinBugsData\MCMC.20170221.hdf5','..\MCMC\JervisWinBugsData\MCMC.20170223.hdf5']
    SiteNumber,Sites,VirginSites,CoastLength=pickle.load(open("..\MCMC\Jervis\Jervis.pickle","rb"))
    
    B=[(i+.5)/10  for i in range(10)]
    q=[.025,.5,.975]    
    IndexRef=-1
    
    test=ProdBySim_hdf5(hdf5file,Sites,B=B,q=q,IndexRef=IndexRef,NewName=True,TruncMinBmass=True,MaxYear=9999,burn=0,nthin=10000)    
    
    