from NextRelBmass import deriv
from GetNames import *
from GetParamStats import GetParamValues2 as GetParamValues

from scipy.stats.mstats import mquantiles
from numpy import ndarray

def GetRel_Abund_Names(mdbfile,EightSixteenOnly=False):
    pname=GetNames(mdbfile)
    pname=[t for t in pname if t[:9]=='RelBmass_']
    if EightSixteenOnly:
         pname=[t for t in pname if (t[:10]=='RelBmass_3') or (t[:11]=='RelBmass_4') ]
    return(pname)
    
def GetProdFuncParam(mdbfile,burn=0,nthin=1000,EightSixteenOnly=False):
    
    #Parameter names associated with relative biomass    
    rname=GetRel_Abund_Names(mdbfile,EightSixteenOnly=EightSixteenOnly)
    allname=["fmax","a","b","xmax"]+rname
    result=GetParamValues(mdbfile, allname,burn=burn,nthin=nthin)
    return(result)    
 
def CalcQuantileProdFunc(CurBmass,ParamVal,quantile=[.025,.5,.975],minRelBmass=1.e-3,maxRelBmass=1-1e-3,ApplyTrunc=False):

    #Number of parameters and numberiterations in MCMC
    nparam=len(ParamVal)
    niter=len(ParamVal[0])
    if isinstance(CurBmass,(float)):
        CalcVal=[]
        for i  in range(niter):
          fmax,a,b,xmax=ParamVal[0][i],  ParamVal[1][i],ParamVal[2][i],ParamVal[3][i]
          
          #Productivity as estimated from model
          CalcVal+=[deriv(CurBmass,a,b,xmax,fmax,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)[1]]
          
          #Check to see if range of CurBmass is beyond range observed in iteration of MCMC
          if ApplyTrunc:          
            #Relative biomass values that occured in the MCMC
            ObsBmass=[ ParamVal[j][i]    for j in range(4,nparam)]
            lowObsBmass=min(ObsBmass)#Lowest relative biomass that occurs in iteration of MCMC
            lowObsBmass=max([lowObsBmass,minRelBmass])#Lowest relative biomass with non-zero productivity
            
            for i, dummy in enumerate(CurBmass):
                if CurBmass[i]<lowObsBmass:
                    CalcVal[i]=0.
         
        result=mquantiles(CalcVal,prob=quantile)
        return (result)

    #CurBmass is a list or array 
    #For effiency, go throught the MCMC once and then transform results
    CalcVal=[]
    for i  in range(niter):
      fmax,a,b,xmax=ParamVal[0][i],  ParamVal[1][i],ParamVal[2][i],ParamVal[3][i]
          
      #Relative biomass values that occured in the MCMC
      ObsBmass=[ ParamVal[j][i]    for j in range(4,nparam)]
      lowObsBmass=min(ObsBmass)#Lowest relative biomass that occurs in iteration of MCMC
      lowObsBmass=max([lowObsBmass,minRelBmass])#Lowest relative biomass with non-zero productivity
      
      if ApplyTrunc:
         CalcVal+=[[deriv(cb,a,b,xmax,fmax,minRelBmass=lowObsBmass,maxRelBmass=maxRelBmass)[1]   for cb in CurBmass]]
         
      else:
         CalcVal+=[[deriv(cb,a,b,xmax,fmax,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)[1]   for cb in CurBmass]]
         
    QuantAtCB=[]
    for i,cb in enumerate( CurBmass):
      Curval=[t[i]   for t in CalcVal ]
      QuantAtCB+=[mquantiles(Curval,prob=quantile)]
        
    
    #More convenient to have results in plotable lists
    result={}
    result['Relative Biomass']=CurBmass
    for i,q in enumerate(quantile):
        result[q]=[t[i] for t in QuantAtCB]
    
    return(result)
       
   
if __name__ == "__main__":
   
   mdbfile='D:\Analyses\CukeEFAProd\MCMC\Jervis\MCMC.hdf5' 
   quantile=[.025,.5,.975]
   minRelBmass=1.e-3
   maxRelBmass=1-1e-3

   x=GetRel_Abund_Names(mdbfile)
   for t in x[:10]:
       print(t)
   ParamVal=GetProdFuncParam(mdbfile,burn=0,nthin=10)   
   x=CalcQuantileProdFunc(0.5,ParamVal,quantile=quantile,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)
   print(x)
   print()
   RB=[(i+.5)/10.  for i in range(10)]
   x=CalcQuantileProdFunc(RB,ParamVal,quantile=quantile,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)
   print(x)
    