# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 14:53:07 2022

@author: HajasW
"""
from numpy import log,exp,average,ndarray
from scipy.stats.mstats import mquantiles
import sys
#sys.path.append('C:/Analyses/EFA_Productivity_Model.20180929/pyfunctions')
#from HarvestAmountUtility import EstPreB as NoTruncEstPreB


def EstPreB(Bpost,t,a,b,xmax,fmax,MinRelAbund):
    if isinstance(a, (list,ndarray)):
        result=[  EstPreB(Bpost,t,a[i],b[i],xmax[i],fmax[i],MinRelAbund[i]) for i in range(len(a))]
        return(result)
    
    if Bpost<=MinRelAbund:
        return(0)
    result=NoTruncEstPreB(Bpost,t,a,b,xmax,fmax)
    return(result)
    
def EstFracVirgin(Bpost,t,a,b,xmax,fmax,MinRelAbund):
    if isinstance(a, (list,ndarray)):
        result=[  EstFracVirgin(Bpost,t,a[i],b[i],xmax[i],fmax[i],MinRelAbund[i]) for i in range(len(a))]
        return(result)
    Bpre= EstPreB(Bpost,t,a,b,xmax,fmax,MinRelAbund)
    if Bpre<=Bpost:
        return(0)
    result=Bpre-Bpost
    return(result)
    
def EstFracBpre(Bpost,t,a,b,xmax,fmax,MinRelAbund):
    if isinstance(a, (list,ndarray)):
        result=[  EstFracBpre(Bpost,t,a[i],b[i],xmax[i],fmax[i],MinRelAbund[i]) for i in range(len(a))]
        return(result)
    Bpre= EstPreB(Bpost,t,a,b,xmax,fmax,MinRelAbund)
    if Bpre<=Bpost:
        return(0)
    result=(Bpre-Bpost)/Bpre
    return(result)

def qfunc(func,Bpost,t,a,b,xmax,fmax,MinRelAbund,prob=[.025,.500,.975] ):
    if isinstance(Bpost, (list,ndarray)):
        result=[ qfunc(func,r,t,a,b,xmax,fmax,MinRelAbund,prob=prob ) for r in Bpost]
        return(result)
    values=func(Bpost,t,a,b,xmax,fmax,MinRelAbund )
    result=list(mquantiles(values,prob=prob))
    return(result)
    
def CurvesFunc(func,Bpost,t,a,b,xmax,fmax,MinRelAbund,prob=[.025,.500,.975] ):
    qvals=qfunc(func,Bpost,t,a,b,xmax,fmax,MinRelAbund,prob=prob)
    result={}
    result['Bpost']=Bpost
    for i in range(len(prob)):
        result[prob[i]]=[t[i]  for t in qvals]
    return(result)
    

    
if __name__ == "__main__":
    
    from GetMinRelAbund import GetMinRelAbund
    from GetParamStats import GetParamValues2 as GetParamValues    
    sys.path.append('../')
    from hdf5file import hdf5file,burn,nthin,PickleFile
    nthin=1000
    hdf5file='C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea/seed.20180824.hdf5'
    nthin=10
    
    a=GetParamValues(hdf5file, 'a',burn=burn,nthin=nthin)
    b=GetParamValues(hdf5file, 'b',burn=burn,nthin=nthin)
    fmax=GetParamValues(hdf5file, 'fmax',burn=burn,nthin=nthin)
    xmax=GetParamValues(hdf5file, 'xmax',burn=burn,nthin=nthin)
    MinRelAbund=GetMinRelAbund(hdf5file,burn=burn,nthin=nthin,minRelBmass=1.e-3,maxRelBmass=1-1e-3,maxTime=9999,root='Rel_Abund')
    

    print('finished reading data')
    Bpost=[(i+.5)/10 for i in range(10)]
    test= CurvesFunc(EstFracVirgin,Bpost,1,a,b,xmax,fmax,MinRelAbund,prob=[.025,.500,.975] )
    Bpost,t,a,b,xmax,fmax,MinRelAbund=.5,1,a[0],b[0],xmax[0],fmax[0],MinRelAbund[0]
