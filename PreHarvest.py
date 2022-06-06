from NextRelBmass import CurRelBmass

from numpy import ceil
from scipy.stats.mstats import mquantiles

def SimHarvest(fmax,a,b,xmax,MinRelBiomass,RelHarv, HarvestInterval,nharvest,minRelBmass=1.e-3,maxRelBmass=1-1e-3,OldBiomass=1):
    
    
    
    #Time before harvest to use
    HarvestHist=[[i*HarvestInterval,RelHarv]  for i in range(1,nharvest+1)]
    UseTime=[0]+[-0.25/365.25+t[0] for t in HarvestHist]
    
    NewBmass=[CurRelBmass(NewTime=t,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=1,HarvestHist=HarvestHist,minRelBmass=1.e-3,maxRelBmass=1-1e-3) for t in UseTime]
         
    result={'biomass':NewBmass,'time':UseTime}  
    return(result)

def qSimHarvest(fmax,a,b,xmax,MinRelBiomass,RelHarv, HarvestInterval,nharvest,quantile=[.025,.5,.975],minRelBmass=1.e-3,maxRelBmass=1-1e-3,OldBiomass=1):
       '''Similar to SimHarvest
       fmax,a,b,xmax are all vectors
       Generates quantiles for relative biomass at time
       '''
       
       byIter=[ SimHarvest(fmax[i],a[i],b[i],xmax[i],MinRelBiomass[i],RelHarv, HarvestInterval,nharvest,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass,OldBiomass=OldBiomass )['biomass']    for i,f in enumerate(fmax)]
       
       i=0
       UseTime=SimHarvest(fmax[i],a[i],b[i],xmax[i],MinRelBiomass[i],RelHarv, HarvestInterval,nharvest,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass,OldBiomass=OldBiomass )['time'] 
       
       ntime=len(byIter[0])
       QuantAtTime=[]
       for i in range(ntime):
           curBmass=[t[i] for t in byIter]
           QuantAtTime+=[mquantiles(curBmass,prob=quantile)]
           
    
       #More convenient to have results in plotable lists
       result={}
       result['time']=UseTime
       for i,q in enumerate(quantile):
        result[q]=[t[i] for t in QuantAtTime]

       return(result)           

if __name__ == "__main__":
    
    fmax,a,b,xmax,MinRelBiomass=.1,1,1,.5,.25
    RelHarv=0.2
    HarvestInterval,nharvest,SimInterval=3,5,.8
    minRelBmass=1.e-3
    maxRelBmass=1-1e-3
    nsim=1000
    
    test=SimHarvest(fmax,a,b,xmax,MinRelBiomass,RelHarv, HarvestInterval,nharvest,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)
    fmax,a,b,xmax,MinRelBiomass=nsim*[.1],nsim*[1],nsim*[1],nsim*[.5],nsim*[.25]
    test2=qSimHarvest(fmax,a,b,xmax,MinRelBiomass,RelHarv, HarvestInterval,nharvest,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)

