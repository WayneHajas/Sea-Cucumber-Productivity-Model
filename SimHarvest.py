from NextRelBmass import NextRelBmass

from numpy import ceil
from scipy.stats.mstats import mquantiles

def SimHarvest(fmax,a,b,xmax,MinRelBiomass,RelHarv, HarvestInterval,nharvest,SimInterval,minRelBmass=1.e-3,maxRelBmass=1-1e-3,OldBiomass=1):
    
    #Simulation intervals per harvest interval
    deltaRatio=   int(ceil(HarvestInterval/SimInterval) )
    #Force time-increments for simulation to be a fraction of the harvest interval   
    UseDeltaT=HarvestInterval/deltaRatio
    
    #Time after harvest to use
    TimeAfter=[.25/365.25]+[i*UseDeltaT  for i in range(1,deltaRatio+1)]
    UseMin=max([minRelBmass,MinRelBiomass])
    #Simulate to first harvest
    NewBmass=NextRelBmass(OldBiomass=OldBiomass,OldTime=0,NewTime=TimeAfter,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=1,CurHarvest=0,minRelBmass=UseMin,maxRelBmass=maxRelBmass)
    UseTime=[t for t in TimeAfter]
    
    for i in range(nharvest):
        NewTime=[UseTime[-1]+t  for t in TimeAfter]



        NewBmass+=NextRelBmass(OldBiomass=NewBmass[-1],OldTime=UseTime[-1],NewTime=NewTime,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=1,CurHarvest=RelHarv,minRelBmass=UseMin,maxRelBmass=maxRelBmass)
        UseTime+=NewTime        
        
        
    result={'biomass':NewBmass,'time':UseTime}  
    return(result)

def qSimHarvest(fmax,a,b,xmax,MinRelBiomass,RelHarv, HarvestInterval,nharvest,SimInterval,quantile=[.025,.5,.975],minRelBmass=1.e-3,maxRelBmass=1-1e-3,OldBiomass=1):
       '''Similar to SimHarvest
       fmax,a,b,xmax are all vectors
       Generates quantiles for relative biomass at time
       '''
       
       byIter=[ SimHarvest(fmax[i],a[i],b[i],xmax[i],MinRelBiomass[i],RelHarv, HarvestInterval,nharvest,SimInterval,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass,OldBiomass=OldBiomass )['biomass']    for i,f in enumerate(fmax)]
       
       i=0
       UseTime=SimHarvest(fmax[i],a[i],b[i],xmax[i],MinRelBiomass[i],RelHarv, HarvestInterval,nharvest,SimInterval,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass,OldBiomass=OldBiomass )['time'] 
       
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
    
    test=SimHarvest(fmax,a,b,xmax,MinRelBiomass,RelHarv, HarvestInterval,nharvest,SimInterval,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)
    n=len(test['time'])
    for i in range(n):
        print(test['time'][i], test['biomass'][i])
    print()
    
    fmax=[.05,.1,.5]
    a=[.5,1.,1.5]
    b=[1.5,1,.5]
    xmax=[.25,.75,.5]
    MinRelBiomass=[.05,.15,.35]
    test2=qSimHarvest(fmax,a,b,xmax,MinRelBiomass,RelHarv, HarvestInterval,nharvest,SimInterval,quantile=[.025,.5,.975],minRelBmass=1.e-3,maxRelBmass=1-1e-3)
    
    for i,t in enumerate(test2['time']):
        print(t , test2[.025][i]  , test2[.5][i]  , test2[.975][i] )