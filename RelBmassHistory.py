
from numpy import ndarray
from scipy.stats.mstats import mquantiles
from ToDecimalYear import ToDecimalYear
from NextRelBmass import NextRelBmass
from GetParamStats import GetParamValues2 as GetParamValues
from ReadWinBugs import  ReadWinBugs,GetParamName

def RelBmassHistory(OldBiomass,NewTime,a,b,xmax,fmax,VBiomass,HarvestHist,minRelBmass=1.e-3,maxRelBmass=1-1e-3, TrustNewTimeForm=False):
    '''
    Generate a probabilistic history of relative biomass    
    
    * OldBiomass is the relative-biomass at the beginnin of the time-interval
    * OldTime is the decimal year at the beginning of the time-interval
    * NewTime is the decimal year at the end       of the time-interval
    * a,b,xmax,fmax are the parameters of the WCH productivity model 
    * VBiomass is the virgin biomass
    * HarvestHist is a list with the same structure as  Harvest.hdata   
        [  [decimal year, harvest amount],
           [decimal year, harvest amount],
           ...,
           [decimal year, harvest amount]]

    
    * minRelBmass and maxRelBmass are the minimum and maximum values of relative-biomass that will be considered.
    * TrustNewTimeForm indicates that NewTime values
        - are decimal
        - include dates of harvest
        - are sorted
    '''
    #Multiple iterations from MCMC
    if isinstance(a,(list,ndarray)):
        niter=len(a)
        UseTime=NewTime
        ntime=len(UseTime)
        
        if not(TrustNewTimeForm):
            htime=[ToDecimalYear(t[0])    for t in HarvestHist if (t[0]> ToDecimalYear(NewTime[0])) and  (t[0]< ToDecimalYear(NewTime[-1])) ]
            ptime=[t +1./365.25 for t in htime] #day after harvest
            stime=[ToDecimalYear(t) for t in NewTime]
            UseTime1=htime+ptime+stime
            UseTime1=set(UseTime1)
            UseTime2=[t for t in UseTime1]
            UseTime2.sort()
            ntime=len(UseTime2)
        niter=len(a)
        ByIter=[  RelBmassHistory(OldBiomass,UseTime2,a[i],b[i],xmax[i],fmax[i],VBiomass[i],HarvestHist,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass, TrustNewTimeForm=True)  for i in range(niter)]
        ByTime=[ [t['RelBmass'][j] for t in ByIter] for j in range(ntime) ]
        
        result={}        
        result['time']=UseTime2
        result['RelBmass']=ByTime
        return(result)
    #Single iteration from the MCMC
        
    UseTime2=NewTime
    ntime=len(UseTime2)
    
    if not(TrustNewTimeForm):
        htime=[ToDecimalYear(t[0])    for t in HarvestHist if (t[0]> ToDecimalYear(NewTime[0])) and  (t[0]< ToDecimalYear(NewTime[-1])) ]
        ptime=[t +1./365.25 for t in htime] #day after harvest
        stime=[ToDecimalYear(t) for t in NewTime]
        UseTime1=set(stime+ptime+htime  )
        UseTime2=[t for t in UseTime1]
        UseTime2.sort()
    ntime=len(UseTime2)
    result={}
    result['RelBmass']=[OldBiomass]
    OldTime=UseTime2[0]
    for t in UseTime2[1:]:
        CurHarvest=GetHarvest(OldTime,HarvestHist)
        result['RelBmass']+=[NextRelBmass(OldBiomass=result['RelBmass'][-1],OldTime=OldTime,NewTime=t,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBiomass,CurHarvest=CurHarvest,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)]
        OldTime=t
    result['time']=UseTime2
    return(result)
 
def QuantRelBmassHistory_hdf5(hdf5file,site,VBmass_name='VBmass_16',NewTime=[1995+i for i in range(16)],q=[.005,.025,.05,.25,.5,.75,.95,.975,.9995],minRelBmass=1.e-3,maxRelBmass=1-1e-3,burn=0,nthin=10000):
    ParamName=['a','b','fmax','xmax',VBmass_name]
    a     =GetParamValues(hdf5file,      'a',burn=burn,nthin=nthin)
    b     =GetParamValues(hdf5file,      'b',burn=burn,nthin=nthin)
    fmax  =GetParamValues(hdf5file,      'fmax',burn=burn,nthin=nthin)
    xmax  =GetParamValues(hdf5file,      'xmax',burn=burn,nthin=nthin)
    VBmass=GetParamValues(hdf5file,VBmass_name ,burn=burn,nthin=nthin)
    HarvestHist=site.SiteHarv.hdata

    PredRelBmass=RelBmassHistory(1,NewTime,a,b,xmax,fmax,VBmass,HarvestHist,minRelBmass=1.e-3,maxRelBmass=1-1e-3, TrustNewTimeForm=False)
    UseTime=PredRelBmass['time']
    ByTime=[mquantiles(t,q) for t in PredRelBmass['RelBmass']  ]
    nq=len(q)
    
    result={}
    result['time']=UseTime
    for i in range(nq):
        result[q[i]]=[t[i] for t in ByTime]
    return(result)
 
def QuantRelBmassHistory_WinBugs(IndFile,OutFile,site,VBmass_name='VBmass[5]',NewTime=[1995+i for i in range(16)],q=[.005,.025,.05,.25,.5,.75,.95,.975,.9995],minRelBmass=1.e-3,maxRelBmass=1-1e-3,burn=0,nthin=10000):
    a     =ReadWinBugs(IndFile,OutFile,      'a',burn=burn,nthin=nthin)
    b     =ReadWinBugs(IndFile,OutFile,      'b',burn=burn,nthin=nthin)
    fmax  =ReadWinBugs(IndFile,OutFile,      'fmax',burn=burn,nthin=nthin)
    xmax  =ReadWinBugs(IndFile,OutFile,      'xmax',burn=burn,nthin=nthin)
    VBmass=ReadWinBugs(IndFile,OutFile,VBmass_name ,burn=burn,nthin=nthin)
    HarvestHist=site.SiteHarv.hdata

    PredRelBmass=RelBmassHistory(1,NewTime,a,b,xmax,fmax,VBmass,HarvestHist,minRelBmass=1.e-3,maxRelBmass=1-1e-3, TrustNewTimeForm=False)
    UseTime=PredRelBmass['time']
    ByTime=[mquantiles(t,q) for t in PredRelBmass['RelBmass']  ]
    nq=len(q)
    
    result={}
    result['time']=UseTime
    for i in range(nq):
        result[q[i]]=[t[i] for t in ByTime]
    return(result)
   
def GetHarvest(tdate, HarvestHist):
    occur=[t for t in HarvestHist if t[0]==tdate]
    if occur:
        return(occur[0][1])
    return(0.0)
    
    

if __name__ == "__main__":  

    import pickle  
        
    
    a,b,fmax,xmax=0.386202457391,0.749317072977,0.0855760509984,0.373301217716
    VBMass=27013.6482306,12641.6722443,30496.9893597,25826.9854016,19335.8963094
    SiteNumber,Sites,VirginSites,CoastLength=pickle.load(open("..\MCMC\Jervis\Jervis.pickle","rb"))
    
        
    hdf5file=['..\MCMC\JervisWinBugsData\MCMC.hdf5','..\MCMC\JervisWinBugsData\MCMC.20170221.hdf5','..\MCMC\JervisWinBugsData\MCMC.20170223.hdf5']

    CurSite=Sites[-1]
    HarvestHist=CurSite.SiteHarv.hdata
    
    OldBiomass=1.0
    VBiomass=VBMass[-1]
    NewTime=[1996+i for i in range(20)]    
    x01=RelBmassHistory(OldBiomass,NewTime,a,b,xmax,fmax,VBiomass,HarvestHist,minRelBmass=1.e-3,maxRelBmass=1-1e-3, TrustNewTimeForm=False)
    
    
    x02=QuantRelBmassHistory_hdf5(hdf5file,CurSite,VBmass_name='VBmass_16',NewTime=[1995+i for i in range(16)],q=[.025,.5000,0.975],minRelBmass=1.e-3,maxRelBmass=1-1e-3,burn=0,nthin=10)


  
    IndFile=['f:\\Archive\\s-drive\\analyses\\CukeExpHarvRecov\\Jervis.beta\\chain1.ind' ,'f:\\Archive\\s-drive\\analyses\\CukeExpHarvRecov\\Jervis.beta\\chain2.ind']       
    OutFile=['f:\\Archive\\s-drive\\analyses\\CukeExpHarvRecov\\Jervis.beta\\chain1.out' ,'f:\\Archive\\s-drive\\analyses\\CukeExpHarvRecov\\Jervis.beta\\chain2.out']     
    x03=QuantRelBmassHistory_WinBugs(IndFile,OutFile,CurSite,VBmass_name='VBmass[5]',NewTime=[1995+i for i in range(16)])
   