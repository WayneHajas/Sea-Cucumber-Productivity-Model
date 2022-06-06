import csv
from numpy import exp
from scipy.optimize import newton,brentq
from GetParamStats import GetParamValues2 as GetParamValues
from NextRelBmass import NextRelBmass,deriv
from GetMinRelAbund import GetVirginBiomass_WithName
from GetCoastLength import GetCoastLength
from GetMinRelAbund import GetRelAbund_WithName,GetVirginBiomass_WithName,IndexToMeanWeight,LinearBiomassDensity,FindMinima,LinearPopDensity,SpatialBiomassDensity,SpatialPopDensity

def FindDeltaB(USR,deltaT,a,b,xmax,fmax,LRP=0):
    def RootFunction(deltaB):
        B0=max([0.001001,USR-deltaB/2,LRP])
        B1=NextRelBmass(OldBiomass=B0,OldTime=0,NewTime=deltaT,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=1,CurHarvest=0)
        result=(B1-B0)-deltaB
        return(result)
        
    InitGuess=deriv(USR,a,b,xmax,fmax)[1]*deltaT
    try:
        result2=newton(RootFunction,x0=InitGuess)
    except:
        result2=brentq(RootFunction,0,deltaT*fmax)
    return(result2)

def RelativeValue(hdf5file,deltaT,burn=0,nthin=None,z=[1,2,3,4,5,6],outcsv='AllowableCatchRelBM.csv',CheckLRP=False):
    sdYear=GetParamValues(hdf5file,'sdYear',burn=burn,nthin=nthin)
    a=GetParamValues(hdf5file,'a',burn=burn,nthin=nthin)
    b=GetParamValues(hdf5file,'b',burn=burn,nthin=nthin)
    xmax=GetParamValues(hdf5file,'xmax',burn=burn,nthin=nthin)
    fmax=GetParamValues(hdf5file,'fmax',burn=burn,nthin=nthin)
    
    #If specified, incorporate the LRP into the calculations    
    LRP=[0 for t in a]
    if(CheckLRP):
        Abund=GetRelAbund_WithName(hdf5file,burn=burn,nthin=nthin,root='Rel_Abund')
        LRP=FindMinima(Abund)
    n=len(sdYear)
    colNames=['z='+str(s)   for s in z]
    csvfile=open(outcsv,'w')
    with csvfile:
        writer=csv.writer(csvfile,lineterminator='\n')
        writer.writerow(colNames)
        
        for i in range(n):
            curUSR=[exp(-t*sdYear[i])  for t in z]
            deltaB=[ FindDeltaB(t,deltaT,a[i],b[i],xmax[i],fmax[i],LRP=LRP[i])      for t in curUSR]  
            writer.writerow(deltaB) 
            if not(i %1000):print(n,i)

def LinearBiomassDensity(hdf5file,deltaT,survey,burn=0,nthin=None,z=[1,2,3,4,5,6],outcsv='AllowableCatchLinearBiomassDensity.csv',CheckLRP=False):
    sdSite=GetParamValues(hdf5file,'sdSite',burn=burn,nthin=nthin)
    sdTransect=GetParamValues(hdf5file,'sdTransect',burn=burn,nthin=nthin)
    sdYear=GetParamValues(hdf5file,'sdYear',burn=burn,nthin=nthin)
    a=GetParamValues(hdf5file,'a',burn=burn,nthin=nthin)
    b=GetParamValues(hdf5file,'b',burn=burn,nthin=nthin)
    xmax=GetParamValues(hdf5file,'xmax',burn=burn,nthin=nthin)
    fmax=GetParamValues(hdf5file,'fmax',burn=burn,nthin=nthin)
    VBiomass=GetVirginBiomass_WithName(hdf5file,burn=burn,nthin=nthin,root='VBmass')
    #If specified, incorporate the LRP into the calculations    
    LRP=[0 for t in a]
    if(CheckLRP):
        Abund=GetRelAbund_WithName(hdf5file,burn=burn,nthin=nthin,root='Rel_Abund')
        LRP=FindMinima(Abund)
    
    CoastLength=GetCoastLength(survey,mdbfile='S:\\analyses\\EFA_Productivity_Model.20180929\\data\\CoastLength.mdb')
    
    n=len(VBiomass[0])
    site=list(VBiomass.keys())
    TotalCoastLength=sum([CoastLength[s]   for s in site])
    TotalVB=[  sum([VBiomass[s][i] for s in site ])     for i in range(n)]
    
    
    colNames=['z='+str(t)   for t in z]
    csvfile=open(outcsv,'w')
    with csvfile:
        writer=csv.writer(csvfile,lineterminator='\n')
        writer.writerow(colNames)
    
        for i in range(n):
            curUSR=[exp(-t*sdYear[i])  for t in z]
            deltaB=[ FindDeltaB(t,deltaT,a[i],b[i],xmax[i],fmax[i],LRP=LRP[i])      for t in curUSR]  
            curAC=[TotalVB[i]/TotalCoastLength*t*exp( (sdTransect[i]**2)/2+(sdSite[i]**2)/2+(sdYear[i]**2)/2)  for t in deltaB]
            writer.writerow(curAC)
            if not(i %1000):print(n,i)

def SpatialBiomassDensity(hdf5file,deltaT,burn=0,nthin=None,z=[1,2,3,4,5,6],outcsv='AllowableCatchSpatialBiomassDensity.csv',CheckLRP=False):
    sdSite=GetParamValues(hdf5file,'sdSite',burn=burn,nthin=nthin)
    sdTransect=GetParamValues(hdf5file,'sdTransect',burn=burn,nthin=nthin)
    sdYear=GetParamValues(hdf5file,'sdYear',burn=burn,nthin=nthin)
    a=GetParamValues(hdf5file,'a',burn=burn,nthin=nthin)
    b=GetParamValues(hdf5file,'b',burn=burn,nthin=nthin)
    xmax=GetParamValues(hdf5file,'xmax',burn=burn,nthin=nthin)
    fmax=GetParamValues(hdf5file,'fmax',burn=burn,nthin=nthin)
    lnGrandMean=GetParamValues(hdf5file, 'lnGrandMean',burn=burn,nthin=nthin) 
    #If specified, incorporate the LRP into the calculations    
    LRP=[0 for t in a]
    if(CheckLRP):
        Abund=GetRelAbund_WithName(hdf5file,burn=burn,nthin=nthin,root='Rel_Abund')
        LRP=FindMinima(Abund)   
    
    n=len(lnGrandMean)
   
    
    colNames=['z='+str(t)   for t in z]
    csvfile=open(outcsv,'w')
    with csvfile:
        writer=csv.writer(csvfile,lineterminator='\n')
        writer.writerow(colNames)
    
        for i in range(n):
            curUSR=[exp(-t*sdYear[i])  for t in z]
            deltaB=[ FindDeltaB(t,deltaT,a[i],b[i],xmax[i],fmax[i],LRP=LRP[i])      for t in curUSR]  
            curAC=[t*exp(lnGrandMean[i]+ (sdTransect[i]**2)/2+(sdSite[i]**2)/2+(sdYear[i]**2)/2)  for t in deltaB]
            writer.writerow(curAC)
            if not(i %1000):print(n,i)

if __name__ == "__main__":  
        
    class Tolmieclass():
        def __init__(self):
            self.hdf5file=[              'S:\\analyses\\EFA_Productivity_Model.20180929\\Tolmie\\NewModel\\seed.20180824.hdf5',\
    			'S:\\analyses\\EFA_Productivity_Model.20180929\\Tolmie\\NewModel\\seed.20180825.hdf5',\
    			'S:\\analyses\\EFA_Productivity_Model.20180929\\Tolmie\\NewModel\\seed.20180826.hdf5',\
    			'S:\\analyses\\EFA_Productivity_Model.20180929\\Tolmie\\NewModel\\seed.20180827.hdf5']
            self.burn=0
            self.nthin=None
    Tolmiehdf5file=Tolmieclass()
    survey='Tolmie Channel'
    
    deltaT=1
    test1= RelativeValue(Tolmiehdf5file.hdf5file,deltaT,burn=Tolmiehdf5file.burn,nthin=Tolmiehdf5file.nthin,z=[1,2,3,4,5,6],outcsv='AllowableCatchRelBM.csv')  
    test3= SpatialBiomassDensity(Tolmiehdf5file.hdf5file,deltaT,burn=Tolmiehdf5file.burn,nthin=Tolmiehdf5file.nthin,z=[1,2,3,4,5,6],outcsv='AllowableCatchLinearBiomassDensity.csv')     
    test2= LinearBiomassDensity(Tolmiehdf5file.hdf5file,deltaT,survey,burn=Tolmiehdf5file.burn,nthin=Tolmiehdf5file.nthin,z=[1,2,3,4,5,6],outcsv='AllowableCatchLinearBiomassDensity.csv')     
