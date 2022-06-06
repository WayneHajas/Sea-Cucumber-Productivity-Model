import csv
from numpy import exp
from GetParamStats import GetParamValues2 as GetParamValues
from GetMinRelAbund import GetRelAbund_WithName,IndexToMeanWeight,GetVirginBiomass_WithName,GetMinRelAbund
from GetCoastLength import GetCoastLength

def CalcUSR(hdf5file,burn=0,nthin=None,z=[1,2,3,4,5,6],outcsv='USRbySdYear.csv',MultipleMinRelAbund=None,root='Rel_Abund'):
    sdYear=GetParamValues(hdf5file,'sdYear',burn=burn,nthin=nthin)
    MinRelAbund=GetMinRelAbund(hdf5file,burn=burn,nthin=nthin,root=root)
    n=len(sdYear)
    if(MultipleMinRelAbund):
        LRP=[MultipleMinRelAbund*t for t in MinRelAbund]
    colNames=['z='+str(t)   for t in z]
    csvfile=open(outcsv,'w')
    with csvfile:
        writer=csv.writer(csvfile,lineterminator='\n')
        writer.writerow(colNames)
    
        for i,t in enumerate(sdYear):
            curUSR=[exp(-s*t)  for s in z]
            if(MultipleMinRelAbund):
                curUSR=[max([s,LRP[i]])    for s in curUSR]
            writer.writerow(curUSR)
def ReadUSR(incsv='USRbySdYear.csv'):
    with open(incsv,newline='\n') as infile:
        curreader=csv.reader(infile)
        ColumnNames=next(curreader)
        result=[[float(s) for s in t] for t in curreader]
    return(result)    
    
           
def LinearBiomassDensity(hdf5file,survey,burn=0,nthin=None,incsv='USRbySdYear.csv',outcsv='USRLinearBiomassDensity.csv',MultipleMinRelAbund=None,root='Rel_Abund'):
    CUSR=ReadUSR(incsv=incsv)
    sdSite=GetParamValues(hdf5file,'sdSite',burn=burn,nthin=nthin)
    sdTransect=GetParamValues(hdf5file,'sdTransect',burn=burn,nthin=nthin)
    sdYear=GetParamValues(hdf5file,'sdYear',burn=burn,nthin=nthin)
    VBiomass=GetVirginBiomass_WithName(hdf5file,burn=burn,nthin=nthin,root='VBmass')
    if len(CUSR)!=len(sdSite):
        dummy=1/0
    
    CoastLength=GetCoastLength(survey,mdbfile='S:\\analyses\\EFA_Productivity_Model.20180929\\data\\CoastLength.mdb')
    
    n=len(VBiomass[0])
    site=list(VBiomass.keys())
    TotalCoastLength=sum([CoastLength[s]   for s in site])
    TotalVB=[  sum([VBiomass[s][i] for s in site ])     for i in range(n)]
    
    
    colNames=['z='+str(i+1)   for i in range(len(CUSR[0]))]
    csvfile=open(outcsv,'w')
    with csvfile:
        writer=csv.writer(csvfile,lineterminator='\n')
        writer.writerow(colNames)
    
        for i in range(n):
            bRBM=CUSR[i]
            curUSR=[TotalVB[i]/TotalCoastLength*t*exp( (sdTransect[i]**2)/2+(sdSite[i]**2)/2+(sdYear[i]**2)/2)  for t in bRBM]
            writer.writerow(curUSR)
    
   
            
def SpatialBiomassDensity(hdf5file,burn=0,nthin=None,incsv='USRbySdYear.csv',outcsv='USRSpatialBiomassDensity.csv',MultipleMinRelAbund=None,root='Rel_Abund'):
    byRelBM=ReadUSR(incsv=incsv)
    sdSite=GetParamValues(hdf5file,'sdSite',burn=burn,nthin=nthin)
    sdTransect=GetParamValues(hdf5file,'sdTransect',burn=burn,nthin=nthin)
    sdYear=GetParamValues(hdf5file,'sdYear',burn=burn,nthin=nthin)
    lnGrandMean=GetParamValues(hdf5file,'lnGrandMean',burn=burn,nthin=nthin)
    if len(byRelBM)!=len(sdSite):
        dummy=1/0
    n=len(lnGrandMean)
    
    colNames=['z='+str(i+1)   for i in range(len(byRelBM[0]))]
    csvfile=open(outcsv,'w')
    with csvfile:
        writer=csv.writer(csvfile,lineterminator='\n')
        writer.writerow(colNames)
    
        for i in range(n):
            bRBM=byRelBM[i]
            curUSR=[t*exp(lnGrandMean[i]+ (sdTransect[i]**2)/2+(sdSite[i]**2)/2+(sdYear[i]**2)/2)  for t in bRBM]
            writer.writerow(curUSR)
    
   