'''
Utilities for finding the maximum sustainable harvest amount

20210210

'''

from scipy.optimize import minimize_scalar
from scipy.stats.mstats import mquantiles
from numpy import log,exp,average
from numpy.random import normal,seed
import csv

from HarvestAmountUtility import Prod,EstPreB
from NextRelBmass import NextRelBmass
from GetParamStats import GetParamValues2 as GetParamValues
from GetMinRelAbund import GetMinRelAbund



### Harvest amount is a fraction of virgin biomass
def minFunc(Bpost,HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund):
    '''
    Evaluate Bpost according to sustainable harvest amount
    '''

    #Calculate the pre-harvest relative biomass
    Bpre=EstPreB(Bpost,HarvestInterval,a,b,xmax,fmax)

    RelHarv=Bpre-Bpost
    #Take negative because scipy utilities want to find the minimum
    return(-RelHarv)

def MaxRelHarvest(HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund):
    '''
    Find the post-harvest relative biomass corresponding to the maximum sustainable harvest amount
    '''
    
    #If there are lists, do operation for each iteration in markov chain 
    if isinstance(fmax,list):
        n=len(fmax)
        result=[  MaxRelHarvest(HarvestInterval,SimInterval,fmax[i],a[i],b[i],xmax[i],MinRelAbund[i],MaxRelAbund[i])   for i in range(n)]
        return(result)
    
    #One iteration of markov chain
   
    #Check for trivial results 
    if xmax<(MinRelAbund):
            return(MinRelAbund)
    
    args=(HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund)
    Bpost=minimize_scalar(minFunc,args=args,method='Brent',bracket=(.001,0.999)).x
    
    
    #Special Cases
    if not(Bpost):
        return(0)
    if Bpost<MinRelAbund:
        return(MinRelAbund)
    if Bpost>0.999:
        return(0.999)
    # Calculation is OK:
    return(Bpost)


def WriteUSR3(hdf5file,burn,nthin,HarvestInterval,outfilename):
         
    outfile=csv.writer(open(outfilename,'w'),lineterminator='\n')
    
    
    #Get Values
    a=GetParamValues(hdf5file, 'a',burn=burn,nthin=nthin)
    b=GetParamValues(hdf5file, 'b',burn=burn,nthin=nthin)
    fmax=GetParamValues(hdf5file, 'fmax',burn=burn,nthin=nthin)
    xmax=GetParamValues(hdf5file, 'xmax',burn=burn,nthin=nthin)
    MinRelAbund=GetMinRelAbund(hdf5file,burn=burn,nthin=nthin,root='Rel_Abund')
    MaxRelAbund=[0.999 for t in MinRelAbund]
    
    #Write some information so we know exactly how the values were generated
    outfile.writerow(hdf5file)
    outfile.writerow([burn,nthin])
    outfile.writerow(HarvestInterval)
    
    n =len(a)
    for i in range(n):
        #print(i)
        LRP21=[MaxRelHarvest(hi,hi,fmax[i],a[i],b[i],xmax[i],MinRelAbund[i],MaxRelAbund[i]) for hi in HarvestInterval]
        outfile.writerow(LRP21)
        if (LRP21[0]!=LRP21[-1]):
            print(fmax[i],a[i],b[i],xmax[i],MinRelAbund[i],MaxRelAbund[i])
            print(LRP21)
            print()
    
    del outfile
    
def ReadUSR3(hdf5file,burn,nthin,HarvestInterval,infilename):
         
    with open(infilename,newline='\n') as infile:
        csvreader=csv.reader(infile,delimiter=',',quotechar="'")
        
        #Make sure hdf5files are the same
        testhdf5file=csvreader.__next__()
        if testhdf5file!=hdf5file:
            return('hdf5 files do not match')
        
        #Make sure burn and nthin are same as used to create file
        testburnnthin=[int(t) for t in csvreader.__next__()]
        if testburnnthin!=[burn,nthin]:
            return('burn or nthin does not match')
        
        #Make sure harvest interval is included in the file
        testHarvestInterval=[int(t) for t in csvreader.__next__()]
        if HarvestInterval not in testHarvestInterval:
            return('HarvestInterval not included in file')
        
        #Get values for the harvest intervals
        index=testHarvestInterval.index(HarvestInterval)
        result=[ float(row[index])  for row in csvreader]
        return(result)
        
        
    
if __name__ == "__main__":
    RelHarv,HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund=0.05,1,1,.1,.5,.5,.5,.025,.999
    test1=MaxRelHarvest(HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund)
    print(test1)
    hdf5file=['c:\\analyses\\EFA_Productivity_Model.20180929\\Zeballos\\NewModel_wideSdYear_WideSiteArea.2007\\seed.20180834.hdf5',\
              'c:\\analyses\\EFA_Productivity_Model.20180929\\Zeballos\\NewModel_wideSdYear_WideSiteArea.2007\\seed.20180835.hdf5',\
              'c:\\analyses\\EFA_Productivity_Model.20180929\\Zeballos\\NewModel_wideSdYear_WideSiteArea.2007\\seed.20180836.hdf5',\
              'c:\\analyses\\EFA_Productivity_Model.20180929\\Zeballos\\NewModel_wideSdYear_WideSiteArea.2007\\seed.20180838.hdf5',\
              'c:\\analyses\\EFA_Productivity_Model.20180929\\Zeballos\\NewModel_wideSdYear_WideSiteArea.2007\\seed.20180839.hdf5']
    
    hdf5file=[hdf5file[0]]
    nthin=10
    burn,nthin,HarvestInterval,outfilename=0, 1000, [1, 2, 3, 4, 5], 'c:\\scratch\\USR.2.1.csv'
    #WriteUSR21(hdf5file,burn,nthin,HarvestInterval,outfilename)
    
    HarvestInterval=3
    test2=ReadUSR3(hdf5file,burn,nthin,HarvestInterval,outfilename)