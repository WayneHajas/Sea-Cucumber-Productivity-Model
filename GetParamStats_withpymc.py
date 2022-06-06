
import os
import sys
import  numpy as np
#from pymc import  *
#from pylab import *

#2020-12-04
   # problems with 
    #from utils import hpd, quantiles
    #ImportError: cannot import name 'iterable' from 'matplotlib.cbook' (C:\Users\hajasw\Anaconda3\envs\python38\lib\site-packages\matplotlib\cbook\__init__.py)
    
    #I am going to try removing the universal input and just get the functions I need as I find out I need them.
#20210216
    #Remove pymc libraries
    #Sometimes there are installation issues (I think).
    #Problems accessing pymc libaries.
    #I am just going to commit to the tables library

import tables
from copy import copy
#from utils import hpd, quantiles
from scipy.stats.mstats import mquantiles
from GetNames import *



def GetParamValues(FileName, ParamName,burn=0,nthin=None):
     '''Get parameter values saved in a hdf5 file.  
     *FileName is the full name of the file.
     *ParamName is the name(s) of the node(s)
     *first burn of the iterations stored in the file are ignored
     *results to be reduced to nthin evenly spaced values'''
     
     #check for multiple nodes
     if isinstance(ParamName,(list,np.ndarray)):
         result=[ GetParamValues(FileName, p,burn=burn,nthin=nthin)  for p in ParamName]
         return(result)
         
     #Check for multiple input files
     if isinstance(FileName,(list,np.ndarray)):
         x=np.ndarray(0,dtype=float)
         for fname in FileName:
           #burn is applied to both chains
           x=np.append(x,GetParamValues(fname, ParamName,burn=burn,nthin=None))
         
         n=len(x)
     
         #No thinning to do     
         if not(nthin) or (nthin>=n):
           return(x)
     
         result=[ x[int(i*(n-1)/(nthin-1))]      for i in range(nthin)]     
         return(result)
     
     #Single node, single file
     #open file and get all values
     f=pymc.database.hdf5.load(FileName)
     try:
       x=f.trace(ParamName)[:]
     except:
       print('GetParamStats 48')
       print(ParamName, ' ', FileName)
       x=f.trace(ParamName)[:]
     n=len(x)
     
     #burn-in not achieved
     if burn>n:
         return([])
     
    #Apply the burn-in     
     x=x[burn:]     
     n=len(x)
     #No thinning to do     
     if not(nthin) or (nthin>=n):
         return(x)
     
     result=[ x[int(i*(n-1)/(nthin-1))]      for i in range(nthin)]     
     return(result)
     
def GetParamValues2(FileName, ParamName,burn=0,nthin=None):
     '''Get parameter values saved in a hdf5 file.  
     Use Table-library directly and avoid pymc code.
     
     *FileName is the full name of the file.
     *ParamName is the name(s) of the node(s)
     *first burn of the iterations stored in the file are ignored
     *results to be reduced to nthin evenly spaced values'''
     
     #check for multiple nodes
     if isinstance(ParamName,(list,np.ndarray)):
         result=[ GetParamValues2(FileName, p,burn=burn,nthin=nthin)  for p in ParamName]
         return(result)
         
     #Check for multiple input files
     if isinstance(FileName,(list,np.ndarray)):
         x=np.ndarray(0,dtype=float)
         for fname in FileName:
           #burn is applied to both chains
           x=np.append(x,GetParamValues2(fname, ParamName,burn=burn,nthin=None))
         
         n=len(x)
     
         #No thinning to do     
         if not(nthin) or (nthin>=n):
           return(x)
     
         result=[ x[int(i*(n-1)/(nthin-1))]      for i in range(nthin)]     
         return(result)
     
     #Single node, single file
     #open file and get all values
     result=OldGetParamValues(FileName, ParamName,burn=burn,nthin=nthin)
     
     return(result)
     
def OldGetParamValues(FileName, ParamName,burn=0,nthin=None):
	'''Get parameter values saved in a hdf5 file.  Chains are appended one after the other.
       Use Table-library directly and avoid pymc code.'''
	
	f=tables.open_file(FileName,mode='r+')

	nTable=len(f.list_nodes('/'))
	
	try:
	  #Names of tables
	  curName='//chain0//PyMCsamples'
	  curTable=f.get_node(curName)
	  values=curTable.col(ParamName)
	except:
	  print('\nGetParamStats 114')
	  print(ParamName,' ',FileName)
	  #Names of tables
	  curName='//chain0//PyMCsamples'
	  curTable=f.get_node(curName)
	  values=curTable.col(ParamName)

	f.close()
	xarray = np.array(values[burn:], float)
	result=thin(xarray,nthin)
	del (f)
	return (result)

    



def thin(x,nthin):
	if nthin==None:return(x)
	n=len(x)
	if nthin>=n:return(x)
	ithin=range(nthin)
	idelta=float(n)/float(nthin)
	iuse=list(map(lambda i:int(i*idelta),ithin))
	result=list(map(lambda i:x[i],iuse))
	return(result)	

def GetParamStats(FileName, ParamName,burn=0):
	"Generate statistics for a parameter saved in a hdf5 file"
	
	xarray=GetParamValues(FileName, ParamName,burn=burn)
	n = len(xarray)
	if not n:
	    print ('Cannot generate statistics for zero-length xarray in', ParamName)
	    return
	 

	try:
	    xquantiles=mquantiles(xarray,prob=[.025,.52,.5,.75,.975])
	    #xhpd=hpd(xarray, .05)
	    return {
		'n': n,
		'standard deviation': xarray.std(0),
		'mean': xarray.mean(0) ,
		#'%s%s HPD interval' % (int(100*(.95)),'%'): xhpd,
		'mc error': xarray.std(0) / sqrt(n)  ,
		'quantiles': xquantiles
	    }
	except:
	    print ('Could not generate output statistics for', ParamName)

def GetSumParamValues(FileName, PreFix,burn=0,nthin=None):
    nPrefix=len(PreFix)
    pname=GetNames(FileName)
    pname=list(filter(lambda p:p[:nPrefix]==PreFix,pname))
    AllValues=list(map(lambda p: GetParamValues(FileName, p,burn=burn,nthin=nthin)   ,pname))
    nparam=len(pname)
    niter=len(AllValues[0])
    result=list(map(lambda i:sum( list(map(lambda t:t[i],AllValues))) ,range(niter)))
    return(result)



if __name__ == "__main__":
    
    import os
    import sys
    os.chdir('Q:/analyses/CrossDating/ModelsAndFunctions')
    sys.path.append('Q:/analyses/CrossDating/ModelsAndFunctions')
    import tables
    
    FileName="Q:/analyses/CrossDating/Barkley/Ring/AutoCorrRecruit/4plus.140/MCMC.hdf5"
    FileName="Q:/analyses/CrossDating/Barkley/Ring/AutoCorrRecruit/4plus/MCMC.hdf5"
    ParamName="Mort"
    x=GetParamStats(FileName, ParamName)
    print (x)
    
    

