'''
2016-01-25
Recreated this library to more directly use existing pymc-libraries
'''

'''Get the names of the traceable nodes associated with a pymc-MCMC.
Ignore names beginning with 'Metropolis_'
'''

import  numpy as np
import pymc

def GetNames(FileName):
    
    #Check for multiple input files
    if isinstance(FileName,(list,np.ndarray)):
        result=GetNames(FileName[0])
        return(result)
    
    #open file
    db=pymc.database.hdf5.load(FileName)
    
    #Get the names
    AllName=[t for t in  db.trace_names[-1]   if t[:11]!= 'Metropolis_']   
    
    #Sort the names - ignoring case
    AllName.sort(key=lambda s: s.lower())
    del (db)
    return(AllName)

if __name__ == "__main__":

     FileName= 'D:\Analyses\CukeNonParamProd\JervisGS5\MCMC.hdf5'
     test=GetNames(FileName)
     dir(test)
	
