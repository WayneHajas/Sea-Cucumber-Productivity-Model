'''
2016-01-25
Recreated this library to more directly use existing pymc-libraries
'''

'''Get the names of the traceable nodes associated with a pymc-MCMC.
Ignore names beginning with 'Metropolis_'

2021-02-16
    Use tables library instead of pymc library.  Circumvent some installation issues
    
2021-06-02
    stop as soon as there are names
'''

import  numpy as np
import tables

def GetNames(FileName):
    
    #Check for multiple input files
    if isinstance(FileName,(list,np.ndarray)):
        result=GetNames(FileName[0])
        return(result)
    
    f=tables.open_file(FileName,mode='r+')
    
    #Names of tables
    curName='//chain0//PyMCsamples'
    curTable=f.get_node(curName)
    curnames=None
    for t in curTable:
        try:
          curnames=t.table.row.table.colnames
        except:
          dummy=True
        if curnames:
            break
    if not curnames:
        return(None)
    
    curnames=[t for t in curnames if t[:11]!= 'Metropolis_']
    return(curnames)

