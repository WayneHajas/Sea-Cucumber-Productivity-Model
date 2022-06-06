import sys
import pickle
from numpy.random import seed,normal
from random import choices
from numpy import average,exp
from scipy.stats.mstats import mquantiles
from numpy.random import normal,seed


sys.path.append('c:/analyses/EFA_Productivity_Model.20180929/pyfunctions')
from GetParamStats import GetParamValues2 as GetParamValues
from GetCoastLength import GetCoastLength

def GetCoastLength(pickle_file):
    SiteNumber,Sites,VirginSites,CoastLength,timestamp=pickle.load(open(pickle_file,"rb"))
    return(CoastLength[0]+CoastLength[2]+CoastLength[4]+CoastLength[8]+CoastLength[16])

    

class USR2():
    def __init__(self, hdf5file,PickleFile,z=4,nquant=1000,burn=0,nthin=None,seed=None):

        self.sdYear=GetParamValues(hdf5file, 'sdYear',burn=burn,nthin=nthin)
        self.sdSite=GetParamValues(hdf5file, 'sdSite',burn=burn,nthin=nthin)
        n=len(self.sdYear)
        self.CoastLength=GetCoastLength(PickleFile)
        
        vbmass=[GetParamValues(hdf5file, t,burn=burn,nthin=nthin)  for t in ['VBmass_0', 'VBmass_2', 'VBmass_4', 'VBmass_8', 'VBmass_16' ]  ]
        self.vbmass=[ sum([t[i] for t in vbmass] ) for i in range(n)  ]
        
        YearEffect=normal(0,self.sdYear)
        SiteEffect=normal(0,self.sdSite)
        self.USR21=[exp(-z*t)  for t in self.sdYear ]
        self.USR22=[self.USR21[i]*self.vbmass[i]/self.CoastLength*exp(SiteEffect[i]+YearEffect[i])   for i in range(n)]
        
        q=[ (i+.5)/nquant for i in range(nquant)]
        result={}
        self.qUSR21=mquantiles(self.USR21,prob=q)
        self.qUSR22=mquantiles(self.USR22,prob=q)
        
    

if __name__ == "__main__":
    sys.path.append('c:\\analyses\\EFA_Productivity_Model.20180929\\Jervis\\NewModel_wideSdYear_WideSiteArea')
    from hdf5file  import hdf5file,burn,nthin,PickleFile
    hdf5file=hdf5file[0]
    SiteNumber,Sites,VirginSites,CoastLength,timestamp=pickle.load(open(PickleFile,"rb"))
    CL=GetCoastLength(PickleFile)
    
    test=USR2(hdf5file,PickleFile,nquant=1000,burn=burn,nthin=nthin,seed=20200424)
