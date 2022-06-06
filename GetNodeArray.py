''''
2019-07-15 Class to aid in reading and managing values of arrays of nodes'''


from numpy import exp
import os
import sys
sys.path.append('S:\\analyses\\EFA_Productivity_Model.20180929\\pyfunctions')
from GetParamStats import GetParamValues2 as GetParamValues
from GetNames import *
from ToDecimalYear import ToDecimalYear

class ArrayNodes():
    def __init__(self,hdf5file,prefix,burn=0,nthin=None):
        '''
            hdf5file is the name of the output file(s)
            prefix gives the first characters of the node names
            burn is the number of intial iterations in each MCMC to aviod
            nthin is the final number of iterations to consider'''
        self.hdf5file,self.prefix,self.burn,self.nthin=hdf5file,prefix,burn,nthin
        self.GetNodeValues()
    
    def GetNodeNames(self):
      ''' All the node-names that begin with prefix'''
      AllNames=GetNames(self.hdf5file)
      nprefix=len(self.prefix)
      pname=[t for t in AllNames if t[:nprefix]==self.prefix]
      pname.sort()
      return(pname)
    
    def GetNodeValues(self):
        '''Get values corresponding to the prefix
        The values go into a dictionary.  
        Node names are based on node-names
        '''
        self.NodeValues={}
        self.ParamName=self.GetNodeNames()
        for p in self.ParamName:
            self.NodeValues[p]=GetParamValues(self.hdf5file,p,burn=self.burn,nthin=self.nthin)
        self.nthin=len(self.NodeValues[self.ParamName[0]])
    
    def GetRB(self,RBname):
        '''Get values for a particular node-name'''
        try:
            result=self.NodeValues[RBname]
        except:
            print()
            print('GetNodeArray 32')
            print('failed to find data for ' ,RBname )
            result=None
        return(result)

def CalcDecimalYear(NodeName):
    '''Assuming a string is of the form, perfix_site_year_month_day, get the sitenumber and decimal-year'''
    chsymd=NodeName.split('_')[-4:]
    symd=[int(t) for t in chsymd]
    DecimalYear=ToDecimalYear(symd[-3:])
    SiteNumber=int(symd[-4])
    result={}
    result['DecimalYear']=DecimalYear
    result['SiteNumber']=SiteNumber
    return(result)

class MeanWeight(ArrayNodes):
    '''A child-class of ArrayNodes.  
    Specialized for MeanWeight'''
    def __init__(self,hdf5file,burn=0,nthin=None):
        
        self.hdf5file,self.burn,self.nthin=hdf5file,burn,nthin
        self.prefix='MeanWeight'
        self.GetNodeValues()
        self.GenDecimalWeightDate()
        
    def GenDecimalWeightDate(self):
        '''Sites and decimals-years corresponding to mean-weight'''
        self.DecimalWeightDate=[CalcDecimalYear(t)  for t in self.ParamName]
        
    def GetDecimalWeightDate(self,NodeName):
        '''Get the mean-weight posterior distribution corresponding to a node-name.  Same site and closest possible date.
        NodeName is of the form: perfix_site_year_month_day'''
        SurveyDecimalYear=CalcDecimalYear(NodeName)
        SSQ=[ (SurveyDecimalYear['SiteNumber']==t['SiteNumber'])*exp((SurveyDecimalYear['DecimalYear']-t['DecimalYear'])**2)  for t in self.DecimalWeightDate]
        i=SSQ.index(max(SSQ))
        MeanWeightValues=self.NodeValues[self.ParamName[i]]
        return(MeanWeightValues)
        


class YearEffect(ArrayNodes):
    '''A child-class of ArrayNodes.  
    Specialized for year-effects'''
    def __init__(self,hdf5file,burn=0,nthin=None):
        self.hdf5file,self.burn,self.nthin=hdf5file,burn,nthin
        self.prefix='year_'
        self.GetNodeValues()
        
        
    def GetYearEffect(self,Rel_Abund):
        chsymd=Rel_Abund.split('_')[-4:]
        YEname='year_'+chsymd[1]
        try:
            result=self.NodeValues[YEname]
            return(result)
        except:
            print('GetNodeArray 81')
            print('no year effect for ',Rel_Abund )
            return(None)        
        
        


class SiteEffect(ArrayNodes):
    '''A child-class of ArrayNodes.  
    Specialized for site-effects'''
    def __init__(self,hdf5file,burn=0,nthin=None):
        self.hdf5file,self.burn,self.nthin=hdf5file,burn,nthin
        self.prefix='site_'
        self.GetNodeValues()
        
        
    def GetSiteEffect(self,Rel_Abund):
        chsymd=Rel_Abund.split('_')[-4:]
        SEname='site_'+chsymd[0]
        try:
          result=self.NodeValues[SEname]
          return(result)
        except:
            print('GetNodeArray 101')
            print('No site-effect for ',Rel_Abund)
            return(None)

        
if __name__ == "__main__":
    hdf5file=[              'S:\\analyses\\EFA_Productivity_Model.20180929\\Tolmie\\NewModel\\seed.20180824.hdf5',\
			'S:\\analyses\\EFA_Productivity_Model.20180929\\Tolmie\\NewModel\\seed.20180825.hdf5',\
			'S:\\analyses\\EFA_Productivity_Model.20180929\\Tolmie\\NewModel\\seed.20180826.hdf5',\
			'S:\\analyses\\EFA_Productivity_Model.20180929\\Tolmie\\NewModel\\seed.20180827.hdf5']
    PickleFile='S:\\analyses\\EFA_Productivity_Model.20180929\\Tolmie\\NewModel\\Tolmie.pickle'
    burn=0
    nthin=None
    prefix='Rel_Abund_16'  
    test2=ArrayNodes(hdf5file,'lnGrandMean',burn=0,nthin=5)
    #test1=ArrayNodes(hdf5file,prefix,burn=0,nthin=5) 
    test3=MeanWeight(hdf5file[0],burn=0,nthin=5)
    test4=test3.GetDecimalWeightDate('Rel_Abund_16_2011_9_23')
    print()
    test5=YearEffect(hdf5file[0],burn=0,nthin=5)
    test6=test5.GetYearEffect('Rel_Abund_16_2011_9_23')
    print()
    test7=SiteEffect(hdf5file[0],burn=0,nthin=5)
    test8=test7.GetSiteEffect('Rel_Abund_16_2011_9_23')
    print()