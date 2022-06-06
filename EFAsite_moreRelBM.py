'''
20200623
    Specical version of EFAsite with with an extra event one day after harvest.  
    That shold result in post-harvest values of relative biomass
    '''

import sys
from numpy import ndarray
import numpy
from ADO import adoBaseClass as daoBaseClass


from transect import Transect, VirginTransect
from Harvest import Harvest
from constants import ConstantValues
from UniqueDates import UniqueDates
from ToDecimalYear import ToDecimalYear,FromDecimalYear,CombineDateLists
from MeanWeight import MeanWeight

from EFAsite import Site as oldSite
from EFAsite import VirginSite

class Site(oldSite):
  def __init__(self,mdbfile,Survey,SiteNumber,CoastLength,DateByMonth=False):
    '''
    Site(mdbfile,Survey,SiteNumber,CoastLength,DateByMonth=False)
    
    Class to represent an EFA-site
    If DateByMonth, survey-dates will be shifted to the nearest first-day of the month.  That will reduce the number of Survey-dates.
    
    '''
    print('\n ')
    self.mdbfile=mdbfile
    self.Survey=Survey
    self.SiteNumber=SiteNumber
    self.CoastLength=CoastLength
    
    self.SiteHarv=Harvest(self.mdbfile,self.Survey,self.SiteNumber,MinYear=ConstantValues.MinYear,MaxYear=ConstantValues.MaxYear)
    self.ReadTransect(DateByMonth=DateByMonth)
    self.MeanWeight=MeanWeight(mdbfile,Survey,Site=self.SiteNumber)
    self.DateByMonth=DateByMonth
  
        
  def GetEventTimes(self, MaxTimeInc=sys.maxsize):
      '''
      Site.GetEventTimes(MaxTimeInc=sys.maxsize)
      Get time values where survey and/or harvest occurs.

      If there is an increment greater than  MaxTimeInc between events, more time-values
      will be included so the maximum time-increment is never greater than MaxTimeInc    
      '''
      print('EFAsite_moreRelBM 52')
      #Time Values that occur in data
      MaxTime=max(self.DayOfSurvey)#Date of last survey
      raw=[t for t in self.SiteHarv.HarvestDates if t<=MaxTime ]#Harvests before last survey
      raw+=[ToDecimalYear(FromDecimalYear(t+1.5/365.25)) for t in raw]
      raw+=self.DayOfSurvey#Combine harvest and survey dates.  Reduce to sorted unique values
      raw=list(set(raw))
      raw.sort()
      
      #Use these values if the maximum time-increment is set to its maximum
      if MaxTimeInc>=sys.maxsize:
          return(raw)
      if not(raw):
          return(raw)
      if len(raw) ==1:
        return(raw)
      
      #Build a list of event-times.  Where required, add more time-values to 
      #satisfy MaxTimeInc
      oldt=raw[0]
      result=[raw[0]]
      for newt in raw[1:]:
          #Do not need to break down the time increment
          if (newt-oldt)<=MaxTimeInc:
              result+=[newt]
          else:
              nincr=int(numpy.ceil((newt-oldt)/MaxTimeInc))
              delta=(newt-oldt)/nincr
              result+=[ oldt+delta*i for i in range(1, nincr+1)]
          oldt=newt
      return(result)  
  

if __name__ == "__main__":
  mdbfile='D:\Analyses\CukeNonParamProd\SeaCuke_Bio.mdb'
  Survey='Jervis Inlet'
  SiteNum=8

  global ConstantValues
  ConstantValues.MinYear=1995
  
  CoastLength=1000
  s=Site(mdbfile,Survey,SiteNum,CoastLength)
  v=VirginSite(s)

  y=s.GetSurveyYears()
 
  for t in s.SurvTran:
      year=[ty.Summary()['year']  for ty in t.AllQuad]
      ncuke=[ty.Summary()['ncuke']  for ty in t.AllQuad]
      nquad=[ty.Summary()['nquad']  for ty in t.AllQuad]
      print(t.TransectNumber,year,nquad[0],ncuke)
  print()
  s.RemoveLowPopTransect(LowBnd=5)
  for t in s.SurvTran:
      year=[ty.Summary()['year']  for ty in t.AllQuad]
      ncuke=[ty.Summary()['ncuke']  for ty in t.AllQuad]
      nquad=[ty.Summary()['nquad']  for ty in t.AllQuad]
      print(t.TransectNumber,year,nquad[0],ncuke)
  
  
  print ('done')

