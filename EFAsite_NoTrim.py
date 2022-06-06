'''
20200623
    Specical version of EFAsite with with an extra event one day after harvest.  
    That shold result in post-harvest values of relative biomass
    
20210423
    Special version that does not use trimmed transect
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
  
  def ReadTransect(self,track=True,DateByMonth=False):
    global ConstantValues
    query='SELECT DISTINCT Headers.Transect '
    query+='FROM Densities INNER JOIN Headers ON Densities.HKey = Headers.Key '
    query+='WHERE (((Headers.Project)="'
    query+=self.Survey
    query+='") AND ((Headers.EFA) Is Not Null) '
    query+=  ' AND (Headers.Year>= '
    query+=         str(ConstantValues.MinYear)+')'
    query+=  ' AND (Headers.Year<= '
    query+=         str(ConstantValues.MaxYear)+')'
    query+=') AND ((Headers.Site)= '
    query+=str(self.SiteNumber) 
    query+=') '
    query+='order by Headers.Transect ;'

    dataSource=daoBaseClass(self.mdbfile,query)
    
    #Transect numbers
    self.TranNum=dataSource.GetVariable('Transect')
    del dataSource
    
    #Create transects.  Each transect will be truncated to have the same lenght and similar depth-profile every time it is surveyed.    
    self.SurvTran=list(map(lambda t: Transect(self.mdbfile,self.Survey,self.SiteNumber,t,WithTrim=False,track=track,DateByMonth=DateByMonth),  \
                      self.TranNum))
    self.SurvTran=list(filter(lambda t:len(t.AllQuad)>0,self.SurvTran))
    
    #Dates when surveys occur.
    self.SetDayOfSurvey()
    
    #Create virgin-transect corresponding to the first Survey-date.  No truncation
    self.VirgTran=[ VirginTransect(t,vyear=self.VirginYear,DateByMonth=DateByMonth)  for t in self.SurvTran]    
    self.VirgTran=[t  for t in self.VirgTran if len(t.ydata)>0]

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

