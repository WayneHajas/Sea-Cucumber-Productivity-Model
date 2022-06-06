'''2019-03-08
        Virgin-transect chosen according to year.
        Assume there is only one survey during that year.
    2017-04-11
        Fix an error in DayTran.DepthProfInRange
        Get the last(shallowest) quadrat
    2017-08-04
        rewrote DayTran.DepthProfInRange to get the correct number of quadrats after trimming
    2017-08-23
        Create Transect.HaveEnoughCuke(LowBnd=5)
        Checks to see if LowBnd or more cukes were observed during any of the surveys
	2017-09-14
	
		modify DayTran.DepthProfInRange
		Resulting profile is from the first quadrat in the depth-range to the last quadrat depth-range.  
		Quadrats in the middle might not be in range.
    2021-06-24
        special version with no 2009 surveys
		
'''

from operator import  attrgetter
from numpy import ndarray
import sys
from sys import maxsize as maxint
from ADO import adoBaseClass as daoBaseClass
from ToDecimalYear import ToDecimalYear,FromDecimalYear
from constants import ConstantValues

from transect import DayTran,VirginTransect
from transect import Transect as oldTransect

  
    		
class Transect(oldTransect):    		
  def __init__(self,mdbfile,Survey,Site,TransectNumber,WithTrim=False,track=True,DateByMonth=False):
    if track:print(' TransectNumber ',Survey,Site,TransectNumber)
    self.mdbfile=mdbfile
    self.Survey=Survey
    self.Site=Site
    self.TransectNumber=TransectNumber
    self.GetDayOfSurvey()
    try:
        self.AllQuad=[DayTran(self.mdbfile,self.Survey,self.Site,self.TransectNumber,DOS,DateByMonth=DateByMonth) for DOS in self.DayOfSurvey]
   
    except:
        self.AllQuad=[DayTran(self.mdbfile,self.Survey,self.Site,self.TransectNumber,DOS,DateByMonth=DateByMonth) for DOS in self.DayOfSurvey]
    try:
      self.AllQuad=[aq for aq in self.AllQuad   if len(aq.ydata)>0]
    except:
      harvestObject.set_trace()
      self.AllQuad=[aq for aq in self.AllQuad   if len(aq.ydata)>0]
    
    #Reduce Day-of-Survey to remaining surveys    
    self.DayOfSurvey=[{'year':aq.year,'month':aq.month,'day':aq.day}  for aq in self.AllQuad]
    self.DayOfSurvey=sorted( self.DayOfSurvey,key=lambda k:(k['year'],k['month'],k['day']))
    self.DecimalYear=ToDecimalYear(self.DayOfSurvey)

    self.GetCommonDepthRange()
    try:
      if WithTrim:self.Trim()
    except:
      print ('transect 169')
      print ('self.TransectNumber',self.TransectNumber)
      if WithTrim:self.Trim()
      

  def GetDayOfSurvey(self):
    global ConstantValues
    query= 'SELECT Headers.Year, Headers.Month, Headers.Day '
    query+='FROM Headers INNER JOIN Densities ON Headers.Key = Densities.HKey '
    query+='GROUP BY Headers.Project, Headers.Site, Headers.Transect, Headers.Year, Headers.Month, Headers.Day '
    query+='HAVING ( '
    query+=     '(Headers.Project= "'
    query+=         self.Survey
    query+=         '") AND '
    query+=     '(Headers.Site= '
    query+=         str(self.Site)+' ) '
    query+=     ' AND (Headers.Transect= '
    query+=         str(self.TransectNumber)+')'
    query+=     ' AND (Headers.Year>= '
    query+=         str(ConstantValues.MinYear)+')'
    query+=     ' AND (Headers.Year<= '
    query+=         str(ConstantValues.MaxYear)+')'
    query+=     ' AND (Headers.Year<> 2009) '
    query+= ') '
    query+=' order by  Headers.Year,Headers.Month,Headers.Day   '    
    query+=';'
    
    try:
        dataSource=daoBaseClass(self.mdbfile,query)
    except:
        print ('149 ',query)
    self.DayOfSurvey=[ {'year':t[0],'month':t[1],'day':t[2]}  for t in dataSource.GetALL()]
    del dataSource




if __name__ == "__main__":
  mdbfile='t:\SeaCuke_Bio.mdb'
  Survey='Laredo Inlet'
  Site=8
  TransectNumber=5
  YearMonth=[1999,2,18]
  WithTrim=True
  
  global ConstantValues
  ConstantValues.MinYear=1998

  
  
  t=Transect(mdbfile,Survey,Site,TransectNumber,WithTrim=WithTrim)
  print('transect 284')
  Profile=t.DepthProfInRange()
  print('transect 286')
  v=VirginTransect(t)
  print(v.Summary())

  print()
  print(t.GetBestProfile(3))
  print()
  for s in t.AllQuad:
      d=[r[1] for r in s.ydata]
      print(d)
  print()  
  for s in t.AllQuad:
      print(s.Summary()['ncuke'])
  
  
  print ('done')
  
  

