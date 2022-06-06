'''
20200206
Python libraries for subratracting one statistical distribution from another - assuming that both are represented by a random sample of values
'''
import csv
from scipy.stats.mstats import mquantiles
from numpy import round
from operator import __add__, __sub__,__mul__,__truediv__

def ReadColumn(csvfile,column=0,header=True):
    '''
    read the column-the column of values as floating-point values
    '''

    with open(csvfile,newline='\n') as openfile:
        csvreader=csv.reader(openfile)
        #If there is a header-row, just read it and ignore it.
        if header:
            ColumnNames=next(csvreader)
            
        #read teh column-th value and force it to a floating-point 
        result=[float(t[column]) for t in csvreader]
    return(result)
        
    

def Diff_PDF_byFile(csv1,csv2,DeterministicTerm=0,nthin=None,q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99],column=0,header=True):
  '''  
  #csv1 and csv2 are csv-files containing two samples of values representing two distributions.  The first column will be considered.
  #Typically, csv1 will contain a sample representing an estimate of abundance.  csv2 will represent a targent reference point.
  #DeterministicTerm represent a deterministic(single) value that will added to the results.
  #Typically, DeterministicTerm will be negative and represent harvest since the measurement of abundance
  #nthin represents the size of arrays that will be considered in the calculations.
  #nthin=NA indicates that the full arrays from the CSV-files will be used
  #When nthin is set to a value, an array of nthin*nthin values will occur in the calculations.
  #nthin=1000 is likely manageable by a reasonably modern laptop or desktop. It will also give reasonable precision.
  #A smaller value nthin will reduce the computer-requirements but also degrade 
  
  # The resulting distribution will be pdf1-pdf2+DeterministicTerm
  # q gives the quantiles of the resulting distribution that will be reported 
  '''
  
  #Get first columns in files.  Assume there is a header.
  pdf1=ReadColumn(csv1,column=column,header=header)
  pdf2=ReadColumn(csv2,column=column,header=header)
  
  result=Diff_PDF(pdf1,pdf2,DeterministicTerm=DeterministicTerm,nthin=nthin,q=q)
  return(result)

  	
def Diff_PDF(pdf1,pdf2,DeterministicTerm=0,nthin=None,q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99]):
  '''
  #pdf1 and pdf2 are samples of values representing two distributions.
  #Typically, pdf1 will represent an estimate of abundance.  csv2 will represent a  reference point.
  #DeterministicTerm represent a deterministic(single) value that will added to the results.
  #Typically, DeterministicTerm will be negative and represent harvest since the measurement of abundance
  #nthin represents the size of arrays that will be considered in the calculations.
  #nthin=NA indicates that the full arrays from the CSV-files will be used
  #When nthin is set to a value, an array of nthin*nthin values will occur in the calculations.
  #nthin=1000 is likely manageable by a reasonably modern laptop or desktop. It will also give reasonable precision.
  #A smaller value nthin will reduce the computer-requirements but also degrade 
  
  # The resulting distribution will be pdf1-pdf2+DeterministicTerm
  # q gives the quantiles of the resulting distribution that will be reported 
  '''
  
  #Thin the values if required
  if nthin:
    pdf1=thin(pdf1,nthin)
    pdf2=thin(pdf2,nthin)
  
  diffVal=[]
  for x1 in pdf1:
    diffVal+=[x1-x2+DeterministicTerm for x2 in pdf2]
 
  result=mquantiles(diffVal,prob=q)
  return(result)


def thin(x,nthin):
  '''
  reduce the size of a random values to size-n
  '''
  
  #Do nothing if nthin is undefined
  if not(nthin):
      return(x)
  
  #Do nothing if the array is alread smaller than nthin
  n=len(x)
  if (n<=nthin):
      return(x)
  
  #sort the original array
  x.sort()
  
  #index of values to use
  index=[int((i+.5)/nthin*n) for i in range(nthin)]  
  result=[x[i] for i in index]
  return(result)
##################
def Arith_PDF(pdf1,pdf2,arithfunction,nthin=None,q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99]):
  '''
  #pdf1 and pdf2 are samples of values representing two distributions.
  #Typically, pdf1 will represent an estimate of abundance.  pdf2 will represent a  reference point.
  
  #arithfunction will typically be one of; __add__,__sub__,__mul__,__truediv__
  	#Must be able to do something like x.arithfunction(y) to subtract y from x
  
 
  #nthin represents the size of arrays that will be considered in the calculations.
  #nthin=NA indicates that the full arrays from the CSV-files will be used
  #When nthin is set to a value, an array of nthin*nthin values will occur in the calculations.
  #nthin=1000 is likely manageable by a reasonably modern laptop or desktop. It will also give reasonable precision.
  #A smaller value nthin will reduce the computer-requirements but also degrade 
  
  # The resulting distribution will be pdf1-pdf2+DeterministicTerm
  # q gives the quantiles of the resulting distribution that will be reported 
  '''
  
  #Thin the values if required
  if nthin:
    pdf1=thin(pdf1,nthin)
    pdf2=thin(pdf2,nthin)
  
  OperatorResult=[]
  for x1 in pdf1:
    OperatorResult+=[arithfunction(x1,x2) for x2 in pdf2]
  if not(q):
  	return(OperatorResult)
  result=mquantiles(OperatorResult,prob=q)
  return(result)
  
#Apply Arith_PDF for addition, subtraction, multiplication and division
def pdf__add__(pdf1,pdf2,nthin=None,q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99]):
  	result=Arith_PDF(pdf1,pdf2,__add__,nthin=nthin,q=q)
  	return(result)  
def pdf__sub__(pdf1,pdf2,nthin=None,q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99]):
  	result=Arith_PDF(pdf1,pdf2,__sub__,nthin=nthin,q=q)
  	return(result)
  
def pdf__mul__(pdf1,pdf2,nthin=None,q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99]):
  	result=Arith_PDF(pdf1,pdf2,__mul__,nthin=nthin,q=q)
  	return(result)
  
def pdf__truediv__(pdf1,pdf2,nthin=None,q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99]):
  	result=Arith_PDF(pdf1,pdf2,__truediv__,nthin=nthin,q=q)
  	return(result)
  	
############
def Arith_PDF_byFile(csv1,csv2,arithfunction,nthin=None,q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99],column=0,header=True):
  '''  
  #csv1 and csv2 are csv-files containing two samples of values representing two distributions.  The first column will be considered.
  #Typically, csv1 will contain a sample representing an estimate of abundance.  csv2 will represent a targent reference point.
  
  #arithfunction will typically be one of; __add__,__sub__,__mul__,__truediv__
  	#Must be able to do something like x.arithfunction(y) to subtract y from x
  	
  #nthin represents the size of arrays that will be considered in the calculations.
  #nthin=NA indicates that the full arrays from the CSV-files will be used
  #When nthin is set to a value, an array of nthin*nthin values will occur in the calculations.
  #nthin=1000 is likely manageable by a reasonably modern laptop or desktop. It will also give reasonable precision.
  #A smaller value nthin will reduce the computer-requirements but also degrade 
  
  # The resulting distribution will be pdf1-pdf2
  # q gives the quantiles of the resulting distribution that will be reported 
  '''
  
  #Get  columns in files.  Assume there is a header.
  pdf1=ReadColumn(csv1,column=column,header=header)
  pdf2=ReadColumn(csv2,column=column,header=header)
  
  result=Arith_PDF(pdf1,pdf2,arithfunction,nthin=nthin,q=q)
  return(result)  	

#Apply Arith_PDF for addition, subtraction, multiplication and division
def pdf_ByFile__add__(csv1,csv2,nthin=None,q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99],column=0,header=True):
  	result=Arith_PDF_byFile(csv1,csv2,__add__,nthin=nthin,q=q,column=column,header=header)
  	return(result)  
def pdf_ByFile__sub__(csv1,csv2,nthin=None,q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99],column=0,header=True):
  	result=Arith_PDF_byFile(csv1,csv2,__sub__,nthin=nthin,q=q,column=column,header=header)
  	return(result)    
def pdf_ByFile__mul__(csv1,csv2,nthin=None,q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99],column=0,header=True):
  	result=Arith_PDF_byFile(csv1,csv2,__mul__,nthin=nthin,q=q,column=column,header=header)
  	return(result)     
def pdf_ByFile__truediv__(csv1,csv2,nthin=None,q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99],column=0,header=True):
  	result=Arith_PDF_byFile(csv1,csv2,__truediv__,nthin=nthin,q=q,column=column,header=header)
  	return(result)  

if __name__ == "__main__":  
    ### Below are examples how the above functions can be used. ###
    csv1='S:\\analyses\\EFA_Productivity_Model.20180929\\rfunctions\\pdf1.csv'
    csv2='S:\\analyses\\EFA_Productivity_Model.20180929\\rfunctions\\pdf2.csv'

    from scipy.stats import norm
    from numpy import array,sqrt
    nsample=10000
    x1=norm.rvs(size=nsample,loc=0,scale=1)
    x2=thin(x1,10)
    q=[.05,.15,.25,.35,.45,.55,.65,.75,.85,.95]

    x3=norm.rvs(size=nsample,loc=10,scale=1)
    x4=norm.rvs(size=nsample,loc=5,scale=2)
    q=[.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99]
    x5=Diff_PDF(x3,x4,DeterministicTerm=-1.5,nthin=2000,q=q)
    
    print(norm.isf(1-array(q),loc=10-5-1.5,scale=sqrt(1*1+2*2)))
    print(x5)
    x6=Diff_PDF_byFile(csv1,csv2,DeterministicTerm=-1,nthin=1000,q=q)
    print(x6)
    
    x7=pdf_ByFile__add__(csv1,csv2,nthin=1000,q=q)
