'''
Utilities to help work with reference points

20210226

Rewrite to make better use of built-in scipy functions

'''


import csv
from GetMinRelAbund import GetMinRelAbund 
from GetParamStats import GetParamValues2 as GetParamValues

def WritecsvMinRelAbund(MinRelAbund,csvMinRelAbund):
  outfile=csv.writer(open(csvMinRelAbund,'w'),lineterminator='\n')
  for t in MinRelAbund:
      outfile.writerow([t])
  del outfile 
      
def ReadcsvMinRelAbund(csvMinRelAbund):
  result=[]
  with open(csvMinRelAbund) as csv_file:
    csv_reader=csv.reader(csv_file,delimiter=',')
    for row in csv_reader:
        result+= [float(row[0])]
  return(result)

def WritecsvqMinRelAbund(q,qMinRelAbund,csvMinRelAbund):
  outfile=csv.writer(open(csvMinRelAbund,'w'),lineterminator='\n')
  for i,t in enumerate( qMinRelAbund):
      outfile.writerow([q[i], t])
  del outfile 
  
  
def ReadcsvqMinRelAbund(csvqMinRelAbund):
  q,MinRelAbund=[],[]
  with open(csvqMinRelAbund) as csv_file:
    csv_reader=csv.reader(csv_file,delimiter=',')
    for row in csv_reader:
        q+=[float(row[0])]
        MinRelAbund+=[float(row[1])]
  return(q,MinRelAbund)

def WritecsvUSRmaxHarvest(USRmaxHarvest,csvUSRmaxHarvest):
    HarvestIntervals=list(USRmaxHarvest.keys())
    HarvestIntervals.sort()
    n=len(USRmaxHarvest[HarvestIntervals[0]])
    outfile=csv.writer(open(csvUSRmaxHarvest,'w'),lineterminator='\n')
    outfile.writerow(HarvestIntervals)
    for i in range(n):
      currow=[USRmaxHarvest[t][i]  for t in HarvestIntervals]
      outfile.writerow(currow)
    del outfile 
    

def ReadcsvUSRmaxHarvest(csvUSRmaxHarvest):
  result={}
  result[1],result[2],result[3],result[4],result[5]=[],[],[],[],[]
  with open(csvUSRmaxHarvest) as csvfile:
    csv_reader=csv.DictReader(csvfile)
    for row in csv_reader:
        result[1]+=[float(row['1'])]
        result[2]+=[float(row['2'])]
        result[3]+=[float(row['3'])]
        result[4]+=[float(row['4'])]
        result[5]+=[float(row['5'])]
  return(result)

def WriteqUSRmaxHarvest(qUSRmaxHarvest,csvqUSRmaxHarvest):
  fieldnames=list(qUSRmaxHarvest.keys())
  n=len(qUSRmaxHarvest[fieldnames[0]])
  outfile=csv.writer(open(csvqUSRmaxHarvest,'w'),lineterminator='\n')
  outfile.writerow(fieldnames)
  for i in range(n):
      currow=[qUSRmaxHarvest[t][i]  for t in qUSRmaxHarvest]
      outfile.writerow(currow)
  del outfile 
  
def inRP(csvqMinRelAbund,q=[.5,.75,.90,.95,.99]):
    if isinstance(csvqMinRelAbund,list):
        result=[]
        for t in csvqMinRelAbund:
            result+=inRP(t,q=q)
        return(result)
    
    with open(csvqMinRelAbund) as csv_file:
      csv_reader=csv.reader(csv_file,delimiter=',')
      result=[]
      for row in csv_reader:
        try: 
            if float(row[0]) in q:
                result+=[[ float(t)  for t in row]+[csvqMinRelAbund]]
        except:
            dummy=True
      del csv_file,csv_reader
      return(result)

def parseRP(infile,outfile,q=[.5,.75,.90,.95,.99]):
    data=inRP(infile,q=q)
    outfile=csv.writer(open(outfile,'w'),lineterminator='\n')
    for t in data:
      outfile.writerow(t)
    del outfile 
  

if __name__ == "__main__":    
    
  q=[.5,.75,.90,.95,.99]
  csvqMinRelAbund='C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqMinRelAbund.csv'
  test1=inRP(csvqMinRelAbund,q=q)
  
  csvqMinRelAbund=[\
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqMinRelAbund.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRenv.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRmaxHarvest.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqComboRP.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqPreBmass.csv',
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqMinRelAbund.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRenv.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRmaxHarvest.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqComboRP.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqPreBmass.csv',
                   
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqMinRelAbund.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRenv.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRmaxHarvest.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqComboRP.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqPreBmass.csv',
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqMinRelAbund.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRenv.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRmaxHarvest.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqComboRP.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqPreBmass.csv',
                   
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqMinRelAbund.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRenv.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRmaxHarvest.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqComboRP.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqPreBmass.csv',
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqMinRelAbund.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRenv.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRmaxHarvest.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqComboRP.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqPreBmass.csv',
                   
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqMinRelAbund.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRenv.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRmaxHarvest.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqComboRP.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqPreBmass.csv',
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqMinRelAbund.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRenv.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRmaxHarvest.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqComboRP.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqPreBmass.csv',
                   \
                   \
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqMinRelAbundLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRenvLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRmaxHarvestLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqComboRPLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqPreBmassLinearBiomass.csv',
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqMinRelAbundLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRenvLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRmaxHarvestLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqComboRPLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Jervis/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqPreBmassLinearBiomass.csv',
                   
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqMinRelAbundLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRenvLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRmaxHarvestLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqComboRPLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqPreBmassLinearBiomass.csv',
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqMinRelAbundLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRenvLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRmaxHarvestLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqComboRPLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Laredo/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqPreBmassLinearBiomass.csv',
                   
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqMinRelAbundLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRenvLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRmaxHarvestLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqComboRPLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqPreBmassLinearBiomass.csv',
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqMinRelAbundLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRenvLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRmaxHarvestLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqComboRPLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Tolmie/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqPreBmassLinearBiomass.csv',
                   
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqMinRelAbundLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRenvLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqUSRmaxHarvestLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqComboRPLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea.2007/ReferencePoints/csvqPreBmassLinearBiomass.csv',
                   
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqMinRelAbundLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRenvLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqUSRmaxHarvestLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqComboRPLinearBiomass.csv',
                   'C:/Analyses/EFA_Productivity_Model.20180929/Zeballos/NewModel_wideSdYear_WideSiteArea/ReferencePoints/csvqPreBmassLinearBiomass.csv'
                   \
                   \
                   
                   
                   ]
  outfile='c:/scratch/parseRP.csv'
  test2=inRP(csvqMinRelAbund,q=q)
  parseRP(csvqMinRelAbund,outfile,q=q)
  
  
  