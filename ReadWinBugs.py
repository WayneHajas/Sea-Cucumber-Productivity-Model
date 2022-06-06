import csv
from numpy import ndarray
from GetParamStats import thin

def ReadWinBugs(IndFile,OutFile,ParamName,burn=0,nthin=None):
    if isinstance(ParamName,(list,ndarray)):
        result={}
        for p in ParamName:
            result[p]=ReadWinBugs(IndFile,OutFile,p,burn=burn,nthin=nthin)
        return(result)
  
    if isinstance(IndFile,(list,ndarray)):     
        x=[]
        for i,ifile in enumerate(IndFile):
            x+=ReadWinBugs(ifile,OutFile[i],ParamName,burn=burn,nthin=None)
        if nthin:
            x=thin(x,nthin)
        return(x)
        
    #Single parameter and single ind-out pair
    with open(IndFile,newline='\n') as f:
        reader=csv.reader(f,delimiter='\t')
        for row in reader:
           if row[0]==ParamName:
               index1,index2=int(row[1]),int(row[2])
               break
    
    with open(OutFile,newline='\n') as f:           
        reader=csv.reader(f,delimiter='\t')
        i=1
        for row in reader:
            if i==index1:
                result=[float(row[1])]
                break
            i+=1
        for row in reader:
            result+=[float(row[1])]
            i+=1
            if i>=index2:
               break
           
    result=result[burn:]
    n=len(result)
    if nthin and (n>nthin):
        result=thin(result,nthin)
    return(result)

def GetParamName(IndFile):

    with open(IndFile,newline='\n') as f:
        reader=csv.reader(f,delimiter='\t')
        result=[]
        for row in reader:
            result+=[row[0]]
    return(result)
    
if __name__ == "__main__":
      burn=100
      nthin=1500
    
      IndFile=['I:\\Archive\\s-drive\\analyses\\CukeExpHarvRecov\\Jervis.beta\\chain1.ind' ,'I:\\Archive\\s-drive\\analyses\\CukeExpHarvRecov\\Jervis.beta\\chain2.ind']       
      OutFile=['I:\\Archive\\s-drive\\analyses\\CukeExpHarvRecov\\Jervis.beta\\chain1.out' ,'I:\\Archive\\s-drive\\analyses\\CukeExpHarvRecov\\Jervis.beta\\chain2.out']     
      ParamName='deviance'
      
      test1=ReadWinBugs(IndFile,OutFile,ParamName,burn=burn,nthin=nthin)
      print(len(test1))
      
      test2=GetParamName(IndFile[0])
      for t in test2:
          print(t)