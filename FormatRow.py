# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:31:36 2020

@author: HajasW
Format a list of values to a string.  
Useful for making csv files.

"""
from numpy import ndarray

def FormatRow(x,colwidth=20,sigdigit=8):
    if isinstance(x,(list,ndarray)):
        result=''
        for t in x:
            result+=(FormatRow(t,colwidth=colwidth,sigdigit=sigdigit)+',')
        return([result])
        
    if isinstance(x,str):
        trunc=x[:colwidth-1]
        fmt='{:>'+str(colwidth)+'}'
        result=fmt.format(x)
        return(result)
    if isinstance(x,int):
        result=FormatRow(str(x),colwidth=colwidth,sigdigit=sigdigit)
        return(result)    
    if isinstance(x,float):
        fmt='{:'+str(colwidth-1)+'.'+str(sigdigit)+'}'
        xchr=fmt.format(x)
        result=FormatRow(xchr,colwidth=colwidth,sigdigit=sigdigit)
        return(result)    
  
if __name__ == "__main__":
    print(FormatRow(1.5454545,colwidth=15,sigdigit=8))       
    print(FormatRow(1.54545  ,colwidth=15,sigdigit=8))       
    print(FormatRow(15454545,colwidth=15,sigdigit=8))        
    print(FormatRow('15454545',colwidth=15,sigdigit=8))  

    print(FormatRow([15454.545,'15xx4545',15454545,1.5454545,10.5454545,1.54545],colwidth=15,sigdigit=8))       
    print(FormatRow([1.5454545,10.5454545,1.54545,15454.545,'154xxx45',15454545],colwidth=15,sigdigit=8))       
    print(FormatRow([1.54545,15454.545,1.5454545,10.5454545,'1xxxx545',15454545],colwidth=15,sigdigit=8))       

        