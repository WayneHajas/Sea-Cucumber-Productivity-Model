# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 09:40:27 2021

@author: HajasW

Tools for applying logit and expit.
Make useful transformamtions for appling bounds with root-functions and minimizations.
"""
from numpy import log,exp

def wchlogistic(x,xmid,steep,ymin,ymax):
    ynorm=1/( 1+exp(-steep*(x-xmid)) )
    
    result=ymin+(ymax-ymin)*ynorm
    return(result)

def wchlogit(y,xmid,steep,ymin,ymax):
    ynorm=(y-ymin)/(ymax-ymin)
    result=xmid-log((1-ynorm)/ynorm)/steep
    return(result)

if __name__ == "__main__":
    ymin,ymax,steep=0,1,1
    xmid=0
    
    x=-3
    ytest=wchlogistic(x,xmid,steep,ymin,ymax)
    xtest=wchlogit(ytest,xmid,steep,ymin,ymax)
    print(x,ytest,xtest)
    print()
    x=3
    ytest=wchlogistic(x,xmid,steep,ymin,ymax)
    xtest=wchlogit(ytest,xmid,steep,ymin,ymax)
    print(x,ytest,xtest)
    print()
    
    
    print()
    n=1000
    xmin,xmax,steep=-10,10,.25
    ymin,ymax,xmid=.35,.9,(xmin+xmax)/2
    x=[xmin+(xmax-xmin)*(t+.5)/n for t in range(n)]
    ytest=[ wchlogistic(t,xmid,steep,ymin,ymax)  for t in x]    
    xtest=[ wchlogit(t,xmid,steep,ymin,ymax)  for t in ytest]
    print(ytest[0],ytest[-1])
    
    #for i in range(n):
        #print(x[i],ytest[i],xtest[i])
    
    import matplotlib.pyplot as plt
    plt.plot(x,ytest,'k-')
    plt.show()