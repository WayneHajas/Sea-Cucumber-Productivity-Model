from scipy.optimize import brentq
from numpy import array, ndarray,exp
from NextRelBmass import NextRelBmass

def ErrAvgBmass_bySigma(Harvest,a,b,xmax,fmax,sigma,nsigma,deltaT):
    TargAvgBmass=exp(-nsigma*sigma)
    B0=TargAvgBmass-Harvest/2
    B1=NextRelBmass(OldBiomass=B0,OldTime=0,NewTime=deltaT,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=1)
    result=(B1+B0)/2-TargAvgBmass
    return(result)

def ErrAvgBmass(Harvest,a,b,xmax,fmax,deltaT,TargAvgBmass):
    B0=TargAvgBmass-Harvest/2
    B1=NextRelBmass(OldBiomass=B0,OldTime=0,NewTime=deltaT,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=1)
    result=(B1+B0)/2-TargAvgBmass
    return(result)

def FindHarvestAmount_bySigma(a,b,xmax,fmax,sigma,nsigma,deltaT,minGuess=0, maxGuess=100):
    if isinstance (a, (list,ndarray)):
        n=len(a)
        result=[ FindHarvestAmount_bySigma(a[i],b[i],xmax[i],fmax[i],sigma[i],nsigma,deltaT,minGuess=minGuess, maxGuess=maxGuess) for i in range(n)]
        return(result)
    
    try:
        result=brentq(ErrAvgBmass_bySigma,minGuess,maxGuess,args=(a,b,xmax,fmax,sigma,nsigma,deltaT))
    except:
        print('Avg50Biomass 27')
        print(a,b,xmax,fmax,deltaT)
        result=brentq(ErrAvgBmass_bySigma,minGuess,maxGuess,args=(a,b,xmax,fmax,sigma,nsigma,deltaT))
    return(result)

def FindHarvestAmount(a,b,xmax,fmax,deltaT,TargAvgBmass=0.5,minGuess=0, maxGuess=100):
    if isinstance (a, (list,ndarray)):
        n=len(a)
        result=[ FindHarvestAmount(a[i],b[i],xmax[i],fmax[i],deltaT,TargAvgBmass=TargAvgBmass,minGuess=minGuess, maxGuess=maxGuess) for i in range(n)]
        return(result)
    
    try:
        result=brentq(ErrAvgBmass,minGuess,maxGuess,args=(a,b,xmax,fmax,deltaT,TargAvgBmass))
    except:
        print('Avg50Biomass 20')
        print(a,b,xmax,fmax,deltaT)
        result=brentq(ErrAvgBmass,minGuess,maxGuess,args=(a,b,xmax,fmax,deltaT,TargAvgBmass))
    return(result)
    
if __name__ == "__main__":  
    Harvest=0.05
    a,b,xmax,fmax,deltaT=0.9326984812362484, 0.12112333934010579, 0.8850627905257634, 0.06238420792442966, 0.019230769230769232
    TargAvgBmass=0.5
    minGuess,maxGuess=0.001,.95
    
    print()
    test1=    ErrAvgBmass(Harvest,a,b,xmax,fmax,deltaT,TargAvgBmass)
    print(test1)
    
    print()
    test2=    FindHarvestAmount(a,b,xmax,fmax,deltaT,TargAvgBmass=0.5,minGuess=0, maxGuess=100)
    print(test2)
    
    print()
    deltaT=3
    test3=    FindHarvestAmount(a,b,xmax,fmax,deltaT,TargAvgBmass=0.5,minGuess=0, maxGuess=100)
    print(test3)
    
    a=[.1,.3,.5,.7,.9]
    b=[1-t  for t in a]
    n=len(a)
    xmax=[a[i]/(a[i]+b[i]) for i in range(n)]
    fmax=[t/4 for t in a]
    print()
    test4=FindHarvestAmount(a,b,xmax,fmax,deltaT,TargAvgBmass=0.5,minGuess=0, maxGuess=100)
    print(test4)
    
     
    Harvest=0.05
    a,b,xmax,fmax,deltaT,sigma=0.9326984812362484, 0.12112333934010579, 0.8850627905257634, 0.06238420792442966, 0.99,0.188
    nsigma=4
    TargAvgBmass=0.5
    minGuess,maxGuess=0.001,.95
    
    print()
    test6=    FindHarvestAmount_bySigma(a,b,xmax,fmax,sigma,nsigma,deltaT,minGuess=minGuess, maxGuess=maxGuess)
    print('test6',test6)
    
    a=[.1,.3,.5,.7,.9]    
    b=[1-t  for t in a]
    n=len(a)
    xmax=[a[i]/(a[i]+b[i]) for i in range(n)]
    fmax=[t/4 for t in a]
    sigma=[t for t in a]
    print()
    test7=    FindHarvestAmount_bySigma(a,b,xmax,fmax,sigma,nsigma,deltaT,minGuess=minGuess, maxGuess=maxGuess)
    print('test7',test7)
    