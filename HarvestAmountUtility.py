'''
Utilities to document the impact of constant harvest amounts

20210127

Rewrite to make better use of built-in scipy functions

'''

from scipy.integrate import odeint
from scipy.optimize import root_scalar
from scipy.stats.mstats import mquantiles
from scipy.optimize import minimize_scalar
from numpy import log,exp,average,ndarray
from numpy.random import normal,seed


#from scipy.special import expit,logit
from wchLogistic import wchlogit as logit
from wchLogistic import wchlogistic as expit
from NextRelBmass import NextRelBmass

def Prod(B,t,a,b,xmax,fmax):
    if B<0.0001:
        return(0)
    if B>0.9999:
        return(0)
    try:
        result=fmax*((B/xmax)**a) * (((1-B)/(1-xmax))**b)
    except:
        print()
        print('HarvestAmountUtility 25')
        print(B,t,a,b,xmax,fmax)
        result=fmax*((B/xmax)**a) * (((1-B)/(1-xmax))**b)
    return(result)

def EstPreB(Bpost,t,a,b,xmax,fmax):
    if Bpost>=0.999:
        return(Bpost)
    args=(a,b,xmax,fmax,)
    try:
       PreB= odeint(Prod,Bpost,[0,t],args=args)[-1][0]
    except:
       print('Harvest AmountUtility 41')
       print(Bpost,t,a,b,xmax,fmax)
       print()
       PreB= odeint(Prod,Bpost,[0,t],args=args)[-1][0]
    return(PreB)

def rootEstPreB(Bpost,target,t,a,b,xmax,fmax):
    Bpre= EstPreB(Bpost,t,a,b,xmax,fmax)
    result=Bpre-Bpost-target
    return(result)

def lminEstPreB(lBpost,target,t,a,b,xmax,fmax,xmid,steep,ymin,ymax):
    Bpost=expit(lBpost,xmid,steep,ymin,ymax)
    test=EstPreB(Bpost,t,a,b,xmax,fmax)
    negsquare=((target)-(test))**2
    return(negsquare)
    
def lrootEstPreB(lBpost,target,t,a,b,xmax,fmax,xmid,steep,ymin,ymax):
    '''Use logistic transform
    '''
    Bpost=expit(lBpost,xmid,steep,ymin,ymax)
    try:
        result=rootEstPreB(Bpost,target,t,a,b,xmax,fmax)
    except:
        print('HarvestAmountUtility 59')
        print(Bpost,target,t,a,b,xmax,fmax)
        print()
    return(result)

def FindBpost(target,t,a,b,xmax,fmax):
    args=(target,t,a,b,xmax,fmax,)
    Bpost=root_scalar(rootEstPreB,args=args, x0=max([.001,xmax-.001]), x1=min([0.999,xmax+.001]) ).root
    if isinstance(Bpost,complex):
        return(0)
    return(Bpost)

#Functions to find Bpost corresponding to maximum productivity over harvest interval
def minFunc(Bpost,HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund):
    '''
    Evaluate Bpost according to sustainable harvest amount
    '''

    #Calculate the pre-harvest relative biomass
    Bpre=EstPreB(Bpost,HarvestInterval,a,b,xmax,fmax)

    RelHarv=Bpre-Bpost
    #Take negative because scipy utilities want to find the minimum
    return(-RelHarv)

def MaxRelHarvest(HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund):
    '''
    Find the post-harvest relative biomass corresponding to the maximum sustainable harvest amount
    '''
    
    #If there are lists, do operation for each iteration in markov chain 
    if isinstance(fmax,(list,ndarray)):
        n=len(fmax)
        result=[  MaxRelHarvest(HarvestInterval,SimInterval,fmax[i],a[i],b[i],xmax[i],MinRelAbund[i],MaxRelAbund[i])   for i in range(n)]
        return(result)
    
    #One iteration of markov chain
   
    #Check for trivial results 
    if xmax<(MinRelAbund):
            return(MinRelAbund)
    y0=EstPreB(MinRelAbund     ,HarvestInterval,a,b,xmax,fmax)  -   MinRelAbund
    y1=EstPreB(MinRelAbund+1e-4,HarvestInterval,a,b,xmax,fmax)  -   MinRelAbund-1e-4
    if y0>y1:
        return(MinRelAbund) #Greatest productivity at minimum value of post-harvest relative-biomass
    
        
        
    #Values for bracket
    Blow=max([MinRelAbund,           xmax-HarvestInterval*fmax/2]) 
    Bupp=min([.999,                  xmax])
    
    args=(HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund)
    Bpost=minimize_scalar(minFunc,args=args,method='Brent',bracket=(Blow,Bupp),bounds=[xmax,0.999],tol=0.001).x
    
    
    #Special Cases
    if not(Bpost):
        return(0)
    if Bpost<MinRelAbund:
        return(0)
    if Bpost>0.999:
        return(0.999)
    # Calculation is OK:
    return(Bpost)

def lFindBpost(target,t,a,b,xmax,fmax,MinRelAbund):
    '''operate on logit-scale so biomass values remain within range'''
    
    
    if isinstance(a, list):
        n=len(a)
        result=[ lFindBpost(target,t,a[i],b[i],xmax[i],fmax[i],MinRelAbund[i])  for i in range(n)]
        return(result)
    #Single iteration
    
    #trivial cases 
    if target>(t*fmax):
        return(0)
        
    x0=EstPreB(.9985,-t,a,b,xmax,fmax)#Post-harvest value that will get you to virgin state
    if (0.9985-x0)>target:
       return(0.9985-target)#virgin state achieved
    
    x1=xmax-fmax*t/2#Approximation to post-harvest value with maximum productivity
    x1=max([x1,MinRelAbund])
    maxHarv=EstPreB(x1,t,a,b,xmax,fmax)-x1
    if maxHarv<target:
        return(0)#target harvest cannot be sustained
    xmid,steep=(x0+x1)/2,0.25
    args=(target,t,a,b,xmax,fmax,xmid,steep,x0,x1,)
 
    

    try:
       lBpost=root_scalar(lrootEstPreB,args=args, x0=-.25,x1=.25 ).root
    except:
       print('HarvestAmountUtility 152')
       print(args,x0,x1)
       print()
       #lBpost=root_scalar(lrootEstPreB,args=args, x0=lx0, x1=lx1 ).root
       lBpost=minimize_scalar(lminEstPreB,args=args,method='Brent',tol=0.001).x

    if isinstance(lBpost,complex):
        return(0)
    Bpost=expit(lBpost,xmid,steep,x0,x1)
    if Bpost<MinRelAbund:
        return(0)
    return(Bpost)

def qFindBpost(RelHarv,HarvestInterval,fmax,a,b,xmax,MinRelAbund,q=[.025,.5,.975]):
    '''
    Generate quantiles/CredibleIntervals for a given relative harvest and harvest interval
    '''
    
    #Get estimates of post-harvest relative abundance 
    PostHarvest=lFindBpost(RelHarv,HarvestInterval,a,b,xmax,fmax,MinRelAbund)
    
    #Get the quantiles
    qPostHarvest=mquantiles(PostHarvest,prob=q)
    return(qPostHarvest)
    
def curvesFindBpost(RelHarv,HarvestInterval,fmax,a,b,xmax,MinRelAbund,q=[.025,.5,.975]):
    '''
    Generate quantiles/CredibleIntervals for a range relative harvest and  a  harvest interval
    RelHarv is a sorted list of values between zero and one
    '''
    result={}
    result['RelHarv']=RelHarv
    
    #quantiles by relativeHarvest
    byrh=[ qFindBpost(rh,HarvestInterval,fmax,a,b,xmax,MinRelAbund,q=q)   for rh in RelHarv]
    
    #Pick out the quantiles and put them in the results
    for i,q2 in enumerate(q):
        result[q2]=[t[i] for t in byrh]
        
    return(result)

def probFindBpost(RelHarv,HarvestInterval,fmax,a,b,xmax,RP):
    '''
    RP is a reference point.
    Returns the probability of violating the reference point for a harvest amount expressed as a fraction of virgin biomass
    '''
    
    #Get estimates of post-harvest relative abundance 
    PostHarvest=lFindBpost(RelHarv,HarvestInterval,a,b,xmax,fmax,RP)
    
    #Get the quantiles
    ProbVioate=average([t==0 for t in PostHarvest])
    return (ProbVioate)
    

    
###Harvest Amount as a fration of current biomass
 
def rootCurFunc(Bpost,   FracCur,HarvestInterval,a,b,xmax,fmax,MinRelAbund):
    '''
    FracCur is the harvest amount as a fraction of pre-harvest abundance
    When the result of this function is zero, relative biomass starts out at Bpost after harvest.
    '''

    #Calculate the pre-harvest relative biomass
    Bpre=EstPreB(Bpost,HarvestInterval,a,b,xmax,fmax)
    
    #Harvest as a fracion of virgin biomass
    RelHarv=Bpre-Bpost
    
    #Harvest as a fraction of pre-harvest abundance
    estFracCur=RelHarv/(Bpre)
    result=log(estFracCur)-log(FracCur)
    return(result)
   
def lrootCurFunc(lBpost,   FracCur,HarvestInterval,a,b,xmax,fmax,MinRelAbund,xmid,steep,x0,x1):
    '''Use logistic transform
    '''
    Bpost=expit(lBpost,xmid,steep,x0,x1)
    result=rootCurFunc(Bpost,   FracCur,HarvestInterval,a,b,xmax,fmax,MinRelAbund)
    return(result)

def minFracCur(Bpost,HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund):
    '''
    Evaluate Bpost according to sustainable harvest as a frac of Bre
    '''

    #Calculate the pre-harvest relative biomass
    Bpre=NextRelBmass(OldBiomass=Bpost,OldTime=0,NewTime=HarvestInterval,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=1,CurHarvest=0,minRelBmass=MinRelAbund-.0001,maxRelBmass=1-1e-3)


    FracCur=(Bpre-Bpost)/Bpre
    #Take negative because scipy utilities want to find the minimum
    return(-log(FracCur))

def lminFracCur(lBpost,HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund,xmid,steep,x0,x1):
    '''
    Evaluate Bpost according to sustainable harvest as a frac of Bpre
    On a logistic scale
    '''
    Bpost=expit(lBpost,xmid,steep,x0,x1)
    result=minFracCur(Bpost,HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund)
    return(result)

def MaxFracCur(HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund):
    '''
    Find the post-harvest relative biomass corresponding to the maximum sustainable harvest as a fraction of Bpre
    '''
    
    #If there are lists, do operation for each iteration in markov chain 
    if isinstance(fmax,(list,ndarray)):
        n=len(fmax)
        result=[  MaxRelHarvest(HarvestInterval,SimInterval,fmax[i],a[i],b[i],xmax[i],MinRelAbund[i],MaxRelAbund[i])   for i in range(n)]
        return(result)
    
    #One iteration of markov chain
   
    #Trivial case
    Bpre1=EstPreB(MinRelAbund+0.000,HarvestInterval,a,b,xmax,fmax)
    Bpre2=EstPreB(MinRelAbund+0.002,HarvestInterval,a,b,xmax,fmax)
    if((Bpre1-MinRelAbund)/Bpre1)>((Bpre2-MinRelAbund-.002)/(Bpre2)):
        return(MinRelAbund)
         
    #Values for bracket
    Blow=MinRelAbund-.0001
    Bupp=EstPreB(0.999,-HarvestInterval,a,b,xmax,fmax)
    if Bupp==0.999:
       Bupp=EstPreB(0.99,-HarvestInterval,a,b,xmax,fmax) 
       if Bupp==0.99:
           Bupp=EstPreB(0.98,-HarvestInterval,a,b,xmax,fmax) 
           if Bupp==0.98:
               Bupp=EstPreB(0.97,-HarvestInterval,a,b,xmax,fmax) 
               if Bupp==0.97:
                   Bupp=EstPreB(0.96,-HarvestInterval,a,b,xmax,fmax) 
                   if Bupp==0.96:
                       Bupp=EstPreB(0.95,-HarvestInterval,a,b,xmax,fmax) 
    
    args=(HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund,0.25,Blow,Bupp)
    try:
        Bpost=minimize_scalar(lminFracCur,args=args,bracket=(-1,1)).x
    except:
        print('HarvestAmountUtility 281')
        print(args)
        print(Blow,Bupp)
        print()
        args=(HarvestInterval,SimInterval,fmax,a,b,xmax,MinRelAbund,MaxRelAbund,(Blow+Bupp)/2,0.25,Blow,Bupp,)
        lBpost=minimize_scalar(lminFracCur,args=args,method='Brent',bracket=(-.25,.25),tol=0.001,options={'maxiter':10}).x
        Bpost=expit(lBpost, (Blow+Bupp)/2,0.25,Blow,Bupp)
          
    #Special Cases
    if not(Bpost):
        return(0)
    if Bpost<MinRelAbund:
        return(0)
    if Bpost>0.999:
        return(0.999)
    # Calculation is OK:
    return(Bpost)

def ConstFracCur(FracCur,HarvestInterval,fmax,a,b,xmax,MinRelAbund):
    '''
    Find the post-harvest relative biomass sustained by constant harvest of a fraction of current stock level
    '''
    
    #If there are lists, do operation for each iteration in markov chain 
    if isinstance(fmax,list):
        n=len(fmax)
        result=[  ConstFracCur(FracCur,HarvestInterval,fmax[i],a[i],b[i],xmax[i],MinRelAbund[i])   for i in range(n)]
        return(result)
    
    #One iteration of markov chain
    
    #Does the stock return to the virgin state, Bpre=1?
    Bpre=EstPreB(1-FracCur,HarvestInterval,a,b,xmax,fmax)
    if Bpre>=(0.999):
        return(1-FracCur)
   
    #Check for trivial results 
    if FracCur>((HarvestInterval*fmax)/xmax):
            return(0)

    #Gotta find a root
    x0=logit(max([.001,1-FracCur-.01]),(MinRelAbund+Bpre)/2,.25,MinRelAbund,Bpre)
    x1=logit(min([.999,1-FracCur]),(MinRelAbund+Bpre)/2,.25,MinRelAbund,Bpre)
    xmid,steep,ymin,ymax=(x0+x1)/2,.25,x0,x1
    args=(FracCur,HarvestInterval,fmax,a,b,xmax,1e-3,xmid,steep,ymin,ymax,)
    try:
       PostAbund=root_scalar(lrootCurFunc,x0=x0,x1=x1,args=args).root
    except:
       #Fraction cannot be achieved
       print('HarvestAmountUtility 158')
       print(FracCur,HarvestInterval,fmax,a,b,xmax,MinRelAbund)
       return(0)
    
    #Special cases
    if not(PostAbund):
        return(0)
    if PostAbund<MinRelAbund:
        return(MinRelAbund)
    #Calculations are unexceptional    
    return(PostAbund)

def lConstFracCur(FracCur,HarvestInterval,a,b,xmax,fmax,MinRelAbund):
    '''operate on logit-scale so biomass values remain within range'''
    
    if isinstance(a, list):
        n=len(a)
        result=[ lConstFracCur(FracCur,HarvestInterval,a[i],b[i],xmax[i],fmax[i],MinRelAbund[i])  for i in range(n)]
        return(result)
    
    #Single iteration
    #trivial cases 
    
    
    #Pre-harvest stock level is virgin state
    x0=EstPreB(.9985,-HarvestInterval,a,b,xmax,fmax)#Post-harvest value that will get you to virgin state
    if (0.9985-x0)>FracCur:
       return(0.9985-FracCur)#virgin state achieved
    
    #Post-harvest abundance where the ratio is maximized
    try:
        x1=MaxFracCur(HarvestInterval,HarvestInterval,fmax,a,b,xmax,MinRelAbund,x0)
        Bpre=EstPreB(x1,HarvestInterval,a,b,xmax,fmax)
        maxfrac=(Bpre-x1)/Bpre
        #Targetted fraction cannot be achieved
        if maxfrac<FracCur:
            return(0)
    except:
        x1=MaxFracCur(HarvestInterval,HarvestInterval,fmax,a,b,xmax,MinRelAbund,x0)
        Bpre=EstPreB(x1,HarvestInterval,a,b,xmax,fmax)
        maxfrac=(Bpre-x1)/Bpre
        #Targetted fraction cannot be achieved
        if maxfrac<FracCur:
            return(0)
    
    xmid,steep=(x0+x1)/2,0.25
    args=(FracCur,HarvestInterval,a,b,xmax,fmax,MinRelAbund,xmid,steep,x0,x1,)
 
    try:
        lBpost=root_scalar(lrootCurFunc,args=args, x0=-.25, x1=.25 ).root
        if isinstance(lBpost,complex):
            return(0)
        Bpost=expit(lBpost,xmid,steep,x0,x1)
        if Bpost<MinRelAbund:
            return(0)
        return(Bpost)
    except:
        print('HarvestAmountUtility 413')
        print(args)
        lBpost=root_scalar(lrootCurFunc,args=args, x0=-.25, x1=.25 ).root
        if isinstance(lBpost,complex):
            return(0)
        Bpost=expit(lBpost,xmid,steep,x0,x1)
        if Bpost<MinRelAbund:
            return(0)
        return(Bpost)



def qConstFracCur(FracCur,HarvestInterval,fmax,a,b,xmax,MinRelAbund,q=[.025,.5,.975]):
    '''
    Generate quantiles/CredibleIntervals for a given relative harvest(fraction of current) and harvest interval
    '''
    
    #Get estimates of post-harvest relative abundance 
    PostHarvest=[ BpostFromCurrentFrac(FracCur,HarvestInterval,a[i],b[i],xmax[i],fmax[i],MinRelAbund[i] )  for i in range(len(a)) ]

    
    #Get the quantiles
    qPostHarvest=mquantiles(PostHarvest,prob=q)
    return(qPostHarvest)
    
def curvesConstFracCur(FracCur,HarvestInterval,fmax,a,b,xmax,MinRelAbund,q=[.025,.5,.975]):
    '''
    Generate quantiles/CredibleIntervals for a range of harvest(fraction of current) and  a  harvest interval
    FracCur is a sorted list of values between zero and one
    '''
    result={}
    result['FracCur']=FracCur
    
    #quantiles by relativeHarvest
    byrh=[ qConstFracCur(rh,HarvestInterval,fmax,a,b,xmax,MinRelAbund,q=q)   for rh in FracCur]
    
    #Pick out the quantiles and put them in the results
    for i,q2 in enumerate(q):
        result[q2]=[t[i] for t in byrh]
        
    return(result)

def probConstFracCur(RelHarv,HarvestInterval,fmax,a,b,xmax,RP):
    '''
    RP is a reference point.
    Returns the probability of violating the reference point for a harvest amount expressed as a fraction of virgin biomass
    '''
    
    #Get estimates of post-harvest relative abundance 
    PostHarvest=lConstFracCur(RelHarv,HarvestInterval,fmax,a,b,xmax,RP)
    
    
    #Get the quantiles
    ProbVioate=average([t==0 for t in PostHarvest])
    return (ProbVioate)
    

        
###### when Harvest amount is a linear biomass density.

def TotalCoastLength(SurveyName):
    if SurveyName=='Jervis Inlet':
        result=sum([12010.51386,13315.0372,10300.03298,13248.02983,11127.21882])
        return(result)
    if SurveyName=='Laredo Inlet':
        result=sum([12655.25756 ,9727.8236 ,11994.60722 ,10642.26167 , 11358.97051])
        return(result)
    if SurveyName=='Tolmie Channel':
        result=sum([10819.13718 , 11425.3497,10631.42025 ,10381.82817 , 10561.15054])
        return(result)
    if SurveyName=='Zeballos':
        result=sum([10951.70743 , 11705.87515, 10026.69119, 11523.11191,10226.53932 ])
        return(result)
    return(None)
def TotalVirginBiomass(SiteVirginBiomass):
    '''
    SiteVirginBiomass is a list.  One element for each site.
    Each site-element is another list.  One member per iteration.
    Sum over sites.  One element per iteration.
    '''
    n=len(SiteVirginBiomass[0])
    result=[sum([t[i] for t in SiteVirginBiomass])  for i in range(n)]
    
    return(result)

def medianLinearBiomassDensity(SurveyName,SiteVirginBiomass):
     CL=TotalCoastLength(SurveyName)
     result=[vb/CL  for vb in SiteVirginBiomass]
     return(result)
def randFacLinearBiomassDensity(sdSite,sdTransect,sdYear,rseed=None):
    if(rseed):
        seed(rseed)
    if isinstance(sdSite,list):
        n=len(sdSite)
        result=[randFacLinearBiomassDensity(sdSite[i],sdTransect[i],sdYear[i]) for i in range(n)]
        return(result)
    epsSite=normal(0,sdSite)
    epsYear=normal(0,sdYear)
    result=exp(sdTransect*sdTransect/2+epsSite+epsYear)
    return(result)
def randLinearBiomassDensity(SurveyName,SiteVirginBiomass,sdSite,sdTransect,sdYear,rseed=None):
    
    #Convert VirginBiomass values to linear biomass denities.
    CL=TotalCoastLength(SurveyName)
    medianBiomassDensity=[t/CL for t in SiteVirginBiomass]
    
    #Generate random effects for site and year.  and average over transect-effects
    rfactor=randFacLinearBiomassDensity(sdSite,sdTransect,sdYear,rseed=rseed)
    
    #Combine median values and effects
    n=len(SiteVirginBiomass)
    result=[medianBiomassDensity[i]*rfactor[i] for i in range(n)]
    return(result)
    
def LinearBiomassHarvestToRelHarvest(LinearBiomassHarvest,VirginLinearBiomassDensity):  

    if isinstance(LinearBiomassHarvest,list):
        result=[LinearBiomassHarvestToRelHarvest(LBH,VirginLinearBiomassDensity) for LBH in LinearBiomassHarvest]
        return(result)
    RelHarvest=[LinearBiomassHarvest/t  for t in VirginLinearBiomassDensity]
    return(RelHarvest)
    

def qVariDensityHarvest(RelHarv,HarvestInterval,fmax,a,b,xmax,MinRelAbund,q=[.025,.5,.975]):
    '''
    Generate quantiles/CredibleIntervals for a variable relative harvest and harvest interval
    '''
     
    #There is a list of relative harvest amounts; each element is another list
    if isinstance(RelHarv,list):
        result=[qVariDensityHarvest(rh,HarvestInterval,fmax,a,b,xmax,MinRelAbund,q=q) for rh in  RelHarv]
        return(result)
    
    #Get estimates of post-harvest relative abundance 
   
    n=len(fmax)
    PostHarvest=[lFindBpost(RelHarv[i],HarvestInterval,a[i],b[i],xmax[i],fmax[i],MinRelAbund[i]) for i in range(n)]
    
    #Get the quantiles
    qPostHarvest=mquantiles(PostHarvest,prob=q)
    return(qPostHarvest)
    
def curvesVariDensHarvest(HarvestDensity,VirginBiomassDensity,HarvestInterval,fmax,a,b,xmax,MinRelAbund,q=[.025,.5,.975]):
    '''
    Generate quantiles/CredibleIntervals for a range relative harvest and  a  harvest interval
    HarvestDensity is a list values
        one element per harvest amount
    '''
 
    
    RelHarv=LinearBiomassHarvestToRelHarvest(HarvestDensity,VirginBiomassDensity)
    result={}
    result['HarvestDensity']=HarvestDensity
    
    #Calculate Post Harvest as relative abundance
    n=len(VirginBiomassDensity)
    PostHarvest=[[lFindBpost(rh[i],HarvestInterval,a[i],b[i],xmax[i],fmax[i],MinRelAbund[i])  for i in range(n)]  for rh in RelHarv]
        
    #Take quantiles
    qPostHarvest=[ mquantiles(ph,prob=q) for ph in PostHarvest]
    
    #Pick out the quantiles and put them in the results
    for i,q2 in enumerate(q):
        result[q2]=[t[i] for t in qPostHarvest]
        
    return(result)
    


def probDensHarvest(HarvestDensity,VirginBiomassDensity,HarvestInterval,fmax,a,b,xmax,MinRelAbund,RP):
    '''
    RP is a reference point.
    Returns the probability of violating the reference point for a harvest density expressed as kg/m
    '''
    
    RelHarv=LinearBiomassHarvestToRelHarvest(HarvestDensity,VirginBiomassDensity)
    
    #Calculate Post Harvest as relative abundance
    n=len(VirginBiomassDensity)
    PostHarvest=[lFindBpost(RelHarv[i],HarvestInterval,a[i],b[i],xmax[i],fmax[i],MinRelAbund[i])  for i in range(n)]

    #Get the quantiles
    ProbVioate=average([t==0 for t in PostHarvest]) 
    return (ProbVioate)
    
 
def FracViolateFracCur(RefPoint,FracCur,HarvestInterval,a,b,xmax,fmax,MinRelAbund):
     ''' If harvest is a fraction of current biomass, how probable is a violation of refernece point'''
     Bpost= lConstFracCur(FracCur,HarvestInterval,a,b,xmax,fmax,MinRelAbund)
     violate=[(Bpost[i]<RefPoint[i])   for i in range(len(Bpost))]
     result=average(violate)
     return(result)
    
#############################
def BpostFromHarvestFrac(HarvestAmount,HarvestInterval,a,b,xmax,fmax,MinRelAbund ):
    Bupp=1-HarvestAmount

    MRH=MaxRelHarvest(HarvestInterval,HarvestInterval,fmax,a,b,xmax,MinRelAbund,0.99)
    Bpre=NextRelBmass(OldBiomass=MRH,OldTime=0,NewTime=HarvestInterval,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=1,CurHarvest=0,minRelBmass=MinRelAbund-.0001,maxRelBmass=1-1e-3)
    MaxHarvest=Bpre-MRH
    if MaxHarvest<HarvestAmount:
        return(0)
    Blow=max([.001,MRH   ])
    args=(HarvestAmount,HarvestInterval,a,b,xmax,fmax,Blow,Bupp)
    lBpost=root_scalar(lrootBpostFromHarvestFrac,args=args,x0=-1,x1=1).root
    Bpost=expit(lBpost,(Blow+Bupp)/2,.25,Blow,Bupp)
    if Bpost>=(1-HarvestAmount):
        print('\nHarvestAmountUtility 588')
        print('HarvestAmount,HarvestInterval,a,b,xmax,fmax,MinRelAbund')
        print(HarvestAmount,HarvestInterval,a,b,xmax,fmax,MinRelAbund)
        print('Bupp,Blow,Bpost')
        print(Bupp,Blow,Bpost)
    return(Bpost)


def rootBpostFromHarvestFrac(Bpost, HarvestAmount,HarvestInterval,a,b,xmax,fmax ):
    '''How close does Bpost come to generating the Harvest Amount'''
    Bpre= NextRelBmass(OldBiomass=Bpost,OldTime=0,NewTime=HarvestInterval,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=1,CurHarvest=0)  
    result=(Bpre-Bpost-HarvestAmount)
    return(result)
def lrootBpostFromHarvestFrac(lBpost, HarvestAmount,HarvestInterval,a,b,xmax,fmax,Blow,Bupp ):
    '''on a logit scale'''
    Bpost=expit(lBpost,(Blow+Bupp)/2,.25,Blow,Bupp)
    result=rootBpostFromHarvestFrac(Bpost, HarvestAmount,HarvestInterval,a,b,xmax,fmax )
    return(result)

#############################
def BpostFromCurrentFrac(HarvestFracCur,HarvestInterval,a,b,xmax,fmax,MinRelAbund ):
    Bupp=1-HarvestFracCur
    Bpre=NextRelBmass(OldBiomass=Bupp,OldTime=0,NewTime=HarvestInterval,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=1,CurHarvest=0,minRelBmass=MinRelAbund-.0001,maxRelBmass=1-1e-4)
    if Bpre>0.999:
        return(Bupp)
    MFC=MaxFracCur(HarvestInterval,HarvestInterval,fmax,a,b,xmax,MinRelAbund,0.99)
    Bpre=NextRelBmass(OldBiomass=MFC,OldTime=0,NewTime=HarvestInterval,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=1,CurHarvest=0,minRelBmass=MinRelAbund-.0001,maxRelBmass=1-1e-3)
    MaxFrac=(Bpre-MFC)/Bpre
    if MaxFrac<HarvestFracCur:
        return(0)
    Blow=max([.001,MFC   ])
    args=(HarvestFracCur,HarvestInterval,a,b,xmax,fmax,Blow,Bupp)
    #lBpost=root_scalar(lrootBpostFromCurrentFrac,args=args,x0=logit(.95*Bupp,(Blow+Bupp)/2,.25,Blow,Bupp),x1=logit(1.05*Blow,(Blow+Bupp)/2,.25,Blow,Bupp)).root
    lBpost=root_scalar(lrootBpostFromCurrentFrac,args=args,x0=-1,x1=1).root
    Bpost=expit(lBpost,(Blow+Bupp)/2,.25,Blow,Bupp)
    if Bpost>=(1-HarvestFracCur):
        print('\nHarvestAmountUtility 626')
        print('HarvestAmount,HarvestInterval,a,b,xmax,fmax,MinRelAbund')
        print(HarvestFracCur,HarvestInterval,a,b,xmax,fmax,MinRelAbund)
        print('Bupp,Blow,Bpost')
        print(Bupp,Blow,Bpost)
    return(Bpost)


def rootBpostFromCurrentFrac(Bpost, HarvestFracCur,HarvestInterval,a,b,xmax,fmax ):
    '''How close does Bpost come to generating the Harvest Amount'''
    Bpre= NextRelBmass(OldBiomass=Bpost,OldTime=0,NewTime=HarvestInterval,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=1,CurHarvest=0)   
    result=log( (Bpre-Bpost)/Bpre/HarvestFracCur)
    return(result)
def lrootBpostFromCurrentFrac(lBpost, HarvestFracCur,HarvestInterval,a,b,xmax,fmax,Blow,Bupp ):
    '''on a logit scale'''
    Bpost=expit(lBpost,(Blow+Bupp)/2,.25,Blow,Bupp)
    result=rootBpostFromCurrentFrac(Bpost, HarvestFracCur,HarvestInterval,a,b,xmax,fmax )
    return(result)


    
if __name__ == "__main__":
    RelHarv,HarvestInterval,fmax,a,b,xmax,MinRelAbund=0.05,11,.1,.5,.5,.5,.025
    test1=lFindBpost(RelHarv,HarvestInterval,fmax,a,b,xmax,MinRelAbund)
    print(test1)
    RelHarv,HarvestInterval,fmax,a,b,xmax,MinRelAbund=0.07,1,[.1,.1],[.5,.25],[.25,.5],[2/3,1/3],[.025,.05]
    test2=lFindBpost(RelHarv,HarvestInterval,fmax,a,b,xmax,MinRelAbund)
    print(test2)
    test3=lFindBpost(RelHarv,2,fmax,a,b,xmax,MinRelAbund)
    print(2,test3)
    test4=lFindBpost(RelHarv,3,fmax,a,b,xmax,MinRelAbund)
    print(3,test4)
    test5=lFindBpost(RelHarv,4,fmax,a,b,xmax,MinRelAbund)
    print(4,test5)
    print()
    
    
    FracCur,HarvestInterval,fmax,a,b,xmax,MinRelAbund=0.05,1,.1,.5,.5,.5,.025
    test1=ConstFracCur(FracCur,HarvestInterval,fmax,a,b,xmax,MinRelAbund)
    print(test1)
    FracCur,HarvestInterval,fmax,a,b,xmax,MinRelAbund=0.07,1,[.1,.1],[.5,.25],[.25,.5],[2/3,1/3],[.025,.05]
    test2=ConstFracCur(FracCur,HarvestInterval,fmax,a,b,xmax,MinRelAbund)
    print(test2)
    test3=ConstFracCur(FracCur,2,fmax,a,b,xmax,MinRelAbund)
    print(2,test3)
    test4=ConstFracCur(FracCur,3,fmax,a,b,xmax,MinRelAbund)
    print(3,test4)
    test5=ConstFracCur(FracCur,4,fmax,a,b,xmax,MinRelAbund)
    print(4,test5)
    
    
    
