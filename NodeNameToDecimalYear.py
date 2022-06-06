from numpy import array,ndarray
from datetime import datetime as dt
from ToDecimalYear import ToDecimalYear,FromDecimalYear

def NodeNameToDecimalYear(NName,NumLead=2,sep='_'):
    '''NName is string.  It is broken up according to occurences of sep.
    The first NumLead substrings are ignored.
    The remaining substrings are converted to strings and then to a decimal year'''
    
    x=NName.split(sep=sep)[NumLead:]
    ymd=[int(i) for i in x][:3]
    result=ToDecimalYear(ymd)
    return(result)

def DecimalYearFromNode(node,NumLead=2,sep='_'):
    if isinstance(node,(list,ndarray)):
        result=[DecimalYearFromNode(n,NumLead=2,sep=sep) for n in node     ]
        return(result)
    result=NodeNameToDecimalYear(node.__name__,NumLead=NumLead,sep=sep)
    return(result)
    
def SSQYearFromNode(node,reftime,NumLead=2,sep='_'):
    dtreftime=ToDecimalYear(reftime)
    if isinstance(node,(list,ndarray)):
        result=[SSQYearFromNode(n,dtreftime,NumLead=NumLead,sep=sep) for n in node     ]
        return(result)
    curdy=DecimalYearFromNode(node,NumLead=2,sep=sep)
    SSQ=(dtreftime-curdy)**2
    result={'node':node,'curdy':curdy,'SSQ':SSQ}
    return(result)
def MatchNodeRefTime(nodelist,reftime,NumLead=2,sep='_'):
    reftime2=FromDecimalYear(reftime)
    listSSQ=SSQYearFromNode(nodelist,reftime2,NumLead=NumLead,sep=sep)
    try:
            minSSQ=min(listSSQ,key=lambda t:t['SSQ'])
    except:
            print('NodeNameToDecimalYear 37 ')
            print(reftime)
            print(NumLead)
            print(nodelist)
            print(listSSQ)
            print()
            minSSQ=min(listSSQ,key=lambda t:t['SSQ'])
    return(minSSQ['node'])
