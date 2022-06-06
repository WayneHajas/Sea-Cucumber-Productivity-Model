import matplotlib.pyplot as plt



def HistoByHalf(plt,x,ParamName=None, LowCol='r',UppCol='b',alpha=0.25):
    
    i=[j for j, t in enumerate(x)]
    n=len(i)
    
    plt.subplot(2,1,1)
    plt.plot(i,x,'k.')
    if ParamName:
        plt.title(ParamName)
    
    xmin=min(x)
    xmax=max(x)
    bins=[ xmin+(xmax-xmin)*j/100 for j in range(101)    ]
    x1=x[:int(n/2)]
    x2=x[-int(n/2):]
    
    plt.subplot(2,1,2)
    plt.hist(x1,bins=bins,color=LowCol,alpha=alpha)
    plt.hist(x2,bins=bins,color=UppCol,alpha=alpha)
    


    
if __name__ == "__main__":
    x=list(range(10))
    prob=[.025,.975]
    
    plt.close()
    HistoByHalf(plt,x,ParamName='test', LowCol='r',UppCol='r',alpha=0.25)
    plt.show()