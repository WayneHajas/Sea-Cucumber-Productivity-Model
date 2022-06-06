from datetime import datetime
from numpy import ndarray,append


from pymc import  *
from pylab import *
from tables import *
import tables
import warnings
warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)

        
def ExecuteModel(Model,Niter=11000,burn=1000,thin=10,name='MCMC',db='hdf5',verbose=2,maxSimplex=0,maxPowell=0,maxBFGS=0):

    N=MAP(Model)
    if maxSimplex>0:
       print('start Simplex ')
       N.fit(iterlim=maxSimplex,tol=1.0,verbose=verbose)
    if maxPowell>0:
       print('start powell ')
       N.fit('fmin_powell',iterlim=maxPowell,verbose=verbose)
    if maxBFGS>0:
       print('start bfgs ')
       N.fit('fmin_l_bfgs_b',iterlim=maxBFGS,verbose=verbose)
       
    # Go through the fit a second time   
    if maxSimplex>0:
       print('start Simplex ')
       N.fit(iterlim=maxSimplex,tol=1.0,verbose=verbose)
    if maxPowell>0:
       print('start powell ')
       N.fit('fmin_powell',iterlim=maxPowell,verbose=verbose)
    if maxBFGS>0:
       print('start bfgs ')
       N.fit('fmin_l_bfgs_b',iterlim=maxBFGS,verbose=verbose)

    
    # Create MCMC object
    M = MCMC(Model,db=db,name=name)

    # Sample
    M.sample(Niter,burn=burn,thin=thin,verbose=verbose)


    #------------------------Export results
    #DIC
    DIC=M.dic
    DICfile=open("dic.txt",'w')
    DICfile.write('%f' %DIC)
    DICfile.close()

    # Create an object to contain output statistics
    results=M.stats()
    # Get paramter names
    parms=results.keys()
    # Get list of summary statistics
    stats=results[parms[0]].keys()

    # Create output file for summary statistics
    output='output.csv'
    h=open(output,'w')
    # Add column headings from names of summary statistics
    h.write(('%s,'*5+'%s\n') %('parameter','mean', 'median', 'sd', '0.025', '0.975'))
    # Interate over output quantities of interest
    for parm in parms:
        h.write(('%s,'*5+'%s\n') %(str(parm),str(results[parm]['mean']), str(results[parm]['quantiles'][50]), str(results[parm]['standard deviation']), str(results[parm]['quantiles'][2.5]), str(results[parm]['quantiles'][97.5])))
    h.close()
        
def Restart(Model,db,Niter=11000,burn=1000,thin=10,verbose=2):
    M=MCMC(Model,db=db)
    M.sample(Niter,burn=burn,thin=thin,verbose=verbose)

       
if __name__ == "__main__":

    from AllQuadTS import AllQuadTS
    import os
    dir = os.path.dirname(__file__)
    path=mdbfile=os.path.join(dir, '..\data\ShowCount.mdb')    
    SurveyTitle='Marina'
    Site=7
    aqt=AllQuadTS( mdbfile,SurveyTitle,Site)
    test=SFModel(aqt,a=1.01,b=1.01)
