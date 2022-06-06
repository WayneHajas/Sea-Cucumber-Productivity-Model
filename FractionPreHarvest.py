'''
Convert the estimated total allowable catch to a fraction of pre-harvest abundance

This would likely be better done as a database operation but I am just trying to minimize the number of software tools.
'''

import csv


def FPH(TargetReferenceCSV,AllowableCatchCSV,FractionPreHarvestCSV):
    '''
    All the files have the same structure.
        A header-row to indicate the z-value.
        One column per z-value
        All files must give abundance in same units; fraction of virgin biomass, linear biomass density or spatial biomass density
    There is the assumption that the ith row in TargetReferenceCSV and the ith row in AllowableCatchCSV will correspond to the same iteration within a markov chain
    '''
    
    #Open file with reference points.  These are post-harvest values.
    with open(TargetReferenceCSV,newline='\n') as TRCSV:
        TRreader=csv.reader(TRCSV)
        ColumnNames=next(TRreader)
        
        #Open file with Allowable catch at reference point
        with open(AllowableCatchCSV,newline='\n') as ACCSV:
            ACreader=csv.reader(ACCSV)
            ColumnNames=next(ACreader)
            
                  
            #Open the results file.  Allowable catch as a fraction of the pre-harvest abundance.
            writefile=open(FractionPreHarvestCSV,'w')
            with writefile:
                writer=csv.writer(writefile,lineterminator='\n')
                writer.writerow(ColumnNames)
                
                for t in TRreader:
                    s=next(ACreader)
                    CurResult=[]
                    for i, r in enumerate(t):
                        CurResult+=[float(s[i])/(float(r)+float(s[i]))]
                    writer.writerow(CurResult)
                    
if __name__ == "__main__":

    TargetReferenceCSV='USRbySDYear.csv'
    AllowableCatchCSV='AC_RelBM_tequal1.csv'
    FractionPreHarvestCSV='FractionPreHarvest_tequal1.csv'