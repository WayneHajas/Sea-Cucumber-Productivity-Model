# Sea-Cucumber-Productivity-Model
Python code to implement the productivity model used in the CSAS Sea Cucumber paper(will provide reference when it is published).

The model was implemented under pymc2(https://github.com/pymc-devs/pymc2).

Markov chains were written as hdf5 files.  The library to read data from these files is GetParamStats.py

The .py files in this repository are sufficient for a Bayesian implementation of the model.  Other files are useful for post-analyses.  There are files that were created for post-analyses that have since been abandoned.  This repository does not conatain all the figures and tables in the report.  Should be cleaned up some day.


An example of code to actually run the model is Jervis.seed.20180824.py
