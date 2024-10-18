This is a README file with instructions on how to run the software associated with my analysis of the KOI-134 system.

# System requirements
Software dependencies: Please see the "python environments" folder for version numbers of all packages associated with the Python environments used to run the REBOUND and REBOUNDx simulations. Generally, the software runs on Python 3.8 or higher, and requires installation of the batman, isochrones, and REBOUND/REBOUNDx packages.

There is no non-standard hardware required, although this code can make use of multiprocessing.

# Installation Guide
1. Install Dependencies
   [will put install links here]
2. Set up Python environment
   Define the python environments using those listed in the "python environments" folder.
3. [insert export expression]

Typical install time: 30 minutes.

# Demo
To use this software on data, download the code and the relevant config ("prior") files. Once the dependencies are installed, simply run the python scripts in the following order:
  1. Batman light curve MCMC
  2. REBOUND grid search (must have an extra argument for the resonance you're probing, in the format '21' for 2:1, '12' for 1:2, etc.)
  3. REBOUND MCMC
  4. REBOUND integration/stability check

The expected outputs are .txt or .csv files with (1) best-fit transit centers at each epoch, (2) best-fit system parameters for a given period resonance, (3) best-fit orbital elements and planet parameters, and (4) planet orbital elements over a timescale of 10 Myr.

On a "normal" desktop computer, this demo may take several days. The integrations and MCMC algorithms are computationally expensive, especially without multiprocessing.

# Instructions for Use
