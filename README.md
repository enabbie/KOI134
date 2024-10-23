This is a README file with instructions on how to run the software associated with my analysis of the KOI-134 system.

# System requirements
Software dependencies: Please see the "python environments" folder for version numbers of all packages associated with the Python environments used to run the REBOUND and REBOUNDx simulations. Generally, the software runs on Python 3.8 or higher, and requires installation of the batman, isochrones, and REBOUND/REBOUNDx packages.

This software has been tested on Linux and MacOS Monterey 12.5 operating systems. Using a Linux environment or HPC is preferred.

There is no non-standard hardware required, although this code can make use of multiprocessing.

# Installation Guide
1. Install Dependencies
   - [batman](https://lkreidberg.github.io/batman/docs/html/installation.html)
   - [Corner](https://corner.readthedocs.io/en/latest/install/)
   - [Emcee](https://emcee.readthedocs.io/en/stable/user/install/)
   - [MultiNest](https://github.com/JohannesBuchner/MultiNest)
   - [Isochrones](https://isochrones.readthedocs.io/en/latest/install.html)
   - [Astroquery](https://astroquery.readthedocs.io/en/latest/#installation)
   - [REBOUND](https://rebound.readthedocs.io/en/latest/quickstart_installation/)
   - [REBOUNDx](https://reboundx.readthedocs.io/en/latest/python_quickstart.html)
3. Set up Python environment
   Define the python environments using those listed in the "python environments" folder.
4. Load environments. The transit fit MCMC requires MultiNest, which must be loaded before running the script (copy/paste the code block below, and modify the path locations).
   ```
   module load multinest
   export PATH=/path/to/py38/bin/:$PATH
   export PYTHONPATH=/path/to/envs/py38/lib/python3.8/site-packages/:$PYTHONPATH

   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/.local/lib
   ```

Typical install time: 30 minutes.

# Demo
To use this software on data, download the code and the relevant config ("prior") files. Once the dependencies are installed, simply run the python scripts in the following order:
  1. Batman light curve MCMC, 'multiplanet_fitting_script.py' (this requires an extra argument for the name of the prior file you are using)
  2. REBOUND grid search (must have an extra argument for the resonance you're probing, in the format '21' for 2:1, '12' for 1:2, etc.)
  3. REBOUND MCMC
  4. REBOUND integration/stability check

The expected outputs are .txt or .csv files with (1) best-fit transit centers at each epoch, (2) best-fit system parameters for a given period resonance, (3) best-fit orbital elements and planet parameters, and (4) planet orbital elements over a timescale of 10 Myr.

On a "normal" desktop computer, this demo may take several days. The integrations and MCMC algorithms are computationally expensive, especially without multiprocessing.

# Instructions for Use
To run the software on data, please follow the steps outlined above. 

The definitions for each file and its functionality are as follows:

- 'bestfit_plotting.py': compiles the best-fit light curve parameters after the light curve MCMC.
- 'bestfit_tdv_Jul17nobias.csv': CSV with best-fit parameters from the REBOUND TTV/TDV joint fit.
- 'fit_individual_transit.py': light curve MCMC where the transit durations are left as a free parameter (used in auxiliary light curve analysis).
- 'gridsearch_rebound_ttv.py': REBOUND grid search code that produces optimal starting points for the REBOUND MCMC (this must be initialized at a specific resonance)
- 'KOI134_seg.csv': the Kepler light curve file for KOI-134.
- 'koi134_prior_notdvs.csv': input CSV for main light curve MCMC fit.
- 'multiplanet_fitting_script.py': script for main light curve MCMC fit.
- 'rebound_int_gr.py': REBOUNDx integration code, which takes into account General Relativity effects.
- 'rebound_mcmc_tdv_v5_nob_ias.py': script for REBOUND TTV/TDV joint fit.
- 'transitduration_fit.csv': CSV with best-fit transit durations.
- 'transittime3.csv': CSV of initial estimates for transit centers, which is used as an input for the light curve MCMC analyses.
