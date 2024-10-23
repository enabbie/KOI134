from multiplanet_fitting_script import batman_model, Variable, Params, read_fit_param_csv, phasefold, bin_data, log_probability
import numpy as np
import matplotlib.pyplot as plt
import batman
import pandas as pd
import corner
import emcee
from IPython.display import display, Math

out_folder = '/path/to/fit/results/directory/'

data = pd.read_csv("KOI134_seg.csv", header=0)  #light curve
t = np.array(data['time'])
flux = np.array(data["flux"])
err = np.array(data["err"]) #error


#read in chains
#reader = emcee.backends.HDFBackend('chains.h5')

chains = np.loadtxt(out_folder+'mcmc_tested_params',skiprows=1)
chain_df = pd.DataFrame(chains).tail(10000)

#plotting chains
system_list, freeparams, fixed_params, true_values = read_fit_param_csv('/path/to/koi134_prior_notdvs.csv')

ndim = len(freeparams.unpack())
fig, axes = plt.subplots(ndim, figsize=(ndim, 10), sharex=True)
#samples = reader.get_chain()

labels = freeparams.unpack()

#############################################################
###########         Obtain results, then          ###########
########### re-pack and store parameters in their ###########
###########      proper planet dictionaries       ###########
#############################################################

#displaying best fit parameters
bestfit_df = pd.DataFrame(columns=['name','value','lower_e','upper_e','1sigma'])

bestfit_values = []
upper_err = []
lower_err = []
sig_arr = []
epochs = []
centers = []
t_err = []
objno_arr = []

for i in range(ndim):
    #mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    mcmc = np.percentile(chain_df[i], [16, 50, 84])
    q = np.diff(mcmc)
    txt = "\mathrm{{{3}}} = {0:.7f}_{{-{1:.7f}}}^{{{2:.7f}}}"
    txt = txt.format(mcmc[1], q[0], q[1], labels[i])

    bestfit_values.append(mcmc[1])
    upper_err.append(q[1])
    lower_err.append(q[0])

    onesig_err = np.sqrt(q[0]*q[1])
    sig_arr.append(onesig_err)


    if labels[i][:2] == 't0':
         word = labels[i]
         var_name = word.split('_')[0]                   #omits everything after '_' -> ex. t01_3 becomes t01
         transit_number = int(eval(word.split('_')[1]))       #saves transit number as an integer (so we still know, for example, that t01_3 corresponds to transit 3)
         var_subscript = int(var_name[-1])                    #do this to get the object number

         centers.append(mcmc[1])
         epochs.append(transit_number)
         t_err.append(onesig_err)
         objno_arr.append(var_subscript)
         

bestfit_df["name"] = labels
bestfit_df["value"] = bestfit_values
bestfit_df["lower_e"] = lower_err
bestfit_df["upper_e"] = upper_err
bestfit_df["1sigma"] = sig_arr       #1 sigma error

bestfit_df.to_csv(out_folder+'bestfit_values.csv')

transittimes_df = pd.DataFrame(columns=['objectno','epoch','time','err'])
transittimes_df['objectno'] = objno_arr
transittimes_df['epoch'] = epochs
transittimes_df['time'] = centers
transittimes_df['err'] = t_err

transittimes_df.to_csv(out_folder+'fitted_transittimes.csv')
