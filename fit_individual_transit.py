#!/usr/bin/env python
import os
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import batman
import pandas as pd
import emcee
from multiprocessing import Pool
data = pd.read_csv("KOI134_seg.csv", comment='#', header=0)

t = np.array(data["time"])
y = np.array(data["flux"])
yerr = np.array(data["err"])
timingfile = pd.read_csv("transittime3.csv")
epochs = np.array(timingfile["epoch"])
t0s = np.array(timingfile["time"])
t0_errs = np.array(timingfile["err"])
    

def batman_model(time, t0, p, tdur, rprs, mstar, rstar, q1, q2):
    
    aAU = mstar**(1/3)*(p/365.25)**(2/3)  #a in AU
    aors = aAU*215/rstar                  #a over r_star
    b = ((1+rprs)**2.-(tdur/p*np.pi*aors)**2.)**0.5
    
    inc = 180/np.pi*np.arccos(b/aors)
       
    u1 = 2*np.sqrt(q1)*q2
    u2 = np.sqrt(q1)*(1 - (2 * q2))

    params = batman.TransitParams()
    params.t0 = t0                          #time of inferior conjunction
    params.per = p                          #orbital period
    params.rp = rprs                        #planet radius (in units of stellar radii)
    params.a = aors                                             #semi-major axis (in units of stellar radii)
    params.inc = inc                                            #orbital inclination (in degrees)
    params.ecc = 0                                            #eccentricity
    params.w = 90                                              #longitude of periastron (in degrees)
    params.u = [u1,u2]                                          #limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"                              #limb darkening model

    m = batman.TransitModel(params, time)                       #initializes model
    flux = m.light_curve(params)
    return flux


def log_likelihood(theta, time, flux, error):

    t0, p, tdur, rprs, mstar, rstar, q1, q2 = theta
    model = batman_model(time, t0, p, tdur, rprs, mstar, rstar, q1, q2)
    res = flux-model
    z = np.polyfit(time, res, 3)
    trend = np.poly1d(z)(time)
    planet_likelihood  = -0.5*np.sum((flux-model-trend)**2./error**2.) 
    return planet_likelihood 


def log_prior(theta, freeparams, priors):
    log_prior = 0
    t0, p, tdur, rprs, mstar, rstar, q1, q2 = theta
    if (q1<0) or (q1>1):
        return -np.inf
    if (q2<0) or (q2>1):
        return -np.inf
    paramdic = {}
    for k, v in zip(freeparams,theta):
        paramdic[k] = v
    
    for key in priors.keys():
        if not key.endswith("err"):
            log_prior+=-0.5*(paramdic[key]-priors[key])**2./priors[key+"_err"]**2.
    
    
    return log_prior
 
####### Log Probability #######

def log_probability(theta, freeparams, time, flux, error, priors):
    lp = log_prior(theta, freeparams, priors)
    if not np.isfinite(lp):
        return -np.inf
    probability = lp + log_likelihood(theta,time, flux, error)

    return probability



if __name__ == '__main__':

    priors_dic = {}
    priors_dic["mstar"] = 1.4135
    priors_dic["mstar_err"] = 0.13
    priors_dic["rstar"] = 1.7273
    priors_dic["rstar_err"] = 0.244
    priors_dic["q1"] = 0.4040685464
    priors_dic["q1_err"] = 0.045
    priors_dic["q2"] = 0.2739
    priors_dic["q2_err"] = 0.043
    priors_dic["p"] = 67.558406
    priors_dic["p_err"] = 0.2
    priors_dic["rprs"] = 0.06646
    priors_dic["rprs_err"] = 0.0004
    for i in range(len(epochs)):
        t0 = t0s[i]
        priors_dic["t0"] = t0
        priors_dic["t0_err"] = t0_errs[i]
        theta = t0, priors_dic["p"], 11/24., priors_dic["rprs"], priors_dic["mstar"], priors_dic["rstar"], priors_dic["q1"], priors_dic["q2"]
        freeparams = ["t0", "p", "tdur", "rprs", "mstar", "rstar", "q1", "q2"]
        mask = (t<t0+25/24.)*(t>t0-25/24.)
        time = t[mask] 
        flux = y[mask]
        error = yerr[mask]
        nwalkers = (len(theta)*2) + 10
        ndim = len(theta)
        pos = []
        spread_arr = [0.001, 0.001, 1.0, 0.001, 0.01, 0.01, 0.05, 0.05]
        for j in range(nwalkers):
            pos_i = theta + spread_arr * np.random.randn(ndim)		   #different spread for each parameter
            prob = log_probability(pos_i, freeparams, time, flux, error, priors_dic)
            while not np.isfinite(prob):
                pos_i = theta + spread_arr * np.random.randn(ndim)
                prob = log_probability(pos_i, freeparams, time, flux, error, priors_dic)     
            pos.append(pos_i)
 

        with Pool(processes=5) as pool:
            max_n = 75000
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(freeparams, time, flux, error, priors_dic),pool=pool)
            master_pos = np.ones([100*nwalkers,ndim+1])
            counter = 0
                
            for sample in sampler.sample(pos, iterations=max_n, store=False):
                position = sample.coords
                probability=sample.log_prob
                
                for k in range(len(position)):
                    master_pos[counter] = np.array(list(position[k])+[probability[k]])
                    counter += 1
                    
                    if counter >= 100 * nwalkers:
                        np.savetxt("testparams_step", master_pos, fmt='%.10f')
                        os.system("cat testparams_step >> mcmc_tested_params_%d" % epochs[i])
                        counter = 0

