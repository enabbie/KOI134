import rebound

import numpy as np

import pandas as pd
from scipy.optimize import minimize
from matplotlib import pyplot as plt
import os
import itertools
from multiprocessing import Pool
import time
import tempfile
import sys
import emcee

os.environ['OMP_NUM_THREADS']="1"
filename = str(sys.argv[1])								#name for final file results will be saved to -> corresponds to the date of the run and period ratio
pratio_name = str(sys.argv[2])


stellar_mass = 1.413524
stellar_radius = 1.727299
Mjup_Msun= 0.000954588

# inclination of the transiting planet (take from the transit fit)

df1 = pd.read_csv("transittime3.csv")
time1 = np.array(df1.time)
epoch1 = np.array(df1.epoch)
err1 = np.array(df1.err)

# compute the linear period of the transiting planet to have a start point for the grid search
fit0 = np.polyfit(epoch1, time1, 1)
linear_P1 = fit0[0]
res1 = 24*60*(time1-fit0[1]-epoch1*linear_P1)

df2 = pd.read_csv("transitduration_fit.csv", delim_whitespace=True)

tdur_obs = np.array(df2["duration"])
tdur_err = np.array(df2["duration_error"])


# maximum number of transits we are computing
N1 = 23 
# maximum number of iteration during bisector
maxiter = 10000

def log_likelihood(theta):
    # the likelihood function using n-body integration
    try:
        # unpack the parameters
        logm1, logm2, p1, p2p1, sesinw1, secosw1, sesinw2, secosw2,Omega1,Omega2,inc2,M1,M2,b1 = theta
        p2 = p2p1*p1	
        aAU = stellar_mass**(1/3)*(p1/365.25)**(2/3)	#a in AU
        aors = aAU*215/stellar_radius

        inc1 = 180/np.pi*np.arccos(b1/aors)

        m1 = 10**(logm1)*1.41/1.17
        m2 = 10**(logm2)*1.41/1.17
        a1 = (p1/365.25)**(2./3.)*(stellar_mass)**(1./3.) 
        a2 = (p2/365.25)**(2./3.)*(stellar_mass)**(1./3.) 
        e1 = sesinw1**2+secosw1**2
        e2 = sesinw2**2+secosw2**2
        w1 = np.arctan2(sesinw1, secosw1)
        w2 = np.arctan2(sesinw2, secosw2)
        Omega1 = Omega1/180.*np.pi
        Omega2 = Omega2/180.*np.pi
        M2 = M2/180.*np.pi
        M1 = M1/180.*np.pi
        
        # initialize the n-body simulation
        sim = rebound.Simulation()

        sim.integrator="ias15"
        sim.add(m=stellar_mass)
        inc1 = (90-inc1)/180.*np.pi 
        inc2 = (90-inc2)/180.*np.pi 

        sim.add(m=m1*Mjup_Msun, a=a1,e=e1, inc = inc1,Omega=Omega1,omega=w1,M=M1)
        sim.add(m=m2*Mjup_Msun, a=a2, e=e2, inc= inc2, Omega=Omega2, omega=w2,M=M2)
        # exit when hit 3 hill radius
        sim.exit_min_distance = 3*((m1+m2)/3.*Mjup_Msun)**(1./3.)*((a1+a2)/2) 
        sim.move_to_com()
        # integrate until N1 transits happened, and record the mid transit times
        transittimes1 = np.zeros(N1)
        durarray = np.zeros(N1)
        dt = p1/365.25*2.*np.pi/10. #0.001
        p = sim.particles
        i = 0
        try:
            while i<N1:
                x_old = p[1].x - p[0].x  # (Thanks to David Martin for pointing out a bug in this line!)
                t_old = sim.t
                sim.integrate(sim.t+dt) # check for transits every dt time units. Note that 0.01 is shorter than one orbit
                t_new = sim.t
                niter = 0
                if x_old*(p[1].x-p[0].x)<0. and p[1].z-p[0].z>0.:		# sign changed (y_old*y<0), planet in front of star (x>0)
                    while t_new-t_old>1e-5:   # bisect until prec of 1e-5 reached
                        if x_old*(p[1].x-p[0].x)<0.:
                            t_new = sim.t
                        else:
                            t_old = sim.t
                        sim.integrate( (t_new+t_old)/2.)
                        niter+=1
                        if niter>maxiter:
                            raise RuntimeError("this is too many iterations")
                    transittimes1[i] = sim.t
                    o = sim.particles[1].calculate_orbit(primary=sim.particles[0]) 
                    ar = o.a*215/stellar_radius
        
                    sep = (1-o.e*np.cos(o.E))#*215/stellar_radius
                    inc1_i = (np.pi/2.-o.inc)
                    b = ar*sep*np.cos(inc1_i)
                    temp = 1-ar**2.*sep**2.*np.cos(inc1_i)**2.
                    dur_i = (o.P/2./np.pi*365.25)/np.pi*(sep**2./np.sqrt(1-o.e**2.))*np.arcsin(np.sqrt(temp)/(ar*sep*np.sin(inc1_i)))
                    durarray[i] = dur_i

                    i += 1
                    sim.integrate(sim.t+dt/10.)				  # integrate 0.001 to be past the transit
        except rebound.Encounter as error:
            print("not stable")
            return -np.inf
        
        # convert time from rebound units to dates and compute residual against observed linear time and epoch

        res_p1_reb = (transittimes1)*(365.25/2./np.pi)-linear_P1*np.array(range(N1))
        # get rid of gapped transits
        res_p1_reb = res_p1_reb[np.in1d(range(N1),epoch1)]
        tdur_reb = durarray[np.in1d(range(N1),epoch1)]

        # convert units to minutes (this is not nesessary, just keep the units to be the same as the err
        res_p1_reb*=60.*24.

        # make sure we offset any constants in the timing of the first transit uning a leastsq
        def func1(x):
            return np.sum((res_p1_reb-res1-x)**2./(err1*24*60)**2.)
        x0 = [0]
        res = minimize(func1, x0)
        baseline1 = res.x

        # the final loglikelihood
        logprob_ttv = func1(baseline1)
        logprob_tdv = np.sum((tdur_reb-tdur_obs)**2./tdur_err**2.)
        logprob = logprob_ttv+logprob_tdv

        print(theta, -0.5*logprob, -0.5*logprob_ttv, -0.5*logprob_tdv)
        return -0.5*logprob   
    
    except RuntimeError:
        return -np.inf
counter = 0

def log_prior(theta):
    logm1, logm2, p1, p2p1, sesinw1, secosw1, sesinw2, secosw2,Omega1,Omega2,inc2,M1,M2, b1 = theta
    log_prior = 0
    p2 = p2p1*p1	

    w1 = np.arctan2(sesinw1,secosw1) * 180 / np.pi
    w2 = np.arctan2(sesinw2,secosw2) * 180 / np.pi
    
    e1 = sesinw1**2 + secosw1**2 
    e2 = sesinw2**2 + secosw2**2
    if inc2>90:
        return -np.inf
    if M2<-180. or M2>180.:
        return -np.inf
    if M1<-180. or M1>180.:
        return -np.inf
    if b1<0 or b1>0.9:
        return -np.inf
    #box priors
    if (sesinw1 > 1) or (sesinw1 < -1) or (sesinw2 > 1) or (sesinw2 < -1):
        #print("sesinw out of bounds")
        return -np.inf
    
    if (logm1 > np.log10(10)) or (logm1 < np.log10(0.5)) or (logm2 > np.log10(10)) or (logm2 < np.log10(.00314558)):	#masses of both planets are bounded between 1 earth mass and 1 jupiter mass
    	#print("mass out of bound")
        return -np.inf
    
    #bounding angles
    if (Omega1 < 0) or (Omega1 > 180.) or (Omega2 < -360.) or (Omega2 > 360.):		#Omega between 0 and pi
        #print("Omega out of bounds")
        return -np.inf
    
    if (w1 > 180) or (w1 < -180.) or (w2 > 360) or (w2 < 0):					#w between 0 and 2pi
        #print("w out of bound")
        return -np.inf
    
    #bounding eccentricity
    if (e1 > 0.8) or (e1 < 0) or (e2 > 0.8) or (e2 < 0):
        #print("e out of bound")
        return -np.inf
    #bounding p1
    if (p1 > 68) or (p1 < 67):
        #print("p1 out of bounds")
        return -np.inf
    
    #bounding period ratio
    pratio = p2/p1
    
    if (pratio_name == '13') :
        if (pratio < .28) or (pratio > .38):   #13 = 1:3, 12 = 1:2, and so on
            return -np.inf
        else:
            pass
    
    elif (pratio_name == '12') :
        if (pratio < .4) or (pratio > .6):
            #print("pratio out of bound")
            return -np.inf
        else:
            pass
    
    elif (pratio_name == '23') :
        if (pratio < .62) or (pratio > .7):
            return -np.inf
        else:
            pass
    
    elif (pratio_name == '34') :
        if (pratio < .71) or (pratio > .79):
            return -np.inf
        else:
            pass
    
    elif (pratio_name == '43') :
        if (pratio < 1.3) or (pratio > 1.38):
            return -np.inf
        else:
            pass
    
    elif (pratio_name == '32') :
        if (pratio < 1.45) or (pratio > 1.55):
            return -np.inf
        else:
            pass
    
    elif (pratio_name == '21') :
        if (pratio < 1.8) or (pratio > 2.2):
            return -np.inf
        else:
            pass
    
    elif (pratio_name == '31') :
        if (pratio < 2.8) or (pratio > 3.2):
            return -np.inf
        else:
            pass
    
    return log_prior

def log_probability(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    if np.isnan(log_likelihood(theta)):
        return -np.inf
    return lp + log_likelihood(theta)

if __name__ == '__main__':
    theta_list = []
    with Pool(processes=60) as pool:

        initial_loc = pd.read_csv('bestfit_tdv_Jul17nobias.csv')
        theta_test = np.array(initial_loc.iloc[0,:])   #takes only the values of the variables (make sure they're in the proper order (m1, m2, p1, p2, sesinw1, secosw1, sesinw2, secosw2, Omega1, Omega2)
        print(theta_test)
        nwalkers = (len(theta_test)*2) + 10
        theta_test[3]/=theta_test[2] 
        ndim = len(theta_test)
        print(nwalkers)
        print(log_likelihood(theta_test))

        spread_arr = [.1,.1,.003,.003,.5,.5,.5,.5,.05,.05,.05, 0.2, 0.2, 0.1]
        spread_arr = np.array(spread_arr)/10.

        pos = []
        for i in range(nwalkers):
            pos_i = theta_test + spread_arr * np.random.randn(ndim)		   #different spread for each parameter
            prob = log_probability(pos_i)
            print(pos_i)
            print(prob)
            while not np.isfinite(prob):
                pos_i = theta_test + spread_arr * np.random.randn(ndim)
                prob = log_probability(pos_i)     

            pos.append(pos_i)
            print(np.array(pos).shape)
        max_n = 150000
            
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, pool=pool)
        os.system("rm -rf rebound_mcmc_tdv_tested_params")
        os.system("rm -rf rebound_tdv_testparams_step")
        master_pos = np.ones([100*nwalkers,ndim+1])
        counter = 0
        pos = np.array(pos)
        print(nwalkers)
        print(ndim)
        print(pos.shape)
        for sample in sampler.sample(pos, iterations=max_n, store=False):
            position = sample.coords
            probability=sample.log_prob
            for i in range(len(position)):
                master_pos[counter] = np.array(list(position[i])+[probability[i]])
                counter += 1
                if counter >= 100 * nwalkers:
                    tf = tempfile.NamedTemporaryFile()
                    np.savetxt(tf.name, master_pos, fmt='%.10f')
                    os.system(f'cat {tf.name} >> rebound_mcmc_tdv_tested_params_'+filename)
                    counter = 0
