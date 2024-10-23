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
import scipy

os.environ['OMP_NUM_THREADS']="1"

bestfit_df = pd.read_csv('/path/to/bestfit_values.csv')

# host star property and constants
stellar_mass = bestfit_df['value'][bestfit_df['name']=='mstar'].item()
stellar_radius = bestfit_df['value'][bestfit_df['name']=='rstar'].item()
Mjup_Msun= 0.000954588

# inclination of the transiting planet (take from the transit fit)
b = bestfit_df['value'][bestfit_df['name']=='b1'].item()
aAU = stellar_mass**(1/3)*(bestfit_df['value'][bestfit_df['name']=='p1'].item()/365.25)**(2/3)  #a in AU
aors = aAU*215/stellar_radius

inc1 = 180/np.pi*np.arccos(b/aors)

# inclination of the non-transiting planet, fixed for the grid search
inc2 = 87.

# read in the transit times we are fitting for
df1 = pd.read_csv("/path/to/fitted_transittimes.csv")
time1 = np.array(df1.time)
epoch1 = np.array(df1.epoch)
err1 = np.array(df1.err)

# compute the linear period of the transiting planet to have a start point for the grid search
fit0 = np.polyfit(epoch1, time1, 1)
linear_P1 = fit0[0]
res1 = 24*60*(time1-fit0[1]-epoch1*linear_P1)



# maximum number of transits we are computing
N1 = 23 
# maximum number of iteration during bisector
maxiter = 10000

def polynomial(t, c0, c1, c2, c3):
    return c3*t**3 + c2*t**2 + c1*t + c0

def normalize(t):
    return (t - np.nanmean(t))/len(t)

def least_sq(theta_p, residual, lc_err):
    c0, c1, c2, c3 = theta_p[0], theta_p[1], theta_p[2], theta_p[3]
    model = polynomial(normalize(t), c0,c1,c2,c3)

    return (residual - model) / lc_err

def log_likelihood(theta):
    # the likelihood function using n-body integration
    try:
        # unpack the parameters
        logm1, logm2, p1, p2, sesinw1, secosw1, sesinw2, secosw2,Omega2,M2 = theta
        m1 = 10**(logm1)
        m2 = 10**(logm2)
        a1 = (p1/365.25)**(2./3.)*(stellar_mass)**(1./3.) 
        a2 = (p2/365.25)**(2./3.)*(stellar_mass)**(1./3.) 
        e1 = sesinw1**2+secosw1**2
        e2 = sesinw2**2+secosw2**2
        w1 = np.arctan2(sesinw1, secosw1)
        w2 = np.arctan2(sesinw2, secosw2)

        # initialize the n-body simulation
        sim = rebound.Simulation()
        sim.integrator="whfast"
        sim.add(m=stellar_mass)

        sim.add(m=m1*Mjup_Msun, a=a1,e=e1, inc =(90-inc1)/180.*np.pi,Omega=0,omega=w1,M=0)
        sim.add(m=m2*Mjup_Msun, a=a2, e=e2, inc=(90-inc2)/180.*np.pi, Omega=Omega2, omega=w2,M=M2)

        sim.exit_min_distance = 3*((m1+m2)/3.*Mjup_Msun)**(1./3.)*((a1+a2)/2) 
        sim.move_to_com()


        # integrate until N1 transits happened, and record the mid transit times
        transittimes1 = np.zeros(N1)
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
                if x_old*(p[1].x-p[0].x)<0. and p[1].z-p[0].z>0.:   # sign changed (y_old*y<0), planet in front of star (x>0)
                    while t_new-t_old>1e-5:   # bisect until prec of 1e-5 reached
                        if x_old*(p[1].x-p[0].x)<0.:
                            t_new = sim.t
                        else:
                            t_old = sim.t
                        sim.integrate( (t_new+t_old)/2.)
                        niter+=1
                        if niter>maxiter:
                            raise RuntimeError("this is too many iterations")
                        #print('#',niter,t_new-t_old)
                    transittimes1[i] = sim.t
                    i += 1
                    sim.integrate(sim.t+dt/10.)       # integrate 0.001 to be past the transit
        except rebound.Encounter as error:
            return -np.inf
 
        # convert time from rebound units to dates and compute residual against observed linear time and epoch
        res_p1_reb = (transittimes1-np.nanmin(transittimes1))*(365.25/2./np.pi)-linear_P1*np.array(range(N1))
    
        # get rid of gapped transits
        res_p1_reb = res_p1_reb[np.in1d(range(N1),epoch1)]

        # convert units to minutes
        res_p1_reb*=60.*24.
        # make sure we offset any constants in the timing of the first transit uning a leastsq
        def func1(x):
            return np.sum((res_p1_reb-res1-x)**2./(err1*24*60)**2.)
        x0 = [0]
        res = minimize(func1, x0)
        baseline1 = res.x
    
        # the final loglikelihood
        logprob = func1(baseline1)
        results = np.array([m1, m2, p1, p2, sesinw1, secosw1, sesinw2, secosw2,Omega2,M2, -0.5*logprob])

        tf = tempfile.NamedTemporaryFile()
        np.savetxt(tf.name, results.reshape(1,results.shape[0]), fmt='%.10f')                                                         #saves each point as a temporary file
        os.system(f"cat {tf.name} >> /path/to/gridsearch/rebound_tested_params")              
 
        return m1, m2, p1, p2, sesinw1, secosw1, sesinw2, secosw2,Omega2,M2,-0.5*logprob 
   
    except RuntimeError:
        results = np.array([m1, m2, p1, p2, sesinw1, secosw1, sesinw2, secosw2,Omega2, M2,-np.inf])
        tf = tempfile.NamedTemporaryFile()
        np.savetxt(tf.name, results.reshape(1,results.shape[0]), fmt='%.10f')

        os.system(f"cat {tf.name} >> /path/to/gridsearch/rebound_tested_params")              

        return m1, m2, p1, p2, sesinw1, secosw1, sesinw2, secosw2,Omega2, M2, -np.inf
counter = 0

if __name__ == '__main__':
    
    p1_array = np.linspace(linear_P1-0.5, linear_P1+0.5, 15)
    
    #### initialize the period ratio grid ####

    # search for 1:2, 2:3, 3:4, 1:3, 4:3, 3:2, 2:1, 3:1
    # create small lists for each resonance
    p12_arr = np.linspace(0.45,0.55,55)  # 1:2
    p23_arr = np.linspace(0.62,0.72,55)  # 2:3
    p34_arr = np.linspace(0.7,0.8,55)    # 3:4
    p13_arr = np.linspace(0.28,0.38,55)  # 1:3
    p43_arr = np.linspace(1.28,1.38,55)  # 4:3
    p32_arr = np.linspace(1.45,1.55,55)  # 3:2
    p21_arr = np.linspace(1.95,2.05,55)  # 2:1
    p31_arr = np.linspace(2.95,3.05,55)  # 3:1
    # compile into big list
    resonance_arrays = [p12_arr,p23_arr,p34_arr,p13_arr,p43_arr,p32_arr,p21_arr,p31_arr]

    big_list = list(itertools.chain(p12_arr,p23_arr,p34_arr,p13_arr,p43_arr,p32_arr,p21_arr,p31_arr))
    Pratio_array = np.array(big_list)
    
    
    # initialize the log mass, e and w grids
    m1_array = np.linspace(np.log10(0.015), np.log10(10), 4) #5 earth mass as lower limit
    m2_array = np.linspace(np.log10(0.015), np.log10(10), 4) 
    w1_array = np.array([0.,60.,320.])
    w2_array = np.array([0.,60.,120.,180.,240.])
    e1_array = np.array([0.1,0.3,0.5])
    e2_array = np.array([0.1,0.3,0.5])
    Omega2_array = np.linspace(0,np.pi,3)
    M2_array = np.linspace(0, 2*np.pi, 4)
    
    
    os.system('rm /path/to/gridsearch/rebound_tested_params')

    theta_list = []
    with Pool(processes=60) as pool:
        for m1 in m1_array:
            for m2 in m2_array:
                for w1 in w1_array:
                    for w2 in w2_array:
                        for e1 in e1_array:
                            for e2 in e2_array:
                                for p1 in p1_array:
                                    for pratio in Pratio_array:
                                        for Omega2 in Omega2_array:
                                            for M2 in M2_array:
                                                p2 = p1*pratio
                                                sesinw1 = e1**0.5*np.sin(w1/180.*np.pi)
                                                secosw1 = e1**0.5*np.cos(w1/180.*np.pi)
                                                sesinw2 = e2**0.5*np.sin(w2/180.*np.pi)
                                                secosw2 = e2**0.5*np.cos(w2/180.*np.pi)
                                                theta = [m1, m2, p1, p2, sesinw1, secosw1, sesinw2, secosw2,Omega2,M2]
                                                theta_list.append(theta)
                                                

        pool.map(log_likelihood, theta_list)
