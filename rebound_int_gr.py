import rebound
import reboundx
import numpy as np

import pandas as pd

#need reboundx 4.0.0
#need rebound 4.4.2

stellar_mass = 1.406714
stellar_radius = 1.667065

Mjup_Msun= 0.000954588
sim = rebound.Simulation()
sim.integrator = "ias15"

theta = [-0.0912494478,-0.7619896415,67.1231927214,33.9479027712,-0.3790572598,0.2305658055,0.5172276984,-0.1487319809,33.1883719080,225.2154336241,74.0433852429,1.4968345682,-4.8600033180,0.2982954315] 
logm1, logm2, p1, p2, sesinw1, secosw1, sesinw2, secosw2,Omega1,Omega2,inc2,M1,M2,b1 = theta
aAU = stellar_mass**(1/3)*(p1/365.25)**(2/3)	#a in AU
aors = aAU*215/stellar_radius

inc1 = 180/np.pi*np.arccos(b1/aors)
mp1 = 10**(logm1)*1.41/1.17
mp2 = 10**(logm2)*1.41/1.17
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
 
inc1 = (90-inc1)/180.*np.pi 
inc2 = (90-inc2)/180.*np.pi 





sim.add(m=stellar_mass)
print("rebound input mp1 mp2 a1 a2 e1 e2 inc1 inc2 w1 w2 Omega1 Omega2 M1 M2[angles in degrees for easy of comparison]", mp1, mp2, a1, a2, e1, e2, inc1/np.pi*180., inc2/np.pi*180., w1/np.pi*180., w2/np.pi*180., Omega1/np.pi*180., Omega2/np.pi*180., M1/np.pi*180., M2/np.pi*180.)
sim.add(m=mp1*Mjup_Msun, a=a1,e=e1, inc = inc1,Omega=Omega1,omega=w1,M=M1)
sim.add(m=mp2*Mjup_Msun, a=a2, e=e2, inc= inc2, Omega=Omega2, omega=w2,M=M2)
sim.move_to_com()
rebx = reboundx.Extras(sim)
gr = rebx.load_force("gr")
rebx.add_force(gr)
gr.params['c'] = 1.e4
dt = p1/365.25*2.*np.pi/50. 
#nstep = 3e9 # for the long integration
nstep = 3e5 # for the short integration
p = sim.particles
i = 0
while i<nstep:
    sim.integrate(sim.t+dt) # check for transits every dt time units. Note that 0.01 is shorter than one orbit
    i+=1
    if 1: # for the short integration
    #if i%1000==0: # for the long integration
        line = "%f" % sim.t 
        for o in sim.orbits():
            line+=" %f %f %f %f %f %f" % (o.a, o.e, o.inc, o.omega, o.Omega, o.M)
        print(line)
