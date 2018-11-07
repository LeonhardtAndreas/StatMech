#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

##################################################################
# Define constants
k = 1
m = 1
d = 3
N = 10
epsilon = 0.01

delta_t = 1e-1
N_t = int(1e6)

##################################################################
# calculate force for a give configuration #######################
def f(q):
    f = np.zeros_like(q)
    # interaction to the left
    f[:-1] += np.sign(q[1:]-q[:-1]+d)*epsilon/(q[1:]-q[:-1]+d)**2 
    # is minus the interaction to the right
    f[1:] -=f[:-1]
    # force from the potential
    f -= k*q
    return f
                
##################################################################
# Verlet integration #############################################
def time_step(q,q_old):
    q_new = -q_old + 2*q + delta_t**2*f(q)/m
    return [q_new, q]



##################################################################
# run as program #################################################
if __name__ == '__main__':

    # define arrays to be used later and set initial condition
    q_old = np.zeros(N)
    # move the first oscillator out of equilibrium 
    q_old[0] = -np.sqrt(2)
    # total energy is then given by
    E0 = k
    # perform first time step
    q = q_old + 1/2*delta_t**2*f(q_old)/m

    E_kin = np.zeros(N_t)
    E_pot = np.zeros(N_t)
    # further time steps
    for t in range(1,N_t):
        # calculate kinetic energy
        E_kin[t] = ((q-q_old)**2).sum()
        # calculate the potential energy of the left oscillator
        E_pot[t] = q[0]
        # overwrite the state vectors
        q, q_old = time_step(q,q_old)
    E_kin *= m/2/delta_t**2/E0
    E_pot = k/2*E_pot**2

    fig, [ax1,ax2] = plt.subplots(2,1,figsize=[5,8])
    fig.suptitle(r'$N=$'+str(N)+r'$, \varepsilon=$'+str(epsilon))

    ax1.plot(np.arange(N_t)*delta_t,E_kin)
    ax1.set_xlim([0,N_t*delta_t])
    ax1.set_xlabel(r'$t$')
    ax1.set_ylabel(r'$E_{kin}/E_0$')
    
    bins = np.arange(21)*.01
    ax2.hist(E_pot,bins)
    ax2.set_xlim([bins[0],bins[-1]])
    ax2.set_xlabel(r'$E_{pot}$')
    ax2.set_ylabel(r'$p(E_{pot}$')

    plt.show()
