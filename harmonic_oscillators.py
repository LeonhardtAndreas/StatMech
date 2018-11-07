#!/usr/bin/env python3
# # Coupled harmonic oscillators
# 
# At first we get some modules, most importantly numpy for calculations and pyplot for plotting
import numpy as np
from math import factorial, gamma
import matplotlib.pyplot as plt


# # Defining the System
# We define the constants in the Hamiltonian as given in the
# exercise.

# spring constant
k = 1
# mass of one oscillator
m = 1
# distance between two oscillators
d = 3


# # Functions for the Time Evolution
# ## The force for a given configuration
# f takes a state and returns the force acting on it, i.e. the gradient of 
# the potential and interaction energy. 
# $$ f(q) = -\nabla_q ( H_{pot}(q) + H_{int}(q) )$$
# 
# 1. First we add the interaction force from the left neighbour by calculating the distances. 
#     Note the indexing: 
#     - q[1:] runs from the second to the last element, 
#     - q[-1] denotes the last element and 
#     - q[:-1] = q[0:-1] runs from the first to the second last element.
# 2. The second step adds interactions from the right neighbour, which is just minus 
#     the one from that site to the right neighbour as calculated before.
#     Therefore we only shift what we already have and subtract it.
# 3. The last step adds just the force of the harmonic oscillator.

def f(q,epsilon):
    f = np.zeros_like(q)
    # interaction to the left
    f[:-1] += np.sign(q[1:]-q[:-1]+d)*epsilon/(q[1:]-q[:-1]+d)**2 
    # is minus the interaction to the right
    f[1:] -=f[:-1]
    # force from the potential
    f -= k*q
    return f
        


# ## Verlet integration
# Verlet integration is a reformulation of $\frac{d^2q}{dt^2} = \frac{1}{m}f(q)$
# for a discretized second derivative $$\frac{d^2 q}{dt^2} = \frac{ q(t+dt) + q(q-dt) - 2q(t) }{dt^2}$$
# 
# Pay attention of the usage later, when q_new is written to q, and q to q_old.

def time_step(q,q_old,delta_t,epsilon):
    q_new = -q_old + 2*q + delta_t**2*f(q,epsilon)/m
    return [q_new, q]


# # Analyical calculation for interaction free case
# Define the analytical solution from a) for comparison

def E_pot_analytical(E,E0,N):
    return E0**(1-N)*(E0-E)**(N-3/2)/np.sqrt(E*np.pi)*gamma(N)/gamma(N-1/2)


# # Calculating and Plotting
# Now we have everthing ready to put it together. 
# We set initial conditons (move the 0th oscillator out of equilibrium), calculate the first time step manually and 
# start the iteration. For each time step we calculate and save the total kinetic Energy E_kin and the potential energy of the zeroth oscillator.
# 
# The rest is plotting, labeling axes etc.

def make_plots(N=10,epsilon=0.01,delta_t=0.1,N_t=1e4):
    N_t = int(N_t)
    # define arrays to be used later and set initial condition
    q_old = np.zeros(N)
    # move the first oscillator out of equilibrium 
    q_old[0] = -np.sqrt(2)
    # total energy is then given by
    E0 = k
    # perform first time step
    q = q_old + 1/2*delta_t**2*f(q_old,epsilon)/m

    E_kin = np.zeros(N_t)
    E_pot = np.zeros(N_t)
    # further time steps
    for t in range(1,N_t):
        # calculate kinetic energy
        E_kin[t] = ((q-q_old)**2).sum()
        # calculate the potential energy of the left oscillator
        E_pot[t] = q[0]
        # overwrite the state vectors
        q, q_old = time_step(q,q_old,delta_t,epsilon)
    E_kin *= m/2/delta_t**2/E0
    E_pot = k/2*E_pot**2

    fig, [ax1,ax2] = plt.subplots(2,1,figsize=[5,8])
    fig.suptitle(r'$N=$'+str(N)+r'$, \varepsilon=$'+str(epsilon))

    ax1.plot(np.arange(N_t)*delta_t,E_kin)
    ax1.set_xlim([0,N_t*delta_t])
    ax1.set_xlabel(r'$t$')
    ax1.set_ylabel(r'$E_{kin}/E_0$')
    
    bins = np.arange(21)*.01
    ax2.hist(E_pot,bins,density=True)
    ax2.set_xlim([bins[0],bins[-1]])
    
    ax2.plot(np.arange(3,100)*bins[-1]/100,E_pot_analytical(np.arange(3,100)*bins[-1]/100,E0,N),'r')
    
    ax2.set_xlabel(r'$E_{pot}$')
    ax2.set_ylabel(r'$p(E_{pot})$')

    plt.show()


# # Calling the function
# 
# This lines calls the function defined above, with sliders for N and $\epsilon$ 
# as well as the Number of time iterations and the time step dt. 
# The other parameters need to be changed in the code, if needed. Remember to re-evaluate the cells after making changes (Shift+Enter).
# 
# *Time step and duration:*
# The python code is not very efficient, a calculation of $10^7$ time steps takes very long.
# $10^6$ with a time step of 0.1 instead of 0.01 does a good job as well in a reasonable time (about a minute).


if __name__ == '__main__':
    make_plots(N=10,epsilon=0.1,delta_t=1e-2,N_t=1e4)

