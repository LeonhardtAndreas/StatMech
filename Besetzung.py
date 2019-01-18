#!/usr/bin/env python3
#
# # Besetzungszahl-Statistik
# 
# Die Wahrscheinlichkeitsverteilung für die Besetzung eines Zustandes $i$ 
# bei gegebenem Erwartungswert $\overline{n_i}$.
import warnings
import numpy as np
from math import factorial
import matplotlib.pyplot as plt
import argparse

# ## Bosonische Besetzungswahrscheinlichkeit
# $$p_i^{BE} = \frac{\overline{n_i}^n}{(\overline{n_i}+1)^{n+1}}$$
def p_BE(n,n_i):
    return n_i**n/( (n_i+1)**(n+1) )


# ## Maxell-Boltzmann Besetzungswahrscheinlichkeit
# $$p_i^{MB}(n) = \frac{\overline{n_i}^n}{n!} \mathrm{e}^{-\overline{n_i}}$$
def p_MB(n,n_i):
        return n_i**n/factorial(n)*np.exp(-n_i)


# ## Fermi-Dirac Besetzungswahrscheinlichkeit
# $$p_i^{FD}(n) =
#                 \begin{array}{lr} 
#                     1-\overline{n_i} & n=0 \\ 
#                     \overline{n_i} & n =1 \\
#                     0 & n>1 \\
#                  \end{array} 
#              $$
def p_FD(n,n_i):
    if n_i>1:
        return None
    else:
        if n == 0:
            return 1-n_i
        elif n == 1:
            return n_i
        else:
            return None


# ## Interaktiver Plot
# $p^{FD}$ ist dabei natürlich nur für Besetzungszahlen $\le 1$ definiert.
def make_plot(n_i=0.5,N_max=20):
    be = np.zeros(N_max+1)
    mb = np.zeros(N_max+1)
    fd = np.zeros(N_max+1)
    for n in range(N_max+1):
        be[n] = p_BE(n,n_i)
        mb[n] = p_MB(n,n_i)
        fd[n] = p_FD(n,n_i)
    
    fig, ax = plt.subplots(1,2,figsize=(10,5))
    ax[0].plot(be,'bo-',label='BE')
    ax[0].plot(mb,'go-',label='MB')
    ax[0].plot(fd,'ro' ,label='FD')
    ax[0].set_xlim([0,N_max])
    ax[0].set_xlabel(r'$n$')
    ax[0].set_ylabel(r'$p_i(n)$')
    ax[0].legend(bbox_to_anchor=(.7,1),loc=2)
    
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ax[1].plot(np.log(be),'bo-',label='BE')
        ax[1].plot(np.log(mb),'go-',label='MB')
        ax[1].plot(np.log(fd),'ro' ,label='FD')
    ax[1].set_xlim([0,N_max])
    ax[1].set_xlabel(r'$n$')
    ax[1].set_ylabel(r'$\log p_i(n)$')
    plt.show()

if __name__ == '__main__':
    # argument handling
    parser = argparse.ArgumentParser(
                description='Plotte Besetzungswahrscheinlichkeiten'
                            +' bei gegebener mittlerer Besetzung.')
    parser.add_argument('n_i',type=float,default=0.5,
                            help='mittlere Besetzungszahl')
    parser.add_argument('N_max',type=int,nargs='?',default=20,
                            help='berechne p(n) bis zu diesem Wert')
    args = parser.parse_args()

    # call function
    make_plot(n_i=args.n_i,N_max=args.N_max)


# EOF
