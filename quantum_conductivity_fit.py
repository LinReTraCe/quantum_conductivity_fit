#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys

import numpy as np
import matplotlib.pyplot as plt
from   mpmath import zeta
import scipy.optimize

'''
quantum_conductivitiy_fit.py: user script to fit experimental data
  with the 2L conductivitiy equation provided in Eq. (3) in arXiv:2008.05846


In main:
Provide a text file for the experimental data:
  1st column : temperatures [K]
  2nd column : resistivities [arbitrary units]

Estimate Scattering rate [eV], Gap [eV] and Prefactor for
brute force fitting: ranges argument in scipy.optimize.brute
function call.

When publishing results obtained with this script, please cite
  Resistivity Saturation in Kondo Insulators
  M. Pickem, E. Maggio, J. M. Tomczak
  arXiv:2008.05846 (2020)
'''

__author__     = 'Matthias Pickem'
__credits__    = ['Matthias Pickem', 'Emanuele Maggio', 'Jan M. Tomczak']

__license__    = 'MIT'
__version__    = '1.0'
__maintainer__ = 'Matthias Pickem'
__email__      = 'matthias.pickem@gmail.com'
__status__     = 'Production'


def sigma_fit(params, *args):
  ''' Fit function for scipy optimize to fit theortical equation to experimental data.

  Input:
    params: Fitting parameters [log10(Prefactor), log10(scattering rate), gap]
            in form of one variable
            Units:
              Prefactor: arbitrary units
              Scattering rate: eV
              Gap: eV
    args:   experimental temperature and conductivity data
            Units:
              temperature: K
              Conductivity: arbitrary units

  Output: Sum of squares between theoretical equation and experimental data
          on a log scale.
  '''


  ''' Fitting parameters '''
  PF, Gamma, Delta = params  # prefactor, scattering rate, gap
  PF    = 10**PF
  Gamma = 10**Gamma

  ''' Experimental data range '''
  temps = args[0]
  cond  = args[1]

  kB = 8.617333e-5  # Boltzmann constant [eV/K]
  values = []

  for iT in temps:
    beta = 1./(kB*iT)

    '''
    Polygamma function is connected to the generalized Zeta function

    Hurwitz zeta function:
    Zeta(s,z) = \sum_k=0^\inf 1/(z + k)^s

    Polygamma function:
    Psi(m,z) = (-1)^(m+1) m! \sum_k=0^\inf 1/(z+k)^m+1

    -> Psi(1,z) = \sum_k 1/(z+k)^2 = Zeta(2,z)
    -> Psi(2,z) = (-1) * 2! \sum_k 1/(z+k)^3 = Zeta(3,z) * (-2)
    '''

    psi_1 = zeta(2,0.5+beta/2/np.pi*(Gamma+1j*Delta/2.))
    psi_2 = zeta(3,0.5+beta/2/np.pi*(Gamma+1j*Delta/2.)) * (-2)

    ''' Eq. (3) '''
    value = PF * 2 * beta/(4 * np.pi**3 * Gamma) * (psi_1.real - beta*Gamma/(2*np.pi) * psi_2.real)

    values.append(float(value))
  values = np.array(values)

  ''' Minimize the sum of the squares on the log scale '''
  return np.sum((np.log10(values) - np.log10(cond))**2)




def sigma_save(temp, params, plot=True, save=None):
  ''' Plot function for theoretical equation with provided parameters.

  Same function as above with variable temperature range
  and options to plot the result and/or save the result
  to a text file.

  Input:
    temp:   temperature range [K]
    params: [log10(Prefactor), log10(scattering rate), gap]
            in form of one variable
  '''

  PF, Gamma, Delta  = params
  PF    = 10**PF
  Gamma = 10**Gamma

  kB = 8.617333e-5  # Boltzmann constant [eV/K]

  values = []
  for iT in temp:
    beta = 1./(kB*iT)

    psi_1 = zeta(2,0.5+beta/2/np.pi*(Gamma+1j*Delta/2.))
    psi_2 = zeta(3,0.5+beta/2/np.pi*(Gamma+1j*Delta/2.)) * (-2)

    value = PF * 2 * beta/(4 * np.pi**3 * Gamma) * (psi_1.real - beta*Gamma/(2*np.pi) * psi_2.real)
    values.append(value)
  values = np.array(values)

  if plot:
    plt.plot(temp, values, linestyle=':', linewidth=3, label='quantum conductivity Eq. 3')

  if save is not None:
    np.savetxt(str(save), np.transpose([temp, values, 1./values]), \
               header='T[K] | Conductivity | Resistivity')


def main():
  '''
  Read in experimental data.
  Optimize theoretical parameters.
  Plot and save theoretical fit.
  '''

  ''' Read in experimental data '''
  data = np.genfromtxt('resistivity.dat', usecols=(0,1))  # should contain columns: T[K] rho(T) [a.u.]
  temp = data[:,0]
  cond = 1./data[:,1] # convert resistivity --> conductivity

  ''' Sort into ascending temperature order '''
  xargs = np.argsort(temp) # get sort indices foget sort indices forng order
  temp  = temp[xargs] # sort
  cond  = cond[xargs] # sort

  ''' Restrict fitting interval to low temperatures '''
  restrict = True
  if restrict:
    t1 = np.searchsorted(temp, 0.)
    t2 = np.searchsorted(temp, 15.) # define 2 points
    tempfit = temp[t1:t2+1]
    condfit = cond[t1:t2+1]
  else:
    tempfit = temp
    condfit = cond

  ''' Define fine output temperature interval '''
  temp_output = np.linspace(0.001,100,2000) # temperature grid for output of fit to ascii

  ''' Optimization:
  ranges: first argument:  log10[prefactor] arb.units
          second argument: log10[scattering rate] eV
          third argumetn:  gap eV
  '''

  print("optimization in progress")
  resbrute = scipy.optimize.brute(func = sigma_fit, \
                                  ranges = (slice(-5.5,-4.5,0.1), slice(-3.,-2,0.1),slice(0.004,0.010,0.002)), \
                                  args = (tempfit, condfit), \
                                  full_output=True)

  # prefactor can be arbitrarily small or large
  #  ... its best to first fit by hand and then run brute force
  #  prefactor can also be roughly estimated via the T->0 expression

  # scattering rate between 1meV and 10 meV
  # Gap between 4 meV and 10 meV

  print("Results:")
  print("scattering rate [meV], gap [meV]")
  print("  ",10**(resbrute[0][1]) * 1e3,",", resbrute[0][2] * 1e3) # <- parameters
  print("Fit deviation")
  print("  ",resbrute[1])                                          # <- quality measure of fit
  fitparams = [resbrute[0][0],resbrute[0][1],resbrute[0][2]]       # use the brute force result


  ''' display experimental conductivity '''
  plt.plot(temp, cond, linewidth=4, label='Experiment')
  ''' plot and export theoretical fit '''
  sigma_save(temp_output, fitparams, plot=True, save="resistivity_fit.dat")

  plt.yscale('log')
  plt.xlim(0,20)
  plt.ylim(1e-2,1e-1)
  plt.legend()
  plt.xlabel('T [K]')
  plt.ylabel('conductivity')
  plt.tight_layout()
  plt.show()


if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    sys.exit('Exit by KeyboardInterrupt.')
  except BaseException as s:
    print('Exception occured: {}'.format(str(s)))
    raise
