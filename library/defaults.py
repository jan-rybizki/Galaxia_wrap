#!/usr/local/bin/python
'''
@author: Tri L. Astraatmadja, MPIA Heidelberg
@contact: astraatmadja@mpia.de
'''

import math
import matplotlib
matplotlib.use('Agg') ## use a non-interactive Agg background
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoLocator, MultipleLocator, AutoMinorLocator, LinearLocator, LogLocator

import numpy as np

def plot_values(fontsize):
   plt.rcParams['image.cmap'] = 'viridis'
   plt.rc('text', usetex=True)
   plt.rc('font', family='serif')
   plt.rc('font', size=fontsize)
   plt.rc('axes', titlesize='large')
   plt.rc('axes', labelsize='medium')
   plt.rc('xtick', labelsize='x-small')
   plt.rc('ytick', labelsize='x-small')
   plt.rc('legend', fontsize='small')

# Default plot values
def def_plot_values():
   plot_values(12)

def def_plot_values_large():
   plot_values(16)

def def_plot_values_extra_large():
   plot_values(20)

def getLogTickMarks(xmin, xmax, dx=1):
   minTicks = math.floor(math.log10(xmin))
   maxTicks = math.ceil(math.log10(xmax))

   mantissa = np.arange(0,10,dx)
   mantissa[0] = 1

   tickMarks = []
   for i in np.arange(minTicks, maxTicks+1, 1):
      for j in mantissa:
         thisTick = j*math.pow(10., i)
         if ((thisTick > xmin) & (thisTick < xmax)):
            tickMarks.append(thisTick)
   return np.asarray(tickMarks)

def drawCommonLabel(xlabel, ylabel, fig, xPad=15, yPad=20):
   ax = fig.add_subplot(111)    # The big subplot

   # Turn off axis lines and ticks of the big subplot
   ax.spines['top'].set_color('none')
   ax.spines['bottom'].set_color('none')
   ax.spines['left'].set_color('none')
   ax.spines['right'].set_color('none')
   ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
   ax.patch.set_visible(False)
                              
   ax.get_xaxis().set_ticks([])
   ax.get_yaxis().set_ticks([])

   ax.set_xlabel(xlabel, labelpad=xPad)
   ax.set_ylabel(ylabel, labelpad=yPad)
   return ax

def APText(index):
   APText = ['TEFF', 'A0', 'LOGG', 'FEH']
   return APText[index]

# Write AP label in TeX format, plus an additional something as a subscript
def APLabel(index, text, Units=True):
   APLabel1 = [r'$T_{\rm eff', \
               r'$A_{0', \
               r'$\log g_{\rm',\
               r'$[{\rm Fe}/{\rm H}]_{\rm']
   if Units:
      APLabel2 = [r'}\;[{\rm K}]$', \
                  r'}\;[{\rm mag}]$', \
                  r'}$', \
                  r'}$']
   else:
      APLabel2 = [r'}$', \
                  r'}$', \
                  r'}$', \
                  r'}$']

   return APLabel1[index] + text + APLabel2[index]
