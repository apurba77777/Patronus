import os,sys
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
from specscripts.auxfns import *

plt.rc('legend', fontsize=8)    # legend fontsize
mpl.rcParams['font.size']=8
mpl.rcParams['lines.linewidth']=1
mpl.rcParams['axes.labelsize']=8

#   ---------------------------------------------------------------------------------------------------
#
#       Functions to plot stacked emission
#
#   ---------------------------------------------------------------------------------------------------

def stackspecplt(stackspec, velres, pars=None):
    
    #   Plot stacked spectrum

    fig		= plt.figure(figsize=(3.0,2.4))
    ax 		= fig.add_axes([0.18, 0.16, 0.80, 0.82])
    ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)

    plt.axhline(0, c='c',ls='--',lw=0.5)
    plt.plot(stackspec[:,0], stackspec[:,1], 'bo-', markersize=3)	
    plt.plot(stackspec[:,0], stackspec[:,2],'r--',lw=0.5)
    plt.plot(stackspec[:,0], -stackspec[:,2],'r--',lw=0.5)

    plt.xlabel('Velocity (km/s)')
    plt.ylabel('Luminosity density (Jy Mpc$^2$)')
    ax.yaxis.set_label_coords(-0.15, 0.5)
    plt.savefig(f"{pars['WorkDir']}/{pars['ResplotDir']}/{pars['StackName']}_{pars["StackStat"]}_spec_{velres}_kms.pdf", \
                                                        transparent=True, format='pdf')
    plt.close()

    return
#	--------------------------------------------------------------------------------------------------


def stackmaplt(stackmap, velres, pars=None):
    
    #   Plot stacked emission map

    fig		= plt.figure(figsize=(3.0,3.0))
    ax 		= fig.add_axes([0.25, 0.25, 0.7, 0.7])
    ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)

    plt.imshow(stackmap, interpolation='none', origin='lower', cmap=pars['ColMap'])
    plt.contour(stackmap,levels=np.array(pars["CLevels"])*np.nanmax(stackmap), colors='k', linestyles='-')
    plt.contour(stackmap,levels=-np.flip(np.array(pars["CLevels"]))*np.nanmax(stackmap), \
                colors='k', linestyles='--')
    plt.xticks(pars['BSize'] + np.array([-pars['PSize'],-pars['PSize']/2,0,pars['PSize']/2,pars['PSize']]), \
               np.array([-pars['PSize'],-pars['PSize']/2,0,pars['PSize']/2,pars['PSize']]))
    plt.yticks(pars['BSize'] + np.array([-pars['PSize'],-pars['PSize']/2,0,pars['PSize']/2,pars['PSize']]), \
               np.array([-pars['PSize'],-pars['PSize']/2,0,pars['PSize']/2,pars['PSize']]))

    plt.xlim([pars['BSize'] - pars['PSize'], pars['BSize'] + pars['PSize']])
    plt.ylim([pars['BSize'] - pars['PSize'], pars['BSize'] + pars['PSize']])

    plt.xlabel('RA offset (pixels)')
    plt.ylabel('Dec offset (pixels)')
    plt.savefig(f"{pars['WorkDir']}/{pars['ResplotDir']}/{pars['StackName']}_{pars["StackStat"]}_map_{velres}_kms.pdf", \
                                                        transparent=True, format='pdf')
    plt.close()

    return
#	--------------------------------------------------------------------------------------------------