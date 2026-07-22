import os,sys
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
from specscripts.auxfns import *

plt.rc('legend', fontsize=10)    # legend fontsize
mpl.rcParams['font.size']=10
mpl.rcParams['lines.linewidth']=1
mpl.rcParams['axes.labelsize']=12

#   ---------------------------------------------------------------------------------------------------
#
#       Functions to plot sample statistics 
#
#   ---------------------------------------------------------------------------------------------------

def noiseratplt(setrats, pars=None):
    
    #   Plot noise ratios

    fig		= plt.figure(figsize=(3.2,2.8))
    ax 		= fig.add_axes([0.15, 0.15, 0.82, 0.84])
    ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
    ax.hist(setrats, bins=100)
    ax.axvline(x=1.0/np.sqrt(pars['CavgFac']), c='k',ls='--')
    ax.axvline(x=pars['RatLims'][0], c='r',ls='--',lw=0.5)
    ax.axvline(x=pars['RatLims'][1], c='r',ls='--',lw=0.5)

    ax.set_ylabel(r'Number of galaxies', fontsize=10)
    ax.yaxis.set_label_coords(-0.11, 0.5)

    ax.set_xlabel(r'Noise (%d km s$^{-1}$) / Noise (%d km s$^{-1}$)'%(int(10*pars['VelRes']), int(pars['VelRes'])), \
                  fontsize=10)

    plt.savefig(pars['WorkDir']+pars['ResplotDir']+pars['SamPlotDir']+"/noiserat_"+pars['StackName']+"_"+pars['ColType']+".pdf", 
                transparent=True, format='pdf')
    plt.close()

    return
#	--------------------------------------------------------------------------------------------------


def zlsmplt(detcat, pars=None):
    
    #   Plot redshift stellar mass distribution

    carr    = detcat[:,10]
    carr[carr < -3.0] = -3.0

    fig		= plt.figure(figsize=(3.2,2.8))
    ax 		= fig.add_axes([0.15, 0.15, 0.82, 0.84])
    ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
    
    plt.scatter(detcat[:,4], detcat[:,11], c=carr, s=3, marker="o", cmap="managua", alpha=0.6)

    ax.set_ylabel(r'log ($M_*$ / $M_{\odot}$)', fontsize=10)
    ax.yaxis.set_label_coords(-0.11, 0.5)
    ax.set_xlabel(r'Redshift (z)', fontsize=10)

    plt.savefig(pars['WorkDir']+pars['ResplotDir']+pars['SamPlotDir']+"/zlsm_"+pars['StackName']+"_"+pars['ColType']+".pdf", 
                transparent=True, format='pdf')
    plt.close()

    return
#	--------------------------------------------------------------------------------------------------
