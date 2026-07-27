import os,sys
from astropy.io import fits
import numpy as np
import pandas as pd
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

def stackspecplt(stackspec, velres, sampname, mh9, erand, lsm, sfr, nobj=0, pars=None):
    
    #   Plot stacked spectrum

    fig		= plt.figure(figsize=(3.0,2.4))
    ax 		= fig.add_axes([0.18, 0.16, 0.80, 0.82])
    ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)

    plt.axhline(0, c='c',ls='--',lw=0.5)
    plt.plot(stackspec[:,0], stackspec[:,1], 'bo-', markersize=3)	
    plt.plot(stackspec[:,0], stackspec[:,2],'r--',lw=0.5)
    plt.plot(stackspec[:,0], -stackspec[:,2],'r--',lw=0.5)

    plt.ylim([1.5*np.amin(stackspec[:,1]),1.2*np.amax(stackspec[:,1])])

    plt.figtext(0.21, 0.90, fr"log ($M_*$) = {lsm:.1f}", ha="left", fontsize=8)
    plt.figtext(0.21, 0.80, fr"SFR = {10.0**sfr:.2f}", ha="left", fontsize=8)

    plt.figtext(0.95, 0.90, f"N = {nobj}", ha="right", fontsize=8)
    plt.figtext(0.95, 0.80, fr"({mh9:.1f}$\pm${erand:.1f})$\times 10^9$", ha="right", fontsize=8)

    plt.xlabel('Velocity (km/s)')
    plt.ylabel('Luminosity density (Jy Mpc$^2$)')
    ax.yaxis.set_label_coords(-0.15, 0.5)
    plt.savefig(f"{pars['WorkDir']}/{pars['ResplotDir']}/{sampname}_{pars["StackStat"]}_spec_{velres}_kms.pdf", \
                                                        transparent=True, format='pdf')
    plt.close()

    return
#	--------------------------------------------------------------------------------------------------


def stackmaplt(stackmap, velres, sampname, pars=None):
    
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
    plt.savefig(f"{pars['WorkDir']}/{pars['ResplotDir']}/{sampname}_{pars["StackStat"]}_map_{velres}_kms.pdf", \
                                                        transparent=True, format='pdf')
    plt.close()

    return
#	--------------------------------------------------------------------------------------------------


def resultable(gsamp, pars=None):
    
    #   Create a data-frame from the results

    gsampd = pd.DataFrame([{'sampname': gsamp.sampname, 'ngal':gsamp.ngal, \
                            'meanlumint': gsamp.meanlumint, 'medlumint': gsamp.medlumint, \
                            'meanlumintavg': gsamp.meanlumintavg, 'medlumintavg': gsamp.medlumintavg, \
                            'meanmh9': gsamp.meanmh9, 'medmh9': gsamp.medmh9, \
                            'meanmh9avg': gsamp.meanmh9avg, 'medmh9avg': gsamp.medmh9avg, \
                            'emeanmhrand': gsamp.emeanmhrand, 'emedmhrand': gsamp.emedmhrand, \
                            'emeanpln': gsamp.emeanpln, 'emedpln': gsamp.emedpln, \
                            'emeanplnavg': gsamp.emeanplnavg, 'emedplnavg': gsamp.emedplnavg, \
                            'emeanmhjk': gsamp.emeanmhjk, 'emedmhjk': gsamp.emedmhjk, \
                            'emeanmhjkavg': gsamp.emeanmhjkavg, 'emedmhjkavg': gsamp.emedmhjkavg, \
                            'meanlb': gsamp.meanlb, 'medlb': gsamp.medlb, 'meanmb': gsamp.meanmb, 'medmb': gsamp.medmb, \
                            'meanz': gsamp.meanz, 'medz': gsamp.medz, \
                            'meanlsm': gsamp.meanlsm, 'medlsm': gsamp.medlsm, \
                            'meansfr': gsamp.meansfr, 'medsfr': gsamp.medsfr, \
                            'meanssfr': gsamp.meanssfr, 'medssfr': gsamp.medssfr, \
                            'meansfr10': gsamp.meansfr10, 'medsfr10': gsamp.medsfr10, \
                            'meansfr100': gsamp.meansfr100, 'medsfr100': gsamp.medsfr100, \
                            'meanssfr10': gsamp.meanssfr10, 'medssfr10': gsamp.medssfr10, \
                            'meanssfr100': gsamp.meanssfr100, 'medssfr100': gsamp.medssfr100}])    

    return (gsampd)
#	--------------------------------------------------------------------------------------------------

