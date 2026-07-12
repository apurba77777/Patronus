import numpy as np
from scipy.integrate import quad

#	--------------------------------------------------------------------------------------------------
#	This file contains auxiliary functions useful for analysing H1 spectra
#	--------------------------------------------------------------------------------------------------


cc		=	299792458.0							#	c in m/s
hfac	=	0.70								#	Hubble parameter
parsec	=	3.08568e16							#	parsec/meter
omegal	=	0.70								#	Omega_Lambda
omegam	=	0.30								#	Omega_m

#-------------------------------------------------------------------------------------------------
#	Function to calculate the primary beam of GMRT in L band
def	lbpbeam(x):			#	x	=	distance from centre in arcminutes * freq in GHz
	a	=	-2.57888
	b	=	26.83336
	c	=	-12.76881
	d	=	2.29793
	y	=	1.0 + (a/10.0**3)*x**2 + (b/10.0**7)*x**4 + (c/10.0**10)*x**6 + (d/10.0**13)*x**8
	return y
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#	Function to calculate angular seperation between two celestial points in arc seconds
def distang(ra,ra0,dec,dec0):
	dra		=	(ra-ra0)*np.pi/180.0
	ddec	=	(dec-dec0)*np.pi/180.0
	dtheta	=	np.sqrt((dra*np.cos(dec0*np.pi/180.0))**2+ddec**2)
	dtheta	=	dtheta*180*60*60/np.pi
	return (dtheta)
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#	Function to calculate dr/dz as a function of z in 10^25 m
def drdz(z):															
	H0	=	100.0*hfac/3.08568e19
	ol	=	omegal
	om	=	omegam
	c	=	299792458.0
	x	=	(c/H0)/np.sqrt(om*((1+z)**3) + ol)
	return (x/1.0e25)
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#	Function to calculate luminosity distance (d_l) as a function of redshift in 10^25 m
def d_lum(z):
	dl	=	((1.0+z)*(quad(drdz,0.0,z))[0])
	return dl
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#	Function to calculate comoving distance (d_c) as a function of redshift in 10^25 m
def d_com(z):
	dc	=	(quad(drdz,0.0,z))[0]
	return dc
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#	Function to calculate luminosity distance (d_l) as a function of redshift in Gpc
def d_lum_gpc(z):
	dl	=	((1.0+z)*(quad(drdz,0.0,z))[0])
	return (dl*1.0e25/(parsec*1.0e9))
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#	Function to calculate comoving distance (d_c) as a function of redshift in Gpc
def d_com_gpc(z):
	dc	=	(quad(drdz,0.0,z))[0]
	return (dc*1.0e25/(parsec*1.0e9))
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#	Function to calculate relativistic velocity in km/s from frequency
def relvel(f0,f):
	rvel	=	cc * (f0*f0 - f*f)/(f0*f0 + f*f)	
	return (rvel/1.0e3)
#--------------------------------------------------------------------------------------------------

#	-----------------------------------------------------------------------------------------------
#	Function to calculate M_HI in M_sun from flux density in uJy 
def massflux(fluxujy,cwmhz,z):
	dlgpc	=	d_lum_gpc(z)
	mh		=	4.945e7*(dlgpc**2)*fluxujy*cwmhz
	return mh
#	-----------------------------------------------------------------------------------------------

#	-----------------------------------------------------------------------------------------------
#	Function to calculate M_HI in M_sun from luminosity density (flux density * dl^2) in uJy*Gpc^2 
def masslum(lumujygpc2,cwmhz):
	mh		=	4.945e7*lumujygpc2*cwmhz
	return mh
#	-----------------------------------------------------------------------------------------------

#	-----------------------------------------------------------------------------------------------
#	Function to calculate M_HI in M_sun from flux density in uJy and velocity in km/s
def massfluxvel(fluxujy,dvkmps,z):
	dlgpc	=	d_lum_gpc(z)
	mh		=	2.343e5*(dlgpc**2)*fluxujy*dvkmps/(1.0+z)
	return mh
#	-----------------------------------------------------------------------------------------------

#	-----------------------------------------------------------------------------------------------
#	Function to calculate M_HI in M_sun from luminosity density (flux density * dl^2) in uJy*Gpc^2 and velocity in km/s
def masslumvel(lumujygpc2,dvkmps,z):
	mh		=	2.343e5*lumujygpc2*dvkmps/(1.0+z)
	return mh
#	-----------------------------------------------------------------------------------------------

#	-----------------------------------------------------------------------------------------------
#	Function to calculate M_HI in M_sun from redshift weighted luminosity density (flux density * dl^2/(1+z)) in uJy*Gpc^2 and velocity in km/s
def masslumvel_noz(lumujygpc2,dvkmps):
	mh		=	2.343e5*lumujygpc2*dvkmps
	return mh
#	-----------------------------------------------------------------------------------------------

#	-----------------------------------------------------------------------------------------------
#	Function to estimate comoving volume in Gpc^3 between redshifts z0 and z1 and within a given solid angle
def	v_com_gpc3(z0,z1,solang):
	r0	=	d_com_gpc(z0)
	r1	=	d_com_gpc(z1)
	vol	=	(4.0/3.0)*np.pi*(r1**3-r0**3)
	vol	=	solang*vol/(4.0*np.pi)	
	return	vol
#	-----------------------------------------------------------------------------------------------

#	-----------------------------------------------------------------------------------------------
#	Function to return cosmological critical density at a given redshift in M_sun/mpc^3
def rho_crit(z):
	ol	=	omegal
	om	=	omegam
	gg	=	6.67408e-11
	msun=	2.0e30
	H0	=	100.0*hfac/3.08568e19
	rho	=	3*H0*H0*(ol+om*(1.0+z)**3)/(8*np.pi*gg)
	rho	=	rho*(1.0e6*parsec)**3/msun
	return rho
#	-----------------------------------------------------------------------------------------------

#	-----------------------------------------------------------------------------------------------
#	This function converts absolute B magnitude to log of B band luminosity in units of solar L_sun_B
def conmblbsun(mb):
	mbsun	=	5.48	# Binny and Merrifield 1998
	llb		=	(mbsun-mb)/2.5	
	return llb
#	-----------------------------------------------------------------------------------------------

#	-----------------------------------------------------------------------------------------------
#	Returns M_HI/M_sun for a given M_B using relation from Denes et al 2002 for delta < +2 deg
def mh1mb(mb):
	lmh		=	2.89 - 0.34*mb
	return (10.0**lmh)
#	-----------------------------------------------------------------------------------------------

#	-----------------------------------------------------------------------------------------------
#	Returns D_HI in kpc for a given M_HI using relation from Wang et al 2016
def	dh1mh(mh):
	ldh		=	0.506*np.log10(mh) - 3.293
	return (10.0**ldh)
#	-----------------------------------------------------------------------------------------------

#	-----------------------------------------------------------------------------------------------
#	Function to convert observed flux density to 1.4 GHz luminosity density in 10^22 W/Hz
	#	lnu		-	1.4 GHz luminosity density in 10^22 W/Hz
	#	sobs	-	Observed flux density in Jy
	#	z		-	redshift
	#	pl		-	power law index
	#	freq	-	observed frequecny in MHz
def s_lum(sobs,z,pl,freq):
	area	=	4*np.pi*(d_lum(z)**2)
	ffac	=	(freq*(1.0+z)/1420.0)**pl
	lnu		=	1.0e2*area*sobs*ffac/(1.0+z)	
	return lnu
#	-----------------------------------------------------------------------------------------------

#	-------------------------------------------------------------------------------------------------------------------------
#	Function to estimate SFR from stellar mass and redshift using the Main sequence parametrization from Whitaker et al. 2012
#	lsfr	-	log(SFR)
#	lsm		-	log(M*)
#	zz		-	redshift
def msfrw12(lsm, zz):
	az		=	0.70 - 0.13*zz
	bz		=	0.38 + 1.14*zz -0.19*zz*zz
	lsfr	=	az*(lsm-10.5) + bz
	#print np.median(lsm), np.median(zz), np.median(lsfr)
	return lsfr
#	-------------------------------------------------------------------------------------------------------------------------

#	------------------------------------------------------------------------------------------------
def	qbase(x, a, b, c):
	return (a*(x**2) + b*x +c)

def	cbase(x, a, b, c, d):
	return (a*(x**3) + b*(x**2) +c*x +d)
#	------------------------------------------------------------------------------------------------









