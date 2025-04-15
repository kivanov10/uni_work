"""
Created by: Kristian Ivanov
URN: 6534278


This script will only work if there is a folder named LeoV_RT in the directory
with all .fits.gz files inside of it. The template file should be inside the
same directory as this script as seen bellow:
(Folder directory example)

spectra_assign.py    lte05000-3.00-1.5.fits   LeoV_RT <== 8 fits.gz files


Almost all plots have been commented out to save time on closing them manually,
after each iteration. They are all functional however and can be used by simply
uncommenting them.
"""


import os
import numpy as np
from scipy import ndimage as ndi
from statistics import mean, stdev
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits
from astropy.visualization import quantity_support
from astropy.modeling import models, fitting
from astropy.nddata import InverseVariance
from specutils.fitting import fit_generic_continuum, fit_continuum
from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import SplineInterpolatedResampler
from specutils.analysis import equivalent_width, snr_derived

def listdir_nohidden(path):
    """
    Input: path: str; Folder path

    Returns: lst: list of str; filenames following some naming rule

    Use to get the .fits.gz filenames for looping later on
    """
    lst = []
    for f in os.listdir(path):
        if not f.startswith('.') and f.endswith('.fits.gz'):
            lst.append(f)
    return lst


def file_extract(file_name):
    """
    Input: file_name: str; Name of the file that is going to be extracted from

    Return: file_data: class; data of the fits file
            file_head: class; header of the fits file

    Useful on its own, but especially so in loops. Third extension only and
    needs the LeoV_RT folder to exist and have its files in it.
    """
    file_open = fits.open('LeoV-RT/'+ file_name)
    print('file',file_open)
    file_data = file_open[2].data #Third extension
    file_head = file_open[2].header
    return file_data, file_head


def template():
    """
    Inputs: NONE

    Returns: templ_spec_norm : Spectrum1D; the normalised template spectra
             templ_head['CDELT1'] : float; the step in the x axis

    This function is entirely done to section the template creation in its own
    corner of the code. Done mainly for clarity.
    """
    templ = fits.open('lte05000-3.00-1.5.fits')
    templ_data = templ[0].data
    templ_head = templ[0].header
    templ_flux = templ_data

    # Creating lambda from header information
    awave = np.arange(templ_head['CRVAL1'],templ_head['CRVAL1']+\
                      templ_head['NAXIS1']*templ_head['CDELT1'],\
                      templ_head['CDELT1'])

    # Correcting to vacuum
    s = 10**4 / awave
    n = 1. - 0.00008336624212083 - \
    (0.02408926869968 / (130.1065924522 - s**2)) +\
    (0.0001599740894897 / (38.92568793293 - s**2))
    templ_lam = awave*n
    templ_spec = Spectrum1D(spectral_axis = templ_lam * u.AA, flux = templ_flux\
                       * u.Unit('erg cm-2 s-1 cm-1'))

    fit = fitting.LinearLSQFitter()
    line_init = models.Chebyshev1D(7)
    cont = fit(line_init, templ_lam, templ_flux)


    templ_spec_norm = Spectrum1D(spectral_axis = templ_lam * u.AA,\
                                 flux = templ_flux/cont(templ_lam) * u.dimensionless_unscaled)

    return templ_spec_norm, templ_head['CDELT1']




def vel_finder(spec, temp, vr):
    """
    Input: spec: Spectrum1D; noisy normalised spectra

           temp: Spectrum1D; normalised template spectra

           vr:   numpy array; velocity range to be tested with

    Returns: vr[best]: int; the velocity with the lowers chisq

    This function determines what is the best fit velocity from matching the
    template and file spectra
    """

    def resampler(spec, temp_norm, temp_new, delt, templ_delt):
        """
        Input: spec: Spectrum1D; the normalised noisy spectra

               temp_norm: Spectrum1D; the normalised template

               temp_new: numpy array; the shifted wavelength range by some velocity

               templ_delt: float; CDELT value of the template

        Returns: temp_resample; flux of the template resampled to match that of
                                the file spectra

        Exclusively used inside vel_finder. Resamples the spectra based on the
        CDELTs of the spec and template
        """
        fluxcon = SplineInterpolatedResampler()
        smooth = ndi.filters.gaussian_filter(temp_norm.flux,sigma=delt/temp_delt) * u.dimensionless_unscaled
        temp_norm_res = Spectrum1D(spectral_axis=temp_new,flux=smooth)
        temp_resample = fluxcon(temp_norm_res, spec.spectral_axis)
        # plt.plot(temp_resample.spectral_axis, temp_resample.flux)
        # plt.show()
        return temp_resample
    chisq = np.zeros(len(vr))
    for i in range(len(vr)):
        temp_new = temp.spectral_axis * (1 + vr[i]/3e5) # Moving it tothe expected velocity
        # Correct for resolution
        temp_checker = resampler(spec,temp_norm,temp_new,head['CD1_1'],temp_delt)
        chisq[i] = np.sum((spec.flux - temp_checker.flux) ** 2 *spec.uncertainty.array)
    best = np.abs(chisq-np.min(chisq)).argmin()
    return vr[best]


def metal(reg1,reg2,reg3,mag_v):
    """
    Simple metalicity function; returns the metalicity based on the equivalent
    width regions and the visible range magnitude supplied to it
    """
    regs = 0.5*reg1 + reg2 + 0.5*reg3
    a = -2.9
    b = 0.187
    c = 0.422
    d = -0.882
    e = 0.0133
    met = a + b * mag_v + c * regs + d * regs ** (-1.5) + e * mag_v* regs
    return met

def half_light_calc(vel_disp):
    """
    Input: vel_disp: float; the velocity dispersion

    Returns: ratio: float; Ratio of mass and light
             ratio_err: float; error of the mass light ratio

             m_half: float; half-light radius mass
             m_half_err: float; error of the mass
    """
    r_half = 70 * u.pc
    r_half_err = 27 * u.pc

    lumin     = 4.9 * 10**3 * u.solLum
    lumin_err = 2.2 * 10**3 * u.solLum

    mu        = 580 * u.Unit('pc-1 km-2 s2') *u.solMass

    m_half = mu*r_half*(vel_disp*u.Unit('km s-1'))**2
    m_half_err = m_half*(np.sqrt(r_half_err/r_half)**2)
    ratio = m_half/lumin
    ratio_err = ratio*(np.sqrt((m_half_err/m_half)**2+(lumin_err/lumin)**2))

    return ratio,ratio_err, m_half, m_half_err


#=============================Script starts here================================


quantity_support() #Helps with plots and axis units

file_lst = listdir_nohidden('LeoV-RT') #extracts the filenames of from the folder

file_lst.sort() #fix order of files

mag_vs = [-0.17449188,-1.7911892,-1.564951,0.34789276,-0.6503639,-1.231924,\
          -1.2690792,-0.7397213] #Mv of Leo stars in order

#Region limits
lower_lim = 8450
higher_lim  = 8700

#Create the template Spectrum1D and the cdelt
temp_norm, temp_delt = template()


#Master storage for velocities
vel_val = []
vel_std = []
vel_hist = []

v_range = np.arange(-500,500,1) #Test velocities

iterations = 10 #Iterations on Monte Carlo noise addition

for i,file_name in enumerate(file_lst):
    print('now doing',file_name)
    data, head = file_extract(file_name)

    flux_spec = data['SPEC'][0]
    flux_spec *= u.dimensionless_unscaled

    lambda_spec = data['LAMBDA'][0]

    invers = InverseVariance(data['IVAR'][0])

    #Bool array to slice the data based on limits
    region = (lambda_spec > lower_lim) &\
             (lambda_spec < higher_lim)

    spectral_class = Spectrum1D(spectral_axis=lambda_spec[region]*u.AA,\
                                flux=flux_spec[region],\
                                uncertainty=invers[region])

    #Plotting to check the unmodified fits data
    # plt.plot(spectral_class.wavelength, spectral_class.flux)
    # plt.ylabel('Flux')
    # plt.title(f'Unmodified data from file {file_name}')
    # plt.show()


    vel_noise = [] #Storage for noisy vel data
    fit = fitting.LinearLSQFitter()
    line_init = models.Chebyshev1D(10)
    for iter in range(iterations):
        sigma_ivar = 1/np.sqrt(invers.array[region])
        add_noise = np.random.normal(scale = sigma_ivar) * u.dimensionless_unscaled
        spec_mod = Spectrum1D(spectral_axis = spectral_class.wavelength,\
                              flux = spectral_class.flux+add_noise)

        #Plotting to compare noisy with non-noisy data
        # plt.plot(spec_mod.wavelength, spec_mod.flux, label='Noisy')
        # plt.plot(spectral_class.wavelength, spectral_class.flux,\
        #          alpha=0.5,label='Original')
        # plt.title('Comparison of before and after random noise addition')
        # plt.legend()
        # plt.show()


        cont = fit(line_init, spec_mod.spectral_axis.value, spec_mod.flux.value)
        flux_norm = spec_mod.flux/cont(spec_mod.spectral_axis.value) \
                    * u.dimensionless_unscaled


        spec_norm = Spectrum1D(spectral_axis=spec_mod.spectral_axis,\
                               flux=flux_norm * u.dimensionless_unscaled,\
                               uncertainty=spectral_class.uncertainty)

        # plt.plot(spec_norm.wavelength, spec_norm.flux,alpha=0.8,label='Normalised spectra')
        # plt.plot(temp_norm.wavelength, temp_norm.flux,label='Template data')
        # plt.xlim(lower_lim,higher_lim)
        # plt.ylabel('Flux')
        # plt.title('Normalised template and spectra data before velocity calculation')
        # plt.legend()
        # plt.show()

        vel_test = vel_finder(spec_norm,temp_norm,v_range)
        vel_noise.append(vel_test)
        vel_hist.append(vel_test)



    vel_val.append(mean(vel_noise)) #Mean velocity of the system
    vel_std.append(stdev(vel_noise)) #Velocity dispersion
    print('velocity',vel_val[i])
    print('velocity std',vel_std[i])

    new_lam = spec_norm.spectral_axis * (1 - vel_val[i] /(3e5)) #best lam estimate
    big_spec = Spectrum1D(spectral_axis=new_lam, flux=spec_norm.flux,\
                          uncertainty=spec_norm.uncertainty)

    # plt.plot(big_spec.wavelength, big_spec.flux, label='Velocity modified')
    # plt.plot(temp_norm.wavelength, temp_norm.flux,alpha=0.9,\
    #          label='Template data')
    # plt.xlim(lower_lim,higher_lim)
    # plt.ylabel('Flux')
    # plt.title('Normalised template and spectra data after velocity calculation')
    # plt.show()



    region_ca1 = SpectralRegion(8493 * u.AA, 8503 * u.AA)
    region_ca2 = SpectralRegion(8537 * u.AA, 8547 * u.AA)
    region_ca3 = SpectralRegion(8657 * u.AA, 8667 * u.AA)
    ew1 = equivalent_width(big_spec, regions=region_ca1)
    snr1 = snr_derived(big_spec, region=region_ca1)
    ew2 = equivalent_width(big_spec, regions=region_ca2)
    snr2 = snr_derived(big_spec, region=region_ca2)
    ew3 = equivalent_width(big_spec, regions=region_ca3)
    snr3 = snr_derived(big_spec, region=region_ca3)

    vmag = mag_vs[i]

    print('equivalent width:',ew1, ew2, ew3)
    print('Signal to noise:',snr1,snr2,snr3)
    print("[Fe/H] = ",metal(ew1.value,ew2.value,ew3.value,vmag))



print('Velocity mean = ',mean(vel_val))
print('Velocity dispersion = ',mean(vel_std))

# 2 * stds for 95% confidence interval
Ratio,Ratio_Err,M_Half,M_Half_err = half_light_calc(mean(vel_std)*2)

print('Ratio =',Ratio,"+/-",Ratio_Err)
print('M_half = ',M_Half,"+/-",M_Half_err)

plt.hist(vel_val,bins=5)
plt.title('Mean velocity for each Leo star file')
plt.xlabel('Velocity [km/s]')
plt.ylabel('Frequency')
plt.show()


plt.hist(vel_hist,bins=10)
plt.title('Velocity estimates from noisy flux ranges')
plt.xlabel('Velocity [km/s]')
plt.ylabel('Frequency')
plt.show()
