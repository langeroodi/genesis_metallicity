import warnings
import numpy as np
from copy import deepcopy
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
from uncertainties import ufloat
from uncertainties import unumpy as unp

from ..data.lines import lines_dict
from .attenuation import KC13

warnings.filterwarnings('ignore', category=RuntimeWarning, message='divide by zero encountered in double_scalars')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in double_scalars')
warnings.filterwarnings('ignore', category=OptimizeWarning)

##############
# Line Class #
##############

class LINE:

    def __init__(self, line_flux, line_lambda):

        self.line_flux   = line_flux
        self.line_lambda = line_lambda

########################
# Emission Lines Class #
########################

class EMISSION_LINES:

    #---- initiating the class ----#

    # line_flux is array of ufloat
    def __init__(self, object, data_dict, ignore_Ha=False, print_progress=False):

        self.object = object

        if print_progress:
            print('--------------------------')
            print('object: %s' %object)
            print('--------------------------')

        try:
            self.Hd = LINE(data_dict['Hdelta'], lines_dict['Hdelta']['lambda'])
        except:
            self.Hd = LINE(ufloat(np.nan, np.nan), lines_dict['Hdelta']['lambda'])

        self.Hg = LINE(data_dict['Hgamma'], lines_dict['Hgamma']['lambda'])
        self.Hb = LINE(data_dict['Hbeta'], lines_dict['Hbeta']['lambda'])
        self.Ha = LINE(data_dict['Halpha'], lines_dict['Halpha']['lambda'])

        if ignore_Ha:
            self.Ha = LINE(ufloat(np.nan, np.nan), lines_dict['Halpha']['lambda'])

        #---- setting up the flux and fluxerr arrays ----#

        # normalizing the Balmer lines to Hb
        try:
            norm_Hd_flux = self.Hd.line_flux/self.Hb.line_flux
        except:
            norm_Hd_flux = ufloat(np.nan, np.nan)

        try:
            norm_Hg_flux = self.Hg.line_flux/self.Hb.line_flux
        except:
            norm_Hg_flux = ufloat(np.nan, np.nan)

        try:
            norm_Hb_flux = self.Hb.line_flux/self.Hb.line_flux
        except:
            norm_Hb_flux = ufloat(np.nan, np.nan)

        try:
            norm_Ha_flux = self.Ha.line_flux/self.Hb.line_flux
        except:
            norm_Ha_flux = ufloat(np.nan, np.nan)

        # making the arrays
        flux_array    = np.array([norm_Hd_flux.n, norm_Hg_flux.n, norm_Hb_flux.n, norm_Ha_flux.n])

        fluxerr_array = np.array([norm_Hd_flux.s,
                                  norm_Hg_flux.s,
                                  norm_Hb_flux.s,
                                  norm_Ha_flux.s])

        mask          = np.where(~((flux_array == 0) | (np.isnan(flux_array))))[0]

        flux_array    = flux_array[mask]
        fluxerr_array = fluxerr_array[mask]

        # flooring the errors
        for i in range(len(flux_array)):
            if fluxerr_array[i] == 0:
                fluxerr_array[i] = flux_array[i] * 0.05

        #---- function for reddening ----#

        def redden(line_flux, line_lambda, Av):

            lambda_array = np.linspace(line_lambda, 1e+4, 1000)
            ext = KC13(lambda_array, Av, delta=0, Eb=0, return_AxAv=True)
            ext = 10**(-0.4*ext*Av)
            line_flux = line_flux * ext[0]

            return line_flux

        #---- function for dereddening ----#

        def deredden(line_flux, line_lambda, Av):

            lambda_array = np.linspace(line_lambda, 1e+4, 1000)
            ext = KC13(lambda_array, Av, delta=0, Eb=0, return_AxAv=True)
            ext = 10**(-0.4*ext*Av)
            line_flux = line_flux / ext[0]

            return line_flux

        #---- function that models the observed lines ----#

        def model_lines(x, Av):

            # built-in EMISSION_LINESs
            HbHd = 3.86
            HbHg = 2.14
            HaHb = 2.86

            # line fluxes
            Hd_flux = 1.00/HbHd
            Hg_flux = 1.00/HbHg
            Hb_flux = 1.00
            Ha_flux = 1.00*HaHb

            # reddening
            Hd_flux = redden(Hd_flux, self.Hd.line_lambda, Av)
            Hg_flux = redden(Hg_flux, self.Hg.line_lambda, Av)
            Hb_flux = redden(Hb_flux, self.Hb.line_lambda, Av)
            Ha_flux = redden(Ha_flux, self.Ha.line_lambda, Av)

            # output
            if Av < 0:
                output_array = np.array([np.inf, np.inf, np.inf, np.inf])
            else:
                output_array = np.array([Hd_flux, Hg_flux, Hb_flux, Ha_flux])/Hb_flux
            output_array = output_array[mask]

            if False:
                try:
                    print(Av, ';', output_array[1]/output_array[0], output_array[2]/output_array[1])
                except:
                    try:
                        print(Av, ';', output_array[1]/output_array[0])
                    except:
                        print(Av, ';', output_array[2]/output_array[1])

            return output_array

        #---- fitting ----#

        if len(mask) > 1:

            if print_progress:
                print('flux', flux_array)
                print('fluxerr', fluxerr_array)

            popt, pcov = curve_fit(model_lines, np.array([0,0,0]), flux_array, sigma=fluxerr_array, absolute_sigma=True)

            Av       = popt[0]
            Av_sigma = np.sqrt(np.diag(pcov))

            self.Av  = Av

        if len(mask) <= 1:
            self.Av = np.nan

        try:
            if print_progress: print('HbHd (3.86) = ', self.Hb.line_flux/self.Hd.line_flux)
        except: a = 1

        try:
            if print_progress: print('HbHg (2.14) = ', self.Hb.line_flux/self.Hg.line_flux)
        except: a = 1

        try:
            if print_progress: print('HaHb (2.86) = ', self.Ha.line_flux/self.Hb.line_flux)
        except: a = 1

        if print_progress: print('Av          = ', self.Av)

        #---- dereddening all the lines ----#

        self.corrected_dict = {}

        for line in data_dict.keys():
            if line not in ['redshift', 'metallicity', 'red._corr.']:

                try:
                    line_flux   = deepcopy(data_dict[line])
                    line_lambda = lines_dict[line]['lambda']

                    if (not data_dict['red._corr.']) and (self.Av > 0.01):
                        line_flux   = deredden(line_flux, line_lambda, self.Av)

                    self.corrected_dict[line] = line_flux

                except:
                    self.corrected_dict[line] = np.nan

########
# Test #
########

if False:

    for object in data_dict.keys():
        test = EMISSION_LINES(object, data_dict[object])
