import numpy as np
from copy import deepcopy
from uncertainties import ufloat
from uncertainties import unumpy as unp

from .data.lines import lines_dict, backend_lines, print_lines
from .dust.extinction_correction import EMISSION_LINES
from .metallicity.direct_method import METALLICITY
from .metallicity.strong_method import measure_metallicity

#######################
# genesis-metallicity #
#######################

class genesis_metallicity:

    def __init__(self, input_dict, object='default', correct_extinction=True):

        #----------------------------------------#
        #---- verifying the input dictionary ----#
        #----------------------------------------#

        #---- converting the input_dict to data_dict ----#

        data_dict = {}

        # object
        self.object = object

        # reading the redshift if it is provided
        try:
            data_dict['redshift'] = input_dict['redshift']
        except:
            data_dict['redshift'] = np.nan

        # deciding if extinction correction has to be done
        data_dict['red._corr.'] = True
        if correct_extinction:
            data_dict['red._corr.'] = False

        # reading the line fluxes
        for line in lines_dict.keys():
            if line in input_dict.keys():
                line_flux = input_dict[line]
                data_dict[line] = ufloat(line_flux[0], line_flux[1])

        #---- OII ----#

        if 'OII' not in data_dict.keys():
            if ('O3727' not in data_dict.keys()) or ('O3729' not in data_dict.keys()):
                print_lines()
                raise ImportError('[OII]3727,29 flux is required! please provide it under the \'OII\' key (or alternatively under the \'O3727\' and \'O3729\' keys) in the input dictionary')
            else:
                data_dict['OII'] = data_dict['O3727'] + data_dict['O3729']

        #---- O4959 and O5007 ----#

        if ('O4959' not in data_dict.keys()) or ('O5007' not in data_dict.keys()):
            if 'OIII' not in data_dict.keys():
                print_lines()
                raise ImportError('[OIII]4959,5007 flux is required! please provide it under the \'OIII\' key (or alternatively under the \'O4959\' and \'O5007\' keys) in the input dictionary')
            else:
                data_dict['O4959'] = data_dict['OIII']/(1+2.98)
                data_dict['O5007'] = data_dict['OIII']/(1+2.98)*2.98

        #---- Hbeta ----#

        if 'Hbeta' not in data_dict.keys():
            print_lines()
            raise ImportError('Hbeta flux is required! please provide it under the \'Hbeta\' key in the input dictionary')

        #---- EWHb ----#

        if 'Hbeta_EW' not in data_dict.keys():
            print_lines()
            raise ImportError('Hbeta equivalent width is required! please provide it under the \'Hbeta_EW\' key in the input dictionary')

        #---- O7320 and O7330 ----#

        if ('O7320' in data_dict.keys()) and ('O7330' in data_dict.keys()):
            data_dict['OII7320'] = data_dict['O7320'] + data_dict['O7330']

        #---- making sure all the backend lines have some input ----#

        for line in backend_lines:
            if line not in data_dict.keys():
                data_dict[line] = ufloat(np.nan, np.nan, np.nan)

        #-------------------------------#
        #---- extinction correction ----#
        #-------------------------------#

        emission_lines                 = EMISSION_LINES(object, data_dict)
        self.Av                        = emission_lines.Av
        self.reddening_corrected_lines = emission_lines.corrected_dict

        #---------------------#
        #---- metallicity ----#
        #---------------------#

        #---- deciding on the metallicity measurement approach ----#

        metallicity_method = 'strong'

        try:
            if (data_dict['O4363'].n/data_dict['O4363'].s) > 1.0:
                metallicity_method = 'direct'
        except:
            metallicity_method = 'strong'

        self.metallicity_method = metallicity_method

        #---- direct-method metallicity ----#

        if self.metallicity_method == 'direct':

            direct_metallicity = METALLICITY(object, self.reddening_corrected_lines)
            self.metallicity   = direct_metallicity.metallicity
            self.t2            = direct_metallicity.Te_OII
            self.t3            = direct_metallicity.Te_OIII

        #---- strong-line metallicity ----#

        # O3
        def calculate_O2():
            line_ratio = self.reddening_corrected_lines['OII']/self.reddening_corrected_lines['Hbeta']
            return line_ratio

        # O3
        def calculate_O3():
            line_ratio = self.reddening_corrected_lines['O5007']/self.reddening_corrected_lines['Hbeta']
            return line_ratio

        # O32
        def calculate_O32():
            line_ratio = self.reddening_corrected_lines['O5007']/self.reddening_corrected_lines['OII']
            return line_ratio

        # R23
        def calculate_R23():
            line_ratio = (self.reddening_corrected_lines['OII']+self.reddening_corrected_lines['O4959']+self.reddening_corrected_lines['O5007'])/self.reddening_corrected_lines['Hbeta']
            return line_ratio

        if self.metallicity_method == 'strong':

            log_O2       = unp.log10([calculate_O2()])[0]
            log_O3       = unp.log10([calculate_O3()])[0]
            log_Hbeta_EW = unp.log10([self.reddening_corrected_lines['Hbeta_EW']])[0]

            strong_metallicity = measure_metallicity(log_O2.n, log_O2.s,
                                                     log_O3.n, log_O3.s,
                                                     log_Hbeta_EW.n, log_Hbeta_EW.s)

            self.metallicity = strong_metallicity
