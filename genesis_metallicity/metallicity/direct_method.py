import warnings
import numpy as np
from copy import deepcopy
from numpy.ma.core import maximum
import pyneb as pn
from uncertainties import ufloat
from uncertainties import unumpy as unp

from ..temperature.temperature_estimator import measure_temperature

warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in log10')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in sqrt')

pn.atomicData.setDataFile('o_iii_coll_Pal12-AK99.dat')

#####################
# Metallicity Class #
#####################

class METALLICITY:

    def __init__(self, object, data_dict, t2_calibration='L24', global_den=100, print_progress=False):

        #---------------------------------#
        #---- reading the line fluxes ----#
        #---------------------------------#

        self.object = object

        try:
            self.O3727  = data_dict['OII']/2
        except:
            self.O3727  = ufloat(np.nan, np.nan)

        try:
            self.O3729  = data_dict['OII']/2
        except:
            self.O3729  = ufloat(np.nan, np.nan)

        try:
            self.O4363  = data_dict['O4363']
        except:
            self.O4363  = ufloat(np.nan, np.nan)

        try:
            self.Hb     = data_dict['Hbeta']
        except:
            self.Hb     = ufloat(np.nan, np.nan)

        try:
            self.O4959  = data_dict['O4959']
        except:
            self.O4959  = ufloat(np.nan, np.nan)

        try:
            self.O5007  = data_dict['O5007']
        except:
            self.O5007  = ufloat(np.nan, np.nan)

        try:
            self.O7320  = data_dict['OII7320']
        except:
            self.O7320  = ufloat(np.nan, np.nan)

        #--------------------------------------#
        #---- calculating the Te_OII_O7320 ----#
        #--------------------------------------#

        self.Te_OII            = ufloat(np.nan, np.nan)
        self.Te_OII_Izotov     = ufloat(np.nan, np.nan)
        self.Te_OII_Langeroodi = ufloat(np.nan, np.nan)
        self.Te_OII_O7320      = ufloat(np.nan, np.nan)

        try:
            O2           = pn.Atom('O', '2')

            OII_ratio    = (self.O3727+self.O3729) / self.O7320
            OII_ratio    = np.array([OII_ratio.n-OII_ratio.s, OII_ratio.n, OII_ratio.n+OII_ratio.s])
            Te_OII_O7320 = O2.getTemDen(OII_ratio, wave1=3727, wave2=7320, den=global_den)
            Te_OII_O7320 = ufloat(min(Te_OII_O7320[1], 3e+4), np.abs(np.mean(np.diff(Te_OII_O7320))))

            if (Te_OII_O7320.n + Te_OII_O7320.s) > 3e+4:
                Te_OII_O7320 = ufloat(Te_OII_O7320.n, 3e+4-Te_OII_O7320.n)

            self.Te_OII_O7320 = Te_OII_O7320
        except:
            self.Te_OII_O7320 = ufloat(np.nan, np.nan)

        #----------------------------------------------------------------#
        #---- function for calculating the direct method metallicity ----#
        #----------------------------------------------------------------#

        def calculate_direct_metallicity(branch):

            #---- create atom objects for OII and OIII ----#

            O2 = pn.Atom('O', '2')
            O3 = pn.Atom('O', '3')

            #---- calculate Te(OIII) ----#

            OIII_ratio = self.O4363 / self.O5007
            OIII_ratio = np.array([OIII_ratio.n-OIII_ratio.s, OIII_ratio.n, OIII_ratio.n+OIII_ratio.s])
            OIII_ratio = np.clip(OIII_ratio, a_min=None, a_max=0.0465)
            Te_OIII    = O3.getTemDen(OIII_ratio, wave1=4363, wave2=5007, den=global_den)
            Te_OIII    = ufloat(Te_OIII[1], np.mean(np.diff(Te_OIII)))

            if (Te_OIII.n + Te_OIII.s) > 3e+4:
                Te_OIII = ufloat(Te_OIII.n, 3e+4-Te_OIII.n)

            #---- calculate Te(OII) (Izotov+2006) ----#

            def calculate_Te_OII(Te_OIII, branch):

                t = 1e-4 * Te_OIII

                if branch == 'low_Z':
                    t_OII = -0.577 + t * (2.065 - 0.498*t)

                if branch == 'intermediate_Z':
                    t_OII = -0.744 + t * (2.338 - 0.610*t)

                if branch == 'high_Z':
                    t_OII = 2.967 + t * (-4.797 + 2.827*t)

                Te_OII = 1e+4 * t_OII
                return Te_OII

            Te_OII_Izotov = calculate_Te_OII(Te_OIII, branch)

            if Te_OII_Izotov.n > 3e+4:
                Te_OII_Izotov = ufloat(3e+4, 0)
            if (Te_OII_Izotov.n + Te_OII_Izotov.s) > 3e+4:
                Te_OII_Izotov = ufloat(Te_OII_Izotov.n, 3e+4-Te_OII_Izotov.n)

            self.Te_OII_Izotov = Te_OII_Izotov

            #---- calculate Te(OII) (Langeroodi+2024) ----#

            O2_ratio_asymunc = (self.O3727+self.O3729)/self.Hb
            O3_ratio_asymunc = self.O5007/self.Hb

            O2_ratio         = ufloat(O2_ratio_asymunc.n, (O2_ratio_asymunc.s+O2_ratio_asymunc.s)/2)
            O3_ratio         = ufloat(O3_ratio_asymunc.n, (O3_ratio_asymunc.s+O3_ratio_asymunc.s)/2)

            O2_ratio         = unp.log10([O2_ratio])[0]
            O3_ratio         = unp.log10([O3_ratio])[0]

            Te_OII_Langeroodi = measure_temperature(O2_ratio.n, O2_ratio.s,
                                                    O3_ratio.n, O3_ratio.s,
                                                    Te_OIII.n/1e+4, Te_OIII.s/1e+4)
            Te_OII_Langeroodi = 1e+4 * Te_OII_Langeroodi
            self.Te_OII_Langeroodi = Te_OII_Langeroodi

            #---- choosing a Te(OII) measurement ----#

            if ~np.isnan(self.Te_OII_O7320.n):
                Te_OII = deepcopy(self.Te_OII_O7320)

            if np.isnan(self.Te_OII_O7320.n):

                if (t2_calibration == 'L24') and (8500 < Te_OIII.n < 14000):
                    Te_OII = deepcopy(self.Te_OII_Langeroodi)

                else:
                    Te_OII = deepcopy(self.Te_OII_Izotov)

            #---- measuring the abundances ----#

            OIII4959_ratio    = self.O4959 / self.Hb
            OIII4959_ratio    = np.array([OIII4959_ratio.n-OIII4959_ratio.s, OIII4959_ratio.n, OIII4959_ratio.n+OIII4959_ratio.s])
            OPP4959_abundance = O3.getIonAbundance(OIII4959_ratio, tem=[Te_OIII.n-Te_OIII.s, Te_OIII.n, Te_OIII.n+Te_OIII.s], den=[global_den,global_den,global_den], wave=4959, Hbeta=1)
            OPP4959_abundance = ufloat(OPP4959_abundance[1], np.mean(np.abs(np.diff(OPP4959_abundance))))

            OIII5007_ratio    = self.O5007 / self.Hb
            OIII5007_ratio    = np.array([OIII5007_ratio.n-OIII5007_ratio.s, OIII5007_ratio.n, OIII5007_ratio.n+OIII5007_ratio.s])
            OPP5007_abundance = O3.getIonAbundance(OIII5007_ratio, tem=[Te_OIII.n-Te_OIII.s, Te_OIII.n, Te_OIII.n+Te_OIII.s], den=[global_den,global_den,global_den], wave=5007, Hbeta=1)
            OPP5007_abundance = ufloat(OPP5007_abundance[1], np.mean(np.abs(np.diff(OPP5007_abundance))))

            OII3727_ratio     = self.O3727 / self.Hb
            OII3727_ratio     = np.array([OII3727_ratio.n-OII3727_ratio.s, OII3727_ratio.n, OII3727_ratio.n+OII3727_ratio.s])
            OP3727_abundance  = O2.getIonAbundance(OII3727_ratio, tem=[Te_OII.n-Te_OII.s, Te_OII.n, Te_OII.n+Te_OII.s], den=[global_den,global_den,global_den], wave=3727, Hbeta=1)
            OP3727_abundance  = ufloat(OP3727_abundance[1], np.mean(np.abs(np.diff(OP3727_abundance))))

            OII3729_ratio     = self.O3729 / self.Hb
            OII3729_ratio     = np.array([OII3729_ratio.n-OII3729_ratio.s, OII3729_ratio.n, OII3729_ratio.n+OII3729_ratio.s])
            OP3729_abundance  = O2.getIonAbundance(OII3729_ratio, tem=[Te_OII.n-Te_OII.s, Te_OII.n, Te_OII.n+Te_OII.s], den=[global_den,global_den,global_den], wave=3729, Hbeta=1)
            OP3729_abundance  = ufloat(OP3729_abundance[1], np.mean(np.abs(np.diff(OP3729_abundance))))

            # O_abundance = (OPP4959_abundance + OPP5007_abundance + OP3727_abundance + OP3729_abundance)/2
            O_abundance = OPP5007_abundance + (OP3727_abundance + OP3729_abundance)/2

            Z = 12 + unp.log10(O_abundance)

            # for debugging
            if False:
                print('---')
                print(' -> branch %s' %branch)
                print(Te_OII, Te_OIII, OP3727_abundance, OPP5007_abundance, Z)

            return Z, Te_OII, Te_OIII

        #----------------------------------------------------------------------------------------------#
        #---- function for checking if the returned metallicity and the used branch are consistent ----#
        #----------------------------------------------------------------------------------------------#

        def check_branch(Z, Te_OIII, branch, tolerate='no'):

            #---- do NOT account for Z uncertainties in determining the branch ----#

            if tolerate not in ['unc', 'max']:

                if (Z < 7.4) and (branch == 'low_Z'):
                    if print_progress:
                        print('--------------------------------------------------')
                        print(' + object %s' %object)
                        print('---------------')
                        print('Te([OIII])    = %.1f +/- %.1f' %(Te_OIII.n/1000, Te_OIII.s/1000))
                        print('12 + log(O/H) = %.2f +/- %.2f' %(Z.n, Z.s))
                        print('branch        : low_Z')
                    return True

                if (7.4 <= Z < 7.9) and (branch == 'intermediate_Z'):
                    if print_progress:
                        print('--------------------------------------------------')
                        print(' + object %s' %object)
                        print('---------------')
                        print('Te([OIII])    = %.1f +/- %.1f' %(Te_OIII.n/1000, Te_OIII.s/1000))
                        print('12 + log(O/H) = %.2f +/- %.2f' %(Z.n, Z.s))
                        print('branch        : intermediate_Z')
                    return True

                if (7.9 <= Z) and (branch == 'high_Z'):
                    if print_progress:
                        print('--------------------------------------------------')
                        print(' + object %s' %object)
                        print('---------------')
                        print('Te([OIII])    = %.1f +/- %.1f' %(Te_OIII.n/1000, Te_OIII.s/1000))
                        print('12 + log(O/H) = %.2f +/- %.2f' %(Z.n, Z.s))
                        print('branch        : high_Z')
                    return True

            #---- account for Z uncertainties in determining the branch ----#

            if tolerate in ['unc', 'max']:

                Z_tolerance = min(Z.s, 0.15)
                if tolerate == 'max': Z_tolerance = 0.15

                if ((Z-Z_tolerance) < 7.4) and (branch == 'low_Z'):
                    if print_progress:
                        print('--------------------------------------------------')
                        print(' + object %s' %object)
                        print('---------------')
                        print('Te([OIII])    = %.1f +/- %.1f' %(Te_OIII.n/1000, Te_OIII.s/1000))
                        print('12 + log(O/H) = %.2f +/- %.2f' %(Z.n, Z.s))
                        print('branch        : low_Z')
                    return True

                if (7.4 <= (Z-Z_tolerance) < 7.9) and (branch == 'intermediate_Z'):
                    if print_progress:
                        print('--------------------------------------------------')
                        print(' + object %s' %object)
                        print('---------------')
                        print('Te([OIII])    = %.1f +/- %.1f' %(Te_OIII.n/1000, Te_OIII.s/1000))
                        print('12 + log(O/H) = %.2f +/- %.2f' %(Z.n, Z.s))
                        print('branch        : intermediate_Z')
                    return True

                if (7.4 <= (Z+Z_tolerance) < 7.9) and (branch == 'intermediate_Z'):
                    if print_progress:
                        print('--------------------------------------------------')
                        print(' + object %s' %object)
                        print('---------------')
                        print('Te([OIII])    = %.1f +/- %.1f' %(Te_OIII.n/1000, Te_OIII.s/1000))
                        print('12 + log(O/H) = %.2f +/- %.2f' %(Z.n, Z.s))
                        print('branch        : intermediate_Z')
                    return True

                if (7.9 <= (Z+Z_tolerance)) and (branch == 'high_Z'):
                    if print_progress:
                        print('--------------------------------------------------')
                        print(' + object %s' %object)
                        print('---------------')
                        print('Te([OIII])    = %.1f +/- %.1f' %(Te_OIII.n/1000, Te_OIII.s/1000))
                        print('12 + log(O/H) = %.2f +/- %.2f' %(Z.n, Z.s))
                        print('branch        : high_Z')
                    return True

            return False

        #------------------------------------------------#
        #---- function to iterate over the above two ----#
        #------------------------------------------------#

        def iterate():

            branch = 'low_Z'
            Z, Te_OII, Te_OIII = calculate_direct_metallicity(branch)
            if check_branch(Z, Te_OIII, branch):
                return Z, Te_OII, Te_OIII

            branch = 'intermediate_Z'
            Z, Te_OII, Te_OIII = calculate_direct_metallicity(branch)
            if check_branch(Z, Te_OIII, branch):
                return Z, Te_OII, Te_OIII

            branch = 'high_Z'
            Z, Te_OII, Te_OIII = calculate_direct_metallicity(branch)
            if check_branch(Z, Te_OIII, branch):
                return Z, Te_OII, Te_OIII

            branch = 'low_Z'
            Z, Te_OII, Te_OIII = calculate_direct_metallicity(branch)
            if check_branch(Z, Te_OIII, branch, tolerate='unc'):
                return Z, Te_OII, Te_OIII

            branch = 'intermediate_Z'
            Z, Te_OII, Te_OIII = calculate_direct_metallicity(branch)
            if check_branch(Z, Te_OIII, branch, tolerate='unc'):
                return Z, Te_OII, Te_OIII

            branch = 'high_Z'
            Z, Te_OII, Te_OIII = calculate_direct_metallicity(branch)
            if check_branch(Z, Te_OIII, branch, tolerate='unc'):
                return Z, Te_OII, Te_OIII

            branch = 'low_Z'
            Z, Te_OII, Te_OIII = calculate_direct_metallicity(branch)
            if check_branch(Z, Te_OIII, branch, tolerate='max'):
                return Z, Te_OII, Te_OIII

            branch = 'intermediate_Z'
            Z, Te_OII, Te_OIII = calculate_direct_metallicity(branch)
            if check_branch(Z, Te_OIII, branch, tolerate='max'):
                return Z, Te_OII, Te_OIII

            branch = 'high_Z'
            Z, Te_OII, Te_OIII = calculate_direct_metallicity(branch)
            if check_branch(Z, Te_OIII, branch, tolerate='max'):
                return Z, Te_OII, Te_OIII

            raise Exception('no metallicty to return!')
            return

        #---- calculating the metallicity ----#

        try:
            Z, Te_OII, Te_OIII = iterate()
            if (Te_OII.n > -np.inf) and (Te_OIII.n < np.inf):
                self.Te_OII      = Te_OII
                self.Te_OIII     = Te_OIII
                self.metallicity = Z
            else:
                self.Te_OII      = ufloat(np.nan, np.nan)
                self.Te_OIII     = ufloat(np.nan, np.nan)
                self.metallicity = ufloat(np.nan, np.nan)
        except:
            self.Te_OII      = ufloat(np.nan, np.nan)
            self.Te_OIII     = ufloat(np.nan, np.nan)
            self.metallicity = ufloat(np.nan, np.nan)
