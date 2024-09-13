import os
import warnings
import numpy as np
import pickle as pkl
from scipy import stats
from uncertainties import ufloat

warnings.filterwarnings("ignore", message="divide by zero encountered in scalar divide")
warnings.filterwarnings("ignore", message="invalid value encountered in scalar multiply")
warnings.filterwarnings("ignore", message="invalid value encountered in scalar divide")

dimensions = 3
percentile = stats.chi2.cdf(1, df=dimensions)

##########
# Config #
##########

load_presaved = True
test          = False

#####################
# Making the Kernel #
#####################

if load_presaved:

    genesis_metallicity_base_path = os.path.dirname(__file__)
    genesis_metallicity_base_path = os.path.dirname(genesis_metallicity_base_path)
    kernel_metallicity_path       = os.path.join(genesis_metallicity_base_path, 'data', 'kernel_metallicity.pkl')

    with open(kernel_metallicity_path, 'rb') as handle:
        kernel_metallicity = pkl.load(handle)

##########################################
# Function for Measuring the Metallicity #
##########################################

def measure_metallicity(O2, O2_unc,
                        O3, O3_unc,
                        Hbeta_EW, Hbeta_EW_unc,
                        length=3):

    #---- making the O2, O3, EW(Hb), Z matrix ----#

    o2_top = np.linspace(O2, O2+O2_unc, length)
    o2_bot = np.linspace(O2, O2-O2_unc, length)
    o2     = np.concatenate((o2_bot, o2_top))
    o2     = np.unique(o2)

    o3_top = np.linspace(O3, O3+O3_unc, length)
    o3_bot = np.linspace(O3, O3-O3_unc, length)
    o3     = np.concatenate((o3_bot, o3_top))
    o3     = np.unique(o3)

    hb_top = np.linspace(Hbeta_EW, Hbeta_EW+Hbeta_EW_unc, length)
    hb_bot = np.linspace(Hbeta_EW, Hbeta_EW-Hbeta_EW_unc, length)
    hb     = np.concatenate((hb_bot, hb_top))
    hb     = np.unique(hb)

    z      = np.arange(6.00, 10.01, 0.01)

    o2o2, o3o3, hbhb, zz = np.meshgrid(o2, o3, hb, z, indexing='ij')
    grid_array           = np.stack([o2o2, o3o3, hbhb, zz], axis=-1)
    grid_array           = grid_array.reshape(-1, 4)
    grid_array           = grid_array.T

    #---- making the weights matrix ----#

    o2_wht = stats.norm.pdf(o2, loc=O2, scale=O2_unc)
    o2_wht = o2_wht/o2_wht[length-1]

    o3_wht = stats.norm.pdf(o3, loc=O3, scale=O3_unc)
    o3_wht = o3_wht/o3_wht[length-1]

    hb_wht = stats.norm.pdf(hb, loc=Hbeta_EW, scale=Hbeta_EW_unc)
    hb_wht = hb_wht/hb_wht[length-1]

    z_wht  = np.ones(len(z))

    o2o2_wht, o3o3_wht, hbhb_wht, zz_wht = np.meshgrid(o2_wht, o3_wht, hb_wht, z_wht, indexing='ij')
    weight_array = np.stack([o2o2_wht, o3o3_wht, hbhb_wht, zz_wht], axis=-1)
    weight_array = np.prod(weight_array, axis=-1)
    weight_array = weight_array.reshape(-1)

    #---- calculating the PDF ----#

    pdf = kernel_metallicity.evaluate(grid_array)
    pdf = pdf * weight_array

    #---- marginalization ----#

    metallicity_array  = grid_array[-1,:]
    sorted_indices     = np.argsort(metallicity_array)

    metallicity_array  = metallicity_array[sorted_indices]
    pdf                = pdf[sorted_indices]

    pdf_normalized     = pdf / np.sum(pdf)
    pdf_metallicity    = metallicity_array

    #---- maximum likelihood value ----#

    ml_index           = np.argmax(pdf_normalized)
    ml_metallicity     = pdf_metallicity[ml_index]

    #---- 50th percentile value ----#

    cdf                = np.cumsum(pdf_normalized)
    ml_cdf             = cdf[ml_index]

    lo_index           = np.argmin(np.abs(cdf-(ml_cdf-percentile/2)))
    up_index           = np.argmin(np.abs(cdf-(ml_cdf+percentile/2)))

    output_array       = [pdf_metallicity[lo_index], ml_metallicity, pdf_metallicity[up_index]]
    output_metallicity = ufloat(ml_metallicity, np.mean(np.diff(output_array)))
    return output_metallicity
