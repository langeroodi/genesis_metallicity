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
    kernel_temperature_path       = os.path.join(genesis_metallicity_base_path, 'data', 'kernel_temperature.pkl')

    with open(kernel_temperature_path, 'rb') as handle:
        kernel_temperature = pkl.load(handle)

##################################
# Function for Estimating the T2 #
##################################

def measure_temperature(O2, O2_unc,
                        O3, O3_unc,
                        T3, T3_unc,
                        length=3):

    #---- making the O2, O3, t3, t2 matrix ----#

    o2_top = np.linspace(O2, O2+O2_unc, length)
    o2_bot = np.linspace(O2, O2-O2_unc, length)
    o2     = np.concatenate((o2_bot, o2_top))
    o2     = np.unique(o2)

    o3_top = np.linspace(O3, O3+O3_unc, length)
    o3_bot = np.linspace(O3, O3-O3_unc, length)
    o3     = np.concatenate((o3_bot, o3_top))
    o3     = np.unique(o3)

    t3_top = np.linspace(T3, T3+T3_unc, length)
    t3_bot = np.linspace(T3, T3-T3_unc, length)
    t3     = np.concatenate((t3_bot, t3_top))
    t3     = np.unique(t3)

    t2     = np.arange(0.6, 2.3, 0.01)

    o2o2, o3o3, t3t3, t2t2 = np.meshgrid(o2, o3, t3, t2, indexing='ij')
    grid_array             = np.stack([o2o2, o3o3, t3t3, t2t2], axis=-1)
    grid_array             = grid_array.reshape(-1, 4)
    grid_array             = grid_array.T

    #---- making the weights matrix ----#

    o2_wht = stats.norm.pdf(o2, loc=O2, scale=O2_unc)
    o2_wht = o2_wht/o2_wht[length-1]

    o3_wht = stats.norm.pdf(o3, loc=O3, scale=O3_unc)
    o3_wht = o3_wht/o3_wht[length-1]

    t3_wht = stats.norm.pdf(t3, loc=T3, scale=T3_unc)
    t3_wht = t3_wht/t3_wht[length-1]

    t2_wht = np.ones(len(t2))

    o2o2_wht, o3o3_wht, t3t3_wht, t2t2_wht = np.meshgrid(o2_wht, o3_wht, t3_wht, t2_wht, indexing='ij')
    weight_array = np.stack([o2o2_wht, o3o3_wht, t3t3_wht, t2t2_wht], axis=-1)
    weight_array = np.prod(weight_array, axis=-1)
    weight_array = weight_array.reshape(-1)

    #---- calculating the PDF ----#

    pdf = kernel_temperature.evaluate(grid_array)
    pdf = pdf * weight_array

    #---- marginalization ----#

    t2_array       = grid_array[-1,:]
    sorted_indices = np.argsort(t2_array)

    t2_array       = t2_array[sorted_indices]
    pdf            = pdf[sorted_indices]

    pdf_normalized = pdf / np.sum(pdf)
    pdf_t2         = t2_array

    #---- maximum likelihood value ----#

    ml_index       = np.argmax(pdf_normalized)
    ml_t2          = pdf_t2[ml_index]

    #---- 50th percentile value ----#

    cdf               = np.cumsum(pdf_normalized)
    ml_cdf            = cdf[ml_index]

    lo_index          = np.argmin(np.abs(cdf-(ml_cdf-percentile/2)))
    up_index          = np.argmin(np.abs(cdf-(ml_cdf+percentile/2)))

    output_array      = [pdf_t2[lo_index], ml_t2, pdf_t2[up_index]]
    output_t2         = ufloat(ml_t2, np.mean(np.diff(output_array)))
    return output_t2
