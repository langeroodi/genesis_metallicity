import numpy as np
from scipy import stats

class CustomKDE(stats.gaussian_kde):
    def __init__(self, dataset, covariance):
        super().__init__(dataset)
        self.covariance = covariance
        self.inv_cov = np.linalg.inv(covariance)
