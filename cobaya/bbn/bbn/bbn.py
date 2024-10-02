import os
import numpy as np
from scipy.interpolate import interp2d, interp1d

from cobaya.likelihood import Likelihood
from cobaya.log import LoggedError

class BBN(Likelihood):
    """
    Likelihood that compares measured abundances of helium-4 and deuterium to BBN predictions.
    """

    def initialize(self):
        """
        Read the BBN table and initialize the interpolation functions.
        """
        ombh2_size = self.ombh2_size
        deltaneff_size = self.deltaneff_size
        
        if os.path.isabs(self.table_file):
            table_file = self.table_file
            self.path = os.path.dirname(table_file)
        else:
            self.path = self.path or self.get_class_path()
            if not self.path:
                raise LoggedError(self.log,
                                  "No path given for %s. Set the likelihood "
                                  "property 'path' or the common property.",
                                  self.table_file)

            table_file = os.path.normpath(os.path.join(self.path, self.table_file))

        self.log.info("Reading BBN table from %s", table_file)

        self.data = np.loadtxt(table_file)
        
        if len(self.data) != ombh2_size * deltaneff_size:
            raise LoggedError(self.log, "The data file %s does not have the expected format." % table_file)
        
        # params grid
        ombh2_data = np.array(self.data[:ombh2_size, 0])  # omegabh2
        deltaneff_data = np.array(self.data[::ombh2_size, 2])  # DeltaNeff
        ypbbn_data = np.array(self.data[:, 4]).reshape(deltaneff_size, ombh2_size)  # YpBBN
        dh_data = 1.e5 * np.array(self.data[:, 6]).reshape(deltaneff_size, ombh2_size)
        YpBBN = interp2d(ombh2_data, deltaneff_data, ypbbn_data, kind='cubic')
        DH = interp2d(ombh2_data, deltaneff_data, dh_data, kind='cubic')
        
        # uncertainties
        sigma_ypbbn_data = np.array(self.data[:, 5]).reshape(deltaneff_size, ombh2_size)
        sigma_dh_data = 1.e5 * np.array(self.data[:, 7]).reshape(deltaneff_size, ombh2_size)
        sigma_YpBBN = interp2d(ombh2_data, deltaneff_data, sigma_ypbbn_data, kind='cubic')
        sigma_DH = interp2d(ombh2_data, deltaneff_data, sigma_dh_data, kind='cubic')

        self.get_Yp = YpBBN
        self.get_DH = DH

        self.get_Yperr = sigma_YpBBN
        self.get_DHerr = sigma_DH

        self.ombh2_bounds = [ombh2_data[0], ombh2_data[-1]]
        self.dNeff_bounds = [deltaneff_data[0], deltaneff_data[-1]]


    def get_requirements(self):
        return {'omega_b': None, 'dneff': None}


    def logp(self, **params_values):
        ombh2, delta_neff = self.provider.get_param(['omega_b', 'dneff'])
        return self.log_likelihood(ombh2, delta_neff)


    def log_likelihood(self, ombh2, delta_neff):
        """
        Compute the likelihood of the BBN data given the input parameters.
        """
        yp = self.get_Yp(ombh2, delta_neff)
        dh = self.get_DH(ombh2, delta_neff)
        yperr = self.get_Yperr(ombh2, delta_neff)
        dherr = self.get_DHerr(ombh2, delta_neff)

        chi2 = 0.
        if self.include_He:
            chi2 += ((yp - self.yp_mean) ** 2) / (yperr ** 2 + self.yp_err ** 2)
        if self.include_DH:
            chi2 += ((dh - self.dh_mean) ** 2) / (dherr ** 2 + self.dh_err ** 2)

        return -0.5 * chi2