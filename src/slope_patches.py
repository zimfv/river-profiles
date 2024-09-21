import numpy as np
from src.math import *


class SlopePatches:
    def __init__(self, patch_starts, uplift_rates, n):
        """
        Atributes:
        ----------
        patch_starts - float array length N > 0 (let N be the number of patches; It's not used in code itself)
            The times tau_i, when i-th patch starts
        
        uplift_rates - float array length N
            The values nu_i, uplift rate of the i-th patch (correspondes tau between tau_i and tau_{i+1})
        
        n - float
            The exponent on channel slope
        """
        if len(patch_starts) == 0:
            raise ValueError('The patch_starts should be not empty.')
        if len(patch_starts) != len(uplift_rates):
            msg = 'The arrays patch_starts and uplift_rates should be the same length. But:'
            msg += f'\nlen(patch_starts) = {len(patch_starts)}\nlen(uplift_rates) = {len(uplift_rates)}'
            raise ValueError(msg)
            
        self.patch_starts = np.array(patch_starts)
        self.uplift_rates = np.array(uplift_rates)
        self.n = float(n)
        
    
    def count(self):
        # returns the number of patches N
        return len(self.patch_starts)
    
    
    def get_slopes(self):
        # Returns the array of slopes for each patch; 
        # Associated with the equation (10) form the article by Leigh Royden and J. Taylor Perron
        return self.uplift_rates**(1/self.n)
    
    
    def get_rights(self, tau, first_is_infinite=True):
        """
        Returns the right spatial borders of the patches for time moments tau.  
        Associated with the equation (B7) form the article by Leigh Royden and J. Taylor Perron
        
        Comment:
        --------
        It's nan for moments tau less than patch starts
        
        Parameters:
        -----------
        tau : float or float array
            The (dimensionless) time argument
        
        first_is_infinite : bool, default True
            If it's True, the first patch right border will be at infinite
            If it's False, it will be defined proportional to the moment tau
            
        Returns:
        --------
        res: float array shape (N, tau.shape)
            Right borders for each patch and each moment tau
        """
        tau = np.array(tau)
        
        shape_res = np.append(self.count(), tau.shape).astype(int)
        shape_rat = shape_res.copy()
        shape_rat[1:] = 1
        shape_tau = shape_res.copy()
        shape_tau[0] = 1
        
        tau_use = tau.reshape(shape_tau)
        tau_use = tau_use - self.patch_starts.reshape(shape_rat)
        tau_use[tau_use < 0] = np.nan
        
        res = self.n*self.uplift_rates.reshape(shape_rat)**((self.n - 1)/self.n)*tau_use
        if first_is_infinite:
            res[0] = np.inf
        return res
    
    
    def get_lengths(self):
        """
        Returns the spatial lengths of patch.
        It's defined as the product of time lengths of patches and their translation speed
        """
        res = self.n*self.uplift_rates**((self.n - 1)/self.n) 
        res *= np.append((self.patch_starts[1:] - self.patch_starts[:-1]), np.inf)
        return res
    
    
    def get_lefts(self, tau):
        """
        Returns the right spatial borders of the patches for time moments tau.
        It's defune as subtraction of patch right (noninfinite) borders and their lengths
        
        Comment:
        --------
        It's nan for moments tau less than patch starts
        
        Parameters:
        -----------
        tau : float or float array
            The (dimensionless) time argument
        
        Returns:
        --------
        res: float array shape (N, tau.shape)
            Left borders for each patch and each moment tau
        """
        tau = np.array(tau)
        
        rights = self.get_rights(tau, first_is_infinite=False)
        lengths = self.get_lengths()
        
        shape_res = np.array(rights.shape)
        shape_len = shape_res.copy()
        shape_len[1:] = 1
        
        res = rights - lengths.reshape(shape_len)
        res[res < 0] = 0
        return res