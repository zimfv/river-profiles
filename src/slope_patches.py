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
        """
        Returns the number of patches N
        """
        return len(self.patch_starts)
    
    
    def get_slopes(self):
        """
        Returns the array of slopes for each patch; 
        Associated with the equation (10) form the article by Leigh Royden and J. Taylor Perron
        """
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
    
    
    def get_elevations_for_patches(self, tau, chi, index=None, filter_outer=True):
        """
        Returns the elevation for each patch for moments tau in spatial points chi for patches given by index
        Associated with the equation (B4) form the article by Leigh Royden and J. Taylor Perron
        
        If index is None
        Returns the elevation for each patch for moments tau in spatial points chi
        
        Parameters:
        -----------
        tau : float or float array
            The (dimensionless) time argument
        
        chi : float or float array
            The (dimensionless) distance argument
            tau and chi should be the same shape
        
        index : int, int array or None
            The patches indices
            Should be the same shape as tau and chi
            If it's None, then this will be changed to array shape (N, *tau/chi.shape)
            where index[i] == i
            
        filter_outer : bool
            Remove the elements out of theoretical patch borders, if it's True
            
        Returns:
        --------
        lam : float array shape same as index.shape
            The elevation throught the time moments tau, spatial points for the patch from patch_index 
        """
        tau = np.array(tau)
        chi = np.array(chi)
        if index is None:
            index = np.arange(self.count(), dtype=int).reshape(self.count(), *np.ones(max(tau.ndim, chi.ndim), dtype=int))
        index = np.array(index)*np.ones(tau.shape, dtype=int)*np.ones(chi.shape, dtype=int)
        tau = tau*np.ones(index.shape)
        chi = chi*np.ones(index.shape)
        
        if filter_outer:
            lefts = self.get_lefts(tau)
            rights = self.get_rights(tau)
            for i in range(self.count()):
                chi[(index == i)&np.invert(chi >= lefts[i])] = np.nan
                chi[(index == i)&np.invert(chi <= rights[i])] = np.nan
        
        lam_prechi = chi*self.get_slopes()[index]
        lam_pretau = np.zeros(tau.shape)
        for i in range(self.count()):
            stair_values = np.append(0, self.uplift_rates - self.uplift_rates[i])
            lam_pretau[index == i] = step_integral(x0=self.patch_starts[i], 
                                                   x1=tau[index == i], 
                                                   borders=self.patch_starts, 
                                                   values=stair_values)
        lam = lam_prechi + lam_pretau
        
        return lam
    
    