import numpy as np


def lengths_check(borders, values, raise_error=True):
    # Returns True, if lengths are correct: borders array should contain 1 element less, than values array
    # If this is not correct, it will raise error if raise_error or return Flase in other case.
    r = len(values) - len(borders) == 1
    if raise_error and not r:
        msg = f'Wrong arrays lengths: borders array should contain 1 element less, than values array. But their length are {len(borders)} and {len(values)}.'
        raise ValueError(msg)
        
    borders = np.array(borders)
    if (borders[1:] - borders[:-1] <= 0).any():
        warnings.warn("The borders array is not strictly increasing")
    return r


def stair_function(x, borders=[], values=[0]):
    """
    Parameters:
    -----------
    x : float or float array
        The argument of the function
    
    borders : list of floats len N
        The points, when the function changes the value
    
    values: list of floats len N+1
        values[0] corresponds to the function value before borders[0]
        values[i] corresponds to the function value between borders[i-1] and borders[i]
        
    Returns:
    --------
    r : float or float array shape x.shape
    """
    lengths_check(borders, values, raise_error=True)
    x = np.array(x)
    if x.shape == ():
        r = values[0]*np.ones(1)
    else:
        r = values[0]*np.ones(x.shape)
    for i in range(len(borders)):
        r[x >= borders[i]] = values[i+1]
    if x.shape == ():
        return r[0]
    return r


def stair_integral(x0, x1, borders=[], values=[0], negative_backward=True):
    """
    Returns the value of the definite integral of a stair function defined by borders on the interval [x0, x1] 
    
    Parameters:
    -----------
    x0 : float or float array
        The 1st Interval border
        
    x0 : float or float array
        The 2nd Interval border
    
    borders : list of floats len N
        The points, when the stair function changes the value
    
    values: list of floats len N+1
        values[0] corresponds to the stair function value before borders[0]
        values[i] corresponds to the stair function value between borders[i-1] and borders[i]
        
    negative_backward: bool
        Is the integral symmetric: F(x0, x1) = -F(x1, x0)
        
    Returns:
    --------
    res : float or float array shape x.shape
    """
    lengths_check(borders, values, raise_error=True)
    
    x1 = np.array(x1)
    x0 = x0*np.ones(x1.shape)
    
    borders_ = np.concatenate([[-np.inf], borders, [+np.inf]])
    
    ones = np.ones(np.append(len(borders_), x1.shape).astype(int))
    res = ones*borders_.reshape(np.append(len(borders_), np.ones(len(x1.shape))).astype(int))
    res[res < x0] = (x0*ones)[res < x0]
    res[res > x1] = (x1*ones)[res > x1]
    res = res[1:] - res[:-1]
    res = res*np.array(values).reshape(np.append(len(values), np.ones(len(res.shape) - 1)).astype(int))
    res = res.sum(axis=0)
    
    if negative_backward and (x0 > x1).any():
        res[x0 > x1] = -stair_integral(x0=x1[x0 > x1], x1=x0[x0 > x1], borders=borders, values=values)
        
    return res


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
    
    
    def get_uplift_rate(self, tau, chi=np.nan, rate_before=np.nan):
        """
        Returns the upluft rate (nu-value) for each moment tau.
        
        Parameters:
        -----------
        tau : float or float array
            The time argument of the uplift rate function
        
        chi : float or float array
            The spatial argument of the uplift rate function
            The result value should not depend on this argument, 
            but sometimes this function should take 2 paraeters.
        
        rate_before: float
            The uplift rate before the first patch starts
        
        Returns:
        --------
        res : float or float array
            The uplift rates for moments tau
        """
        tau = np.array(tau)
        chi = np.array(chi)
        if tau.ndim > chi.ndim:
            use_shape = tau.shape
        else:
            use_shape = chi.shape
            
        res = stair_function(tau, borders=self.patch_starts, values=np.append(rate_before, self.uplift_rates))
        res = np.ones(use_shape)*res
        return res
    
    
    def get_elevations_for_patches(self, tau, chi):
        """
        Returns the elevation for each patch for moments tau in spatial points chi
        Associated with the equation (B4) form the article by Leigh Royden and J. Taylor Perron
        
        Comment:
        --------
        It's nan in the case, if the point (tau, chi) do not correspond the patch.
        
        Parameters:
        -----------
        tau : float or float array
            The (dimensionless) time argument
        
        chi : float or float array
            The (dimensionless) distance argument
            tau and chi should be the same shape
        
        Returns:
        --------
        lam : float array shape (N, tau/chi.shape)
            The elevation throught the time moments tau, spatial points chi for each patch
        """
        tau = np.array(tau)
        chi = np.array(chi)
        
        chi_left = self.get_lefts(tau)
        chi_right = self.get_rights(tau)
        chi_use = chi*np.ones(np.append(self.count(), chi.shape).astype(int))
        chi_use[chi_use < chi_left] = np.nan
        chi_use[chi_use > chi_right] = np.nan
        chi_use[np.isnan(chi_left)] = np.nan
        chi_use[np.isnan(chi_right)] = np.nan
        
        lam_prechi = chi_use*self.get_slopes().reshape(np.append(self.count(), np.ones(chi.ndim, dtype=int)))
        
        lam_pretau = np.zeros(lam_prechi.shape)
        for i in range(self.count()):
            stair_values = np.append(0, self.uplift_rates - self.uplift_rates[i])
            lam_pretau[i] = stair_integral(x0=self.patch_starts[i], 
                                           x1=tau, 
                                           borders=self.patch_starts, 
                                           values=stair_values)
        lam = lam_prechi + lam_pretau
        return lam
    
    
    def get_right_elevations(self, tau, first_is_infinite=True):
        """
        Returns elevations pon the right borders of the patches for time moments tau.  
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
        lam: float array shape (N, tau.shape)
            Elevations on the right borders for each patch and each moment tau
        """
        tau = np.array(tau)
        chi = self.get_rights(tau, first_is_infinite=first_is_infinite)
        
        lam_prechi = np.array([chi[i]*s for i, s in enumerate(self.get_slopes())])
        lam_pretau = np.zeros(lam_prechi.shape)
        for i in range(self.count()):
            stair_values = np.append(0, self.uplift_rates - self.uplift_rates[i])
            lam_pretau[i] = stair_integral(x0=self.patch_starts[i], 
                                           x1=tau, 
                                           borders=self.patch_starts, 
                                           values=stair_values)
        lam = lam_prechi + lam_pretau
        return lam
    
    
    def get_left_elevations(self, tau):
        """
        Returns elevations on the left borders of the patches for time moments tau.  
        Associated with the equation (B7) form the article by Leigh Royden and J. Taylor Perron
        
        Comment:
        --------
        It's nan for moments tau less than patch starts
        
        Parameters:
        -----------
        tau : float or float array
            The (dimensionless) time argument
        
        Returns:
        --------
        lam: float array shape (N, tau.shape)
            Elevations on the left borders for each patch and each moment tau
        """
        tau = np.array(tau)
        chi = self.get_lefts(tau)
        
        lam_prechi = np.array([chi[i]*s for i, s in enumerate(self.get_slopes())])
        lam_pretau = np.zeros(lam_prechi.shape)
        for i in range(self.count()):
            stair_values = np.append(0, self.uplift_rates - self.uplift_rates[i])
            lam_pretau[i] = stair_integral(x0=self.patch_starts[i], 
                                           x1=tau, 
                                           borders=self.patch_starts, 
                                           values=stair_values)
        lam = lam_prechi + lam_pretau
        return lam
    
    
    def get_neighbour_intersections(self, tau):
        """
        Returns the spatial position and elevation of intersections between neighbours patches
        
        Parameters:
        -----------
        tau : float or float array
            The (dimensionless) time argument
        
        Returns:
        --------
        chi1 : float arrray shape(N-1, tau.shape)
            THe spatial positions of intersection between i-th and (i+1) patches for time moments tau
            
        lam1 : float arrray shape(N-1, tau.shape)
            THe elevations of intersection between i-th and (i+1) patches for time moments tau
        """
        chi0 = self.get_lefts(tau)[:-1]
        chi2 = self.get_rights(tau)[1:]
        lam0 = self.get_left_elevations(tau)[:-1]
        lam2 = self.get_right_elevations(tau)[1:]
        slopes0 = self.get_slopes()[:-1]
        slopes2 = self.get_slopes()[1:]
        
        slopes_shape = np.ones(chi0.ndim, dtype=int)
        slopes_shape[0] = self.count() - 1
        slopes0 = slopes0.reshape(slopes_shape)
        slopes2 = slopes2.reshape(slopes_shape)
        
        chi1 = (lam0 - lam2 - chi0*slopes0 + chi2*slopes2)/(slopes2 - slopes0)
        lam1 = lam0 + (chi1 - chi0)*slopes0 
        return chi1, lam1
    
    
    def get_neighbour_connections(self, tau, chi):
        """
        Returns the elevation over connection between i-th and (i+1)-th slope patches.
        This is stretch zones or consuming knick points
        
        Parameters:
        -----------
        tau : float or float array
            The (dimensionless) time argument
        
        chi : float or float array
            The (dimensionless) distance argument
            tau and chi should be the same shape
        
        Returns:
        --------
        lam: float array shape (N-1, tau/chi.shape)
            The elevation over connection between i-th and (i+1)-th slope patches;
            Stretch zones or consuming knick points
        """
        if self.n == 1:
            # here I returns nans matrix, but maybe its better to return right lams when chi is equal
            if tau.ndim > chi.ndim:
                shape = tau.shape
            else:
                shape = chi.shape
            shape = np.append(self.count() - 1, shape).astype(int)
            return np.nan*np.ones(shape)
        chi0 = self.get_lefts(tau)[:-1]
        chi1 = self.get_rights(tau)[1:]
        lam0 = self.get_left_elevations(tau)[:-1]
        lam1 = self.get_right_elevations(tau)[1:]
        
        chi_min = np.array([chi0, chi1]).min(axis=0)
        chi_max = np.array([chi0, chi1]).max(axis=0)
        chi = np.ones(chi0.shape)*chi
        chi[np.isnan(chi0)] = np.nan
        chi[np.isnan(chi1)] = np.nan
        chi[chi < chi_min] = np.nan
        chi[chi > chi_max] = np.nan
        chi[chi_min == chi_max] = np.nan
        
        
        tau1_shape = np.ones(chi1.ndim, dtype=int)
        tau1_shape[0] = self.count() - 1
        tau1 = self.patch_starts[1:].reshape(tau1_shape)
        
        # equation 16b
        lam = (self.n - 1)/self.n*(chi**self.n/(tau - tau1)/self.n)**(1/(self.n - 1))
        lam += (tau - tau1)*self.uplift_rates[1:].reshape(tau1_shape)
        for i in range(self.count() - 1):
            stair_values = np.append(0, self.uplift_rates - self.uplift_rates[i + 1])
            lam[i] += stair_integral(x0=self.patch_starts[i + 1], 
                                     x1=tau, 
                                     borders=self.patch_starts, 
                                     values=stair_values)
        return lam
    
    
    def get_elevation(self, tau, chi):
        """
        Returns the elevation
        
        Parameters:
        -----------
        tau : float or float array
            The (dimensionless) time argument
        
        chi : float or float array
            The (dimensionless) distance argument
            tau and chi should be the same shape
        
        Returns:
        --------
        lam: float array shape  tau/chi.shape
            The elevation
        """
        tau, chi = np.array(tau), np.array(chi)
        if tau.ndim > chi.ndim:
            chi = chi*np.ones(tau.shape)
        else:
            tau = tau*np.ones(chi.shape)
        
        lam = np.concatenate([self.get_elevations_for_patches(tau, chi), 
                              self.get_neighbour_connections(tau, chi)], axis=0)
        if self.n < 1:
            lam[np.isnan(lam)] = -np.inf
            lam = lam.max(axis=0)
        else:
            lam[np.isnan(lam)] = +np.inf
            lam = lam.min(axis=0)
        return lam