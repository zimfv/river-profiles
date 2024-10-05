import numpy as np
from bguriver.math import step_function, step_integral, solve_bisect
from bguriver.warning_handler import get_runtime_warning_act_decorator

# set ignoring RuntimeWarnings in few methods
runtime_warning_act = get_runtime_warning_act_decorator(runtime_warning_action="ignore")


class SlopePatches:
    r"""
    The slope patches are needed to understand the evolution of river profiles, given by formula
    
    .. math::
        \cfrac{\partial\lambda}{\partial\tau} + (\cfrac{\partial\lambda}{\partial\chi})^n = \nu(\tau, \chi)
    
    where the function :math:`\nu` is :math:`\chi`-independent step function:
    
    .. math::
        \nu(\tau, \chi) = \nu(\tau) = 
        \begin{cases}
        \nu_0, \; \tau_0 \le \tau < \tau_1 \\
        \nu_1, \; \tau_1 \le \tau < \tau_2 \\
        ... \\
        \nu_i, \; \tau_{i} \le \tau < \tau_{i+1} \\
        ... \\
        \nu_N, \; \tau_N \le \tau < \infty \\
        \end{cases}

    we can just say, that :math:`\nu_{N+1} = \infty`.

    The concept of slope patches is described in the article by Leigh Royden and J. Taylor Perron:
    https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrf.20031


    Attributes:
    -----------
    patch_starts - float array length N
        The times :math:`\tau_i`, when :math:`i`-th patch starts
        
        Let N be the number of patches; It's not used in code itself
        
    uplift_rates - float array length N
        The values :math:`\nu_i`, uplift rate of the :math:`i`-th patch 
        (correspondes :math:`\tau` between :math:`\tau_i` and :math:`\tau_{i+1}`)
        
    n - float
        The constant :math:`n` - the exponent on channel slope

    """
    def __init__(self, patch_starts, uplift_rates, n):
        r"""
        Atributes:
        ----------
        patch_starts - float array length N
            The times :math:`\tau_i`, when :math:`i`-th patch starts
        
        uplift_rates - float array length N
            The values :math:`\nu_i`, uplift rate of the :math:`i`-th patch 
            (correspondes :math:`\tau` between :math:`\tau_i` and :math:`\tau_{i+1}`)

        n - float
            The constant :math:`n` - the exponent on channel slope
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
        r"""
        Returns the number of patches N
        """
        return len(self.patch_starts)
    
    
    def get_slopes(self):
        r"""
        Returns the array of slopes for each patch.

        Associated with the equation (13) form the article by Leigh Royden and J. Taylor Perron:

        .. math::
            \sigma(\chi, \tau) = \frac{\partial \lambda}{\partial \chi} = \nu(\tau)^{1/n}
        .. math::
            \sigma_i = \nu_i^{1/n}

        Returns:
        --------
        slopes - float array length N 
            The slopes of the patches
        """
        slopes = self.uplift_rates**(1/self.n)
        return slopes
    
    
    def get_rights(self, tau, first_is_infinite=True):
        r"""
        Returns the right spatial borders of the patches for time moments tau.  

        Associated with the equation (B7) form the article by Leigh Royden and J. Taylor Perron

        .. math::
            \chi_{R, i}(\tau) = n\nu_i^{(n - 1)/n}(\tau_i)


        Comment:
        --------
        It's nan for moments tau less than patch starts
        
        Parameters:
        -----------
        tau : float or float array
            :math:`\tau` - the dimensionless time argument
        
        first_is_infinite : bool, default True
            If it's True, the first patch right border will be at infinite

            If it's False, it will be defined proportional to tau
            
        Returns:
        --------
        res: float array shape (N, tau.shape)
            Right borders :math:`\chi_{R, i}` for each patch and each :math:`\tau` moment in tau
        """
        tau = np.array(tau, dtype=float)
        
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
        r"""
        Returns the spatial lengths of patch.

        It's defined as the product of time lengths of patches and their translation speed:

        .. math::
            \Delta_i = \chi_{R, i}(\tau_{i+1}) = n\nu_i^{(n - 1)/n}(\tau_{i+1} - \tau_i)

        Returns:
        --------
        res: float array shape (N, *tau.shape)
            Lengths of the patches :math:`\Delta_i`
        """
        res = self.n*self.uplift_rates**((self.n - 1)/self.n) 
        res *= np.append((self.patch_starts[1:] - self.patch_starts[:-1]), np.inf)
        return res
    
    
    def get_lefts(self, tau):
        r"""
        Returns the right spatial borders of the patches for :math:`\tau` moments tau.

        It's defune as subtraction of patch right (noninfinite) borders and their lengths:
        
        .. math::
            \chi_{L, i}(\tau) = \chi_{R, i}(\tau) - \Delta_i

        Comment:
        --------
        It's nan for moments tau less than patch starts
        
        Parameters:
        -----------
        tau : float or float array
            :math:`\tau` - the dimensionless time argument
        
        Returns:
        --------
        res: float array shape (N, tau.shape)
            Left borders :math:`\chi_{L, i}` for each patch and each :math:`\tau` moment in tau
        """
        tau = np.array(tau, dtype=float)
        
        rights = self.get_rights(tau, first_is_infinite=False)
        lengths = self.get_lengths()
        
        shape_res = np.array(rights.shape)
        shape_len = shape_res.copy()
        shape_len[1:] = 1
        
        res = rights - lengths.reshape(shape_len)
        res[res < 0] = 0
        return res
    
    
    @runtime_warning_act
    def get_elevations_for_patches(self, tau, chi, index=None, filter_outer=True):
        r"""
        Returns the elevation for each patch for moments tau in spatial points chi for patches given by index.

        Associated with the equation (B4) form the article by Leigh Royden and J. Taylor Perron:

        .. math::
            \lambda_i(\tau, \chi) = \chi\sigma_i + \int\limits_{\tau_i}^\tau [\nu(t) - \nu_i]dt
        
        If index is None
        Returns the elevation for each patch for points :math:`(\tau, \chi)`, which does not correspond the patch
        
        Parameters:
        -----------
        tau : float or float array
            :math:`\tau` - the dimensionless time argument
        
        chi : float or float array
            :math:`\chi` - the dimensionless distance argument
            
            The arguments tau and chi should be the same shape
        
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
            The elevation :math:`\lambda_i(\tau, \chi)` 
            throught the time moments :math:`\tau`, spatial points :math:`\chi` 
            for the patch :math:`i` from patch_index 

        """
        tau = np.array(tau, dtype=float)
        chi = np.array(chi, dtype=float)
        if index is None:
            index = np.arange(self.count(), dtype=int).reshape(self.count(), *np.ones(max(tau.ndim, chi.ndim), dtype=int))
        index = np.array(index)*np.ones(tau.shape, dtype=int)*np.ones(chi.shape, dtype=int)
        tau = np.array(tau*np.ones(index.shape))
        chi = np.array(chi*np.ones(index.shape))
        
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
    
    
    @runtime_warning_act
    def get_right_elevations(self, tau, first_is_infinite=True):
        r"""
        Returns elevations :math:`\lambda_i(\tau, \chi_{R, i})` on the right borders of the patches for time moments :math:`\tau`.

        Comment:
        --------
        It's nan for moments tau less than patch starts
        
        Parameters:
        -----------
        tau : float or float array
            :math:`\tau` - the dimensionless time argument
        
        first_is_infinite : bool, default True
            If it's True, the first patch right border will be at infinite

            If it's False, it will be defined proportional to the moment tau
            
        Returns:
        --------
        lam : float array shape same as index.shape
            The elevation :math:`\lambda_i(\tau, \chi_{R, i})`
            throught the time moments :math:`\tau`
        """
        tau = np.array(tau, dtype=float)
        chi = self.get_rights(tau, first_is_infinite=first_is_infinite)
        
        lam_prechi = np.array([chi[i]*s for i, s in enumerate(self.get_slopes())])
        lam_pretau = np.zeros(lam_prechi.shape)
        for i in range(self.count()):
            step_values = np.append(0, self.uplift_rates - self.uplift_rates[i])
            lam_pretau[i] = step_integral(x0=self.patch_starts[i], 
                                          x1=tau, 
                                          borders=self.patch_starts, 
                                          values=step_values)
        lam = lam_prechi + lam_pretau
        return lam
    
    
    @runtime_warning_act
    def get_left_elevations(self, tau):
        r"""
        Returns elevations :math:`\lambda_i(\tau, \chi_{L, i})` on the left borders of the patches for time moments :math:`\tau`.
        
        Comment:
        --------
        It's nan for moments tau less than patch starts
        
        Parameters:
        -----------
        tau : float or float array
            :math:`\tau` - the dimensionless time argument
        
        Returns:
        --------
        lam : float array shape same as index.shape
            The elevation :math:`\lambda_i(\tau, \chi_{L, i})` 
            throught the time moments :math:`\tau`
        """
        tau = np.array(tau, dtype=float)
        chi = self.get_lefts(tau)
        
        lam_prechi = np.array([chi[i]*s for i, s in enumerate(self.get_slopes())])
        lam_pretau = np.zeros(lam_prechi.shape)
        for i in range(self.count()):
            step_values = np.append(0, self.uplift_rates - self.uplift_rates[i])
            lam_pretau[i] = step_integral(x0=self.patch_starts[i], 
                                          x1=tau, 
                                          borders=self.patch_starts, 
                                          values=step_values)
        lam = lam_prechi + lam_pretau
        return lam
    
    
    @runtime_warning_act
    def get_elevations_for_stretch_zones(self, tau, chi, index=None, filter_outer=True):
        r"""
        Returns the elevation for each patch for moments tau in spatial points chi distinct for each patch.

        Associated with the equation (16b) form the article by Leigh Royden and J. Taylor Perron:

        .. math::
            \lambda_i = 
            \cfrac{n-1}{n}\left(\cfrac{\chi^n}{n(\tau - \tau_{i+1})}\right)^{\cfrac{1}{n-1}} + 
            \left(\tau - \tau_{i+1}\right)\nu_{i+1} + 
            \int\limits_{\tau_{i+1}}^\tau (\nu(t) - \nu_{i+1}) dt
        
        Returns the elevation over connection between :math:`i`-th and :math:`(i+1)`-th slope patches.
        This is stretch zones or consuming knick points
        
        Comment:
        --------
        It's nan in the case, if the point :math:`(\tau, \chi)` do not correspond the stretch zone.
        
        Parameters:
        -----------
        tau : float or float array
            :math:`\tau` - the dimensionless time argument
        
        chi : float or float array
            :math:`\chi` - the dimensionless distance argument
            
            The arguments tau and chi should be the same shape
        
        index : int or int array
            The stretch zones indices

            Should be the same shape as tau and chi

            If it's None, then this will be changed to array shape (N - 1, *tau/chi.shape)
            where index[i] == i
            
        filter_outer : bool
            Remove elevations out of theoretical patch borders, if it's True
            
        Returns:
        --------
        lam : float array shape tau/chi/inidex.shape
            The elevation throught the time moments tau, spatial points chi for the stretch zone from zone_inidex 
        """
        tau = np.array(tau, dtype=float)
        chi = np.array(chi, dtype=float)
        if index is None:
            index = np.arange(self.count() - 1, dtype=int).reshape(self.count() - 1, *np.ones(max(tau.ndim, chi.ndim), dtype=int))
        index = np.array(index)*np.ones(tau.shape, dtype=int)*np.ones(chi.shape, dtype=int)
        tau = np.array(tau*np.ones(index.shape))
        chi = np.array(chi*np.ones(index.shape))
        
        if self.n == 1:
            # here I returns nans matrix, but maybe its better to return right lams when chi is equal
            # lam = np.nans*np.ones(shape); lam[chi1 == chi] = lam1[chi1 == chi]; return lam
            lam = np.nan*tau*chi*index
            return lam
        
        if filter_outer:
            lefts, rights = self.get_rights(tau)[1:], self.get_lefts(tau)[:-1]
            cond = lefts > rights
            lefts[cond], rights[cond] = rights[cond], lefts[cond]
            for i in range(self.count() - 1):
                chi[(index == i)&np.invert(chi >= lefts[i])] = np.nan
                chi[(index == i)&np.invert(chi <= rights[i])] = np.nan
        
        tau1 = self.patch_starts[1:][index]
        lam = np.array((self.n - 1)/self.n*(chi**self.n/(tau - tau1)/self.n)**(1/(self.n - 1)))
        lam += (tau - tau1)*self.uplift_rates[1:][index]
        
        for i in range(self.count() - 1):
            stair_values = np.append(0, self.uplift_rates - self.uplift_rates[i + 1])
            lam[index == i] += step_integral(x0=self.patch_starts[i + 1], 
                                             x1=tau[index == i], 
                                             borders=self.patch_starts, 
                                             values=stair_values)
        return lam
    
    
    def get_elevation(self, tau, chi):
        r"""
        Returns the elevation :math:`\lambda(\tau, \chi)` at the point :math:`\chi` at moment :math:`\tau`

        .. math::
            \lambda(\tau, \chi) = 
            \begin{cases}
            \min \lambda_j(\tau, \chi), \; n \ge 1 \\
            \max \lambda_j(\tau, \chi), \; n < 1
            \end{cases}

        where :math:`\lambda_j(\tau, \chi)` is from the set of functions, 
        corresponing the elevation of patches or neighbour patch connections.
        
        Parameters:
        -----------
        tau : float or float array
            :math:`\tau` - the dimensionless time argument
        
        chi : float or float array
            :math:`\chi` - the dimensionless distance argument
            
            The arguments tau and chi should be the same shape
        
        Returns:
        --------
        lam: float array shape  tau/chi.shape
            The elevation :math:`\lambda(\tau, \chi)`
        """
        lam = np.concatenate([self.get_elevations_for_patches(tau, chi), 
                              self.get_elevations_for_stretch_zones(tau, chi)], axis=0)
        if self.n < 1:
            lam[np.isnan(lam)] = -np.inf
            lam = lam.max(axis=0)
        else:
            lam[np.isnan(lam)] = +np.inf
            lam = lam.min(axis=0)
        return lam
    
    
    def get_elevation_index(self, tau, chi):
        r"""
        Returns the index of patch or stretch_zone corresponding the real elevation
        
        .. math::
            \text{index}(\tau, \chi) = 
            \begin{cases}
            \text{argmin}\; \lambda_j(\tau, \chi), \; n \ge 1 \\
            \text{argmax}\; \lambda_j(\tau, \chi), \; n < 1
            \end{cases}

        where :math:`\lambda_j(\tau, \chi)` is from the set of functions, 
        corresponing the elevation of patches or neighbour patch connections.

        Parameters:
        -----------
        tau : float or float array
            :math:`\tau` - the dimensionless time argument
        
        chi : float or float array
            :math:`\chi` - the dimensionless distance argument
            
            The arguments tau and chi should be the same shape
        
        Returns:
        --------
        index: int array shape  tau/chi.shape
            The indices of patches and stretch zones coresponding the real elevation

            The index :math:`i` less than `N` coresponds the :math:`i`-th patch
            
            The index :math:`i` equal or higher than `N` coresponds the :math:`(i - N)`-th stretch zone
        """
        lam = np.concatenate([self.get_elevations_for_patches(tau, chi), 
                              self.get_elevations_for_stretch_zones(tau, chi)], axis=0)
        if self.n < 1:
            lam[np.isnan(lam)] = -np.inf
            index = lam.argmax(axis=0)
        else:
            lam[np.isnan(lam)] = +np.inf
            index = lam.argmin(axis=0)
        return index
    
    
    @runtime_warning_act
    def get_intersections_of_patches(self, tau, filtration=True):
        r"""
        Returns the spatial position :math:`\chi` of intersections between patches.
        
        Parameters:
        -----------
        tau : float or float array
            :math:`\tau` - the dimensionless time argument
        
        Returns:
        --------
        chi : float arrray shape (N, N, *tau.shape)
            The matrices of spatial positions :math:`\chi` of intersection for :math:`\tau` moments in tau

            The element [i, j] corresponds the intersection between i-th and j-th patches
        """
        shape = np.concatenate([[self.count(), self.count()], tau.shape]).astype(int)
        shape0 = np.ones(len(shape), dtype=int)
        shape2 = np.ones(len(shape), dtype=int)
        shape0[0], shape2[1] = self.count(), self.count()
        permutation = np.arange(len(shape), dtype=int)
        permutation[:2] = [1, 0]
        
        chi0 = np.transpose(self.get_lefts(tau)*np.ones(shape), axes=permutation)
        chi2 = self.get_rights(tau, first_is_infinite=False)*np.ones(shape)
        lam0 = np.transpose(self.get_left_elevations(tau)*np.ones(shape), axes=permutation)
        lam2 = self.get_right_elevations(tau, first_is_infinite=False)*np.ones(shape)
        slopes0 = self.get_slopes().reshape(shape0)*np.ones(shape)
        slopes2 = self.get_slopes().reshape(shape2)*np.ones(shape)
        
        cond = slopes0 != slopes2
        chi = np.nan*np.ones(shape)
        chi[cond] = (lam0[cond] - lam2[cond] - chi0[cond]*slopes0[cond] + chi2[cond]*slopes2[cond])/(slopes2[cond] - slopes0[cond])
        
        return chi
    
    
    @runtime_warning_act
    def get_intersections_of_stretch_zones(self, tau):
        r"""        
        Returns the spatial positions :math:`\chi` of intersections between stretch zones.
        
        Parameters:
        -----------
        tau : float or float array
            :math:`\tau` - the dimensionless time argument
        
        Returns:
        --------
        chi : float arrray shape (N - 1, N - 1, *tau.shape)
            The matrices of spatial positions :math:`\chi` of intersection for :math:`\tau` moments in tau

            The element [i, j] corresponds the intersection between i-th and j-th stretch zone
        """
        tau = np.array(tau, dtype=float)
        shape = np.concatenate([(self.count() - 1, self.count() - 1), tau.shape]).astype(int)
        if self.n == 1:
            return np.nan*np.ones(shape)
        
        start_matrix_j = self.patch_starts[1:]*np.ones([self.count() - 1, self.count() - 1])
        start_matrix_j = start_matrix_j.reshape(np.concatenate([(self.count() - 1, self.count() - 1), 
                                                                np.ones(tau.ndim, dtype=int)]))
        uplift_matrix_j = self.uplift_rates[1:]*np.ones([self.count() - 1, self.count() - 1])
        uplift_matrix_j = uplift_matrix_j.reshape(np.concatenate([(self.count() - 1, self.count() - 1), 
                                                                  np.ones(tau.ndim, dtype=int)]))
        
        aj = self.n*(tau*np.ones(shape) - start_matrix_j)
        aj[aj <= 0] = np.nan
        aj = aj**(1/(1 - self.n)) * (self.n - 1)/self.n
        ak = np.swapaxes(aj, 0, 1)
        
        bj = (tau*np.ones(shape) - start_matrix_j)*uplift_matrix_j
        for i in range(self.count() - 1):
            step_values = np.append(0, self.uplift_rates - self.uplift_rates[i+1])
            bj[:, i] += step_integral(x0=self.patch_starts[i+1], 
                                      x1=tau, 
                                      borders=self.patch_starts, 
                                      values=step_values)
        bk = np.swapaxes(bj, 0, 1)
        chi = ((bk - bj)/(aj - ak))**((self.n - 1)/self.n)
        return chi
    
    
    @runtime_warning_act
    def get_intersections_of_patches_and_stretch_zones(self, tau, xtol=1e-8, maxiter=100):
        r"""
        Returns the spatial positions :math:`\chi` of intersections between patches and stretch zones
        
        Parameters:
        -----------
        tau : float or float array
            :math:`\tau` - the dimensionless time argument
        
        xtol: float
            Accuracy for the bisection solver

        maxiter: int
            Maximal number of iterations for the bisection solver
            
        Returns:
        --------
        chi : float arrray shape (N, N - 1, *tau.shape)
            The matrices of spatial positions :math:`\chi` of intersections for :math:`\tau` moments in tau

            The element [i, j] corresponds the intersection between i-th patch and j-th stretch zone
        """
        tau = np.array(tau, dtype=float)
        tau_big = np.ones([self.count(), self.count() - 1, *tau.shape])*tau
        
        patch_index = np.arange(self.count()).reshape(self.count(), *np.ones(tau.ndim + 1, dtype=int))
        patch_index = patch_index**np.ones(tau_big.shape, dtype=int)
        zone_index = np.arange(self.count() - 1).reshape(1, self.count() - 1, *np.ones(tau.ndim, dtype=int))
        zone_index = zone_index**np.ones(tau_big.shape, dtype=int)
        
        f_full = lambda chi: self.get_elevations_for_patches(tau_big, chi, patch_index) - self.get_elevations_for_stretch_zones(tau_big, chi, zone_index)
        
        patches_rights = self.get_rights(tau).reshape(self.count(), 1, *tau.shape)*np.ones(tau_big.shape)
        patches_lefts = self.get_lefts(tau).reshape(self.count(), 1, *tau.shape)*np.ones(tau_big.shape)
        zones_rights = self.get_lefts(tau)[:-1].reshape(1, self.count() - 1, *tau.shape)*np.ones(tau_big.shape)
        zones_lefts = self.get_rights(tau)[1:].reshape(1, self.count() - 1, *tau.shape)*np.ones(tau_big.shape)
        
        patches_rights[np.isnan(patches_rights)] = +np.inf
        patches_lefts[np.isnan(patches_lefts)] = 0
        zones_rights[np.isnan(zones_rights)] = +np.inf
        zones_lefts[np.isnan(zones_lefts)] = 0
        
        chi0 = np.max([patches_lefts, zones_lefts], axis=0)
        chi1 = np.min([patches_rights, zones_rights], axis=0)
        
        cond = (chi0 < chi1)&(np.sign(f_full(chi0))*np.sign(f_full(chi1)) <= 0)
        
        f_cond = lambda chi: self.get_elevations_for_patches(tau_big[cond], chi, patch_index[cond]) - self.get_elevations_for_stretch_zones(tau_big[cond], chi, zone_index[cond])
        chi_cond = solve_bisect(f_cond, chi0[cond], chi1[cond], xtol=xtol, maxiter=maxiter)
        chi = np.nan*tau_big
        chi[cond] = chi_cond
        return chi
    
    
    def get_patches_relisation_borders(self, tau, xtol=1e-8, maxiter=100):
        r"""
        Returns the borders of patch realisation
        
        The left border is the maximal intersection of the patch 
        with another patch or stretch zone after the current, or the left border of patch itself.

        The right border is the minimal intersection of the patch 
        with another patch or stretch zone before the current, or the right border of patch itself.

        Parameters:
        -----------
        tau: float or float array
            :math:`\tau` - the dimensionless time argument
        
        xtol: float
            Accuracy for the bisection solver

        maxiter: int
            Maximal number of iterations for the bisection solver
            
        Returns:
        --------
        left_borders : float array shape (N, *tau.shape)
            left_borders[i] corresponds the left realisation borders of the :math`i`-th patch
        
        right_borders : float array shape (N, *tau.shape)
            right_borders[i] corresponds the right realisation borders of the :math`i`-th patch
        """
        tau = np.array(tau, dtype=float)
        
        
        lefts = self.get_lefts(tau).reshape(self.count(), 1, *tau.shape)
        rights = self.get_rights(tau).reshape(self.count(), 1, *tau.shape)
        intersections_pp = self.get_intersections_of_patches(tau)
        intersections_pz = self.get_intersections_of_patches_and_stretch_zones(tau, 
                                                                               xtol=xtol, 
                                                                               maxiter=maxiter)
        # define and sort the points, which are potential patches realisation borders
        points = np.concatenate([lefts, intersections_pp, intersections_pz, rights], axis=1)
        points = np.sort(points, axis=1)
        points[np.where(points[:, 1:] == points[:, :-1])] = np.nan
        points = np.sort(points, axis=1)
        
        # change the infinity values to some big enough but finite
        points_finite = points.copy()
        points_finite[points == +np.inf] = 1 + np.max(points[(points != +np.inf)&np.invert(np.isnan(points))])
        
        # check if patches correspond the space between points
        between = 0.5*(points_finite[:, 1:] + points_finite[:, :-1])
        between = self.get_elevation_index(tau*np.ones(between.shape), between)
        between = between == np.arange(self.count(), dtype=int).reshape(self.count(), *np.ones(tau.ndim + 1, dtype=int))
        
        # filter the points, by corresponding the correct patches
        pb_cond = np.logical_or(np.concatenate([np.zeros([self.count(), 1, *tau.shape], dtype=bool),  between], axis=1), 
                                np.concatenate([between, np.zeros([self.count(), 1, *tau.shape], dtype=bool)], axis=1))
        pb_cond = np.invert(pb_cond)
        points[pb_cond] = np.nan
        
        borders_left = points.copy()
        borders_left[np.isnan(borders_left)] = +np.inf
        borders_left = borders_left.min(axis=1)
        
        borders_right = points.copy()
        borders_right[np.isnan(borders_right)] = -np.inf
        borders_right = borders_right.max(axis=1)
        
        nan_cond = (borders_left > borders_right)|np.isnan(borders_left)|np.isnan(borders_right)
        borders_left[nan_cond] = np.nan
        borders_right[nan_cond] = np.nan
        
        return borders_left, borders_right
    

    def get_stretch_zones_relisation_borders(self, tau, xtol=1e-8, maxiter=100):
        r"""
        Returns the borders of stretch zones realisation
        
        The left border is the maximal intersection of the stretch zone 
        with another patch or stretch zone after the current, or the left border of patch itself.

        The right border is the minimal intersection of the stretch zone 
        with another patch or stretch zone before the current, or the right border of patch itself.

        Parameters:
        -----------
        tau: float or float array
            :math:`\tau` - the dimensionless time argument
        
        xtol: float
            Accuracy for the bisection solver

        maxiter: int
            Maximal number of iterations for the bisection solver
            
        Returns:
        --------
        left_borders : float array shape (N - 1, *tau.shape)
            left_borders[i] corresponds the left realisation borders of the :math`i`-th stretch zone
        
        right_borders : float array shape (N - 1, *tau.shape)
            left_borders[i] corresponds the right realisation borders of the :math`i`-th stretch zone
        """
        tau = np.array(tau, dtype=float)
        
        lefts = self.get_rights(tau)[1:].reshape(self.count() - 1, 1, *tau.shape)
        rights = self.get_lefts(tau)[:-1].reshape(self.count() - 1, 1, *tau.shape)
        intersections_zz = self.get_intersections_of_stretch_zones(tau)
        intersections_pz = self.get_intersections_of_patches_and_stretch_zones(tau, 
                                                                               xtol=xtol, 
                                                                               maxiter=maxiter)
        intersections_pz = np.swapaxes(intersections_pz, 0, 1)
        
        # define and sort the points, which are potential stretch zones realisation borders
        points = np.concatenate([lefts, intersections_pz, intersections_zz, rights], axis=1)
        points = np.sort(points, axis=1)
        points[np.where(points[:, 1:] == points[:, :-1])] = np.nan
        points = np.sort(points, axis=1)
        
        # check if stretch zones correspond the space between points
        between = 0.5*(points[:, 1:] + points[:, :-1])
        between = self.get_elevation_index(tau*np.ones(between.shape), between)
        between = between == self.count() + np.arange(self.count() - 1, dtype=int).reshape(self.count() - 1, *np.ones(tau.ndim + 1, dtype=int))
        
        # filter the points, by corresponding the correct stretch zones
        zb_cond = np.logical_or(np.concatenate([np.zeros([self.count() - 1, 1, *tau.shape], dtype=bool),  between], axis=1), 
                                np.concatenate([between, np.zeros([self.count() - 1, 1, *tau.shape], dtype=bool)], axis=1))
        zb_cond = np.invert(zb_cond)
        points[zb_cond] = np.nan
        
        borders_left = points.copy()
        borders_left[np.isnan(borders_left)] = +np.inf
        borders_left = borders_left.min(axis=1)
        
        borders_right = points.copy()
        borders_right[np.isnan(borders_right)] = -np.inf
        borders_right = borders_right.max(axis=1)
        
        nan_cond = (borders_left > borders_right)|np.isnan(borders_left)|np.isnan(borders_right)
        borders_left[nan_cond] = np.nan
        borders_right[nan_cond] = np.nan
        
        return borders_left, borders_right
    
    
    def get_nu_value(self, tau, chi=np.nan, rate_before=None):
        r"""
        Returns the uplift rate (nu-value) for each moment :math:`\tau`.

        Just returns the uplift rate at :math:`\chi=0`, this is the value of the function:

        .. math::
            \nu(\tau, \chi=0) = \nu(\tau) = \nu(\tau, \chi)

        from the equation defining slope patches object:

        .. math::
            \cfrac{\partial\lambda}{\partial\tau} + (\cfrac{\partial\lambda}{\partial\chi})^n = \nu(\tau, \chi)
        
        Parameters:
        -----------
        tau : float or float array
            The time argument of the uplift rate function
        
        chi : float or float array
            The spatial argument of the uplift rate function

            The result value does not depend on this argument, 
            but sometimes this function should take 2 parameters.
        
        rate_before: float or None
            The uplift rate before the first patch starts.

            If it is None, then set rate_before same as rate of 1st patch (patch index 0).
        
        Returns:
        --------
        nu : float or float array
             The uplift rates for moments tau
        """
        if rate_before is None:
            rate_before = self.uplift_rates[0]
        tau = np.array(tau, dtype=float)
        chi = np.array(chi, dtype=float)
        if tau.ndim > chi.ndim:
            use_shape = tau.shape
        else:
            use_shape = chi.shape
            
        nu = step_function(tau, borders=self.patch_starts, values=np.append(rate_before, self.uplift_rates))
        nu = np.ones(use_shape)*nu
        return nu