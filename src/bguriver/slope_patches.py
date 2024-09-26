import numpy as np
from bguriver.math import step_function, step_integral, solve_bisect


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
            step_values = np.append(0, self.uplift_rates - self.uplift_rates[i])
            lam_pretau[i] = step_integral(x0=self.patch_starts[i], 
                                          x1=tau, 
                                          borders=self.patch_starts, 
                                          values=step_values)
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
            step_values = np.append(0, self.uplift_rates - self.uplift_rates[i])
            lam_pretau[i] = step_integral(x0=self.patch_starts[i], 
                                          x1=tau, 
                                          borders=self.patch_starts, 
                                          values=step_values)
        lam = lam_prechi + lam_pretau
        return lam
    
    
    def get_elevations_for_stretch_zones(self, tau, chi, index=None, filter_outer=True):
        """
        Returns the elevation for each patch for moments tau in spatial points chi distinct for each patch
        Associated with the equation (16b) form the article by Leigh Royden and J. Taylor Perron
        
        
        Returns the elevation over connection between i-th and (i+1)-th slope patches.
        This is stretch zones or consuming knick points
        
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
        tau = np.array(tau)
        chi = np.array(chi)
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
        """
        Returns the index of patch or stretch_zone corresponding the real elevation
        
        Parameters:
        -----------
        tau : float or float array
            The (dimensionless) time argument
        
        chi : float or float array
            The (dimensionless) distance argument
            tau and chi should be the same shape
        
        Returns:
        --------
        index: int array shape  tau/chi.shape
            The indices of patches and stretch zones coresponding the real elevation
            The index i less than N coresponds the i-th patch
            The index i equal or higher than N coresponds the (i - N)-th stretch zone
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
    
    
    def get_intersections_of_patches(self, tau, filtration=True):
        """
        Returns the spatial position of intersections between patches
        
        Parameters:
        -----------
        tau : float or float array
            The (dimensionless) time argument
        
        Returns:
        --------
        chi : float arrray shape (N, N, *tau.shape)
            The matrices of spatial positions of intersection  for time moments tau
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
    
    
    def get_intersections_of_stretch_zones(self, tau):
        """
        Returns the spatial position of intersections between neighbour connections
        
        Parameters:
        -----------
        tau : float or float array
            The (dimensionless) time argument
        
        Returns:
        --------
        chi : float arrray shape (N - 1, N - 1, *tau.shape)
            The matrices of spatial positions of intersection  for time moments tau
            The element [i, j] corresponds the intersection between i-th and j-th neighbour connection
        """
        tau = np.array(tau)
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
    
    
    def get_intersections_of_patches_and_stretch_zones(self, tau, xtol=1e-8, maxiter=100):
        """
        Returns the spatial position of intersections between patches and stretch zones
        
        Parameters:
        -----------
        tau : float or float array
            The (dimensionless) time argument
        
        xtol: float
            Accuracy for the bisection solver

        maxiter: int
            Maximal number of iterations for the bisection solver
            
        Returns:
        --------
        chi : float arrray shape (N, N - 1, *tau.shape)
            The matrices of spatial positions of intersection  for time moments tau
            The element [i, j] corresponds the intersection between i-th patch and j-th stretch zone
        """
        tau = np.array(tau)
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
        """
        Returns the borders of patch realisation
        
        Parameters:
        -----------
        tau: float or float array
        
        xtol: float
            Accuracy for the bisection solver

        maxiter: int
            Maximal number of iterations for the bisection solver
            
        Returns:
        --------
        left_borders : float array shape (N, *tau.shape)
        
        right_borders : float array shape (N, *tau.shape)
        """
        tau = np.array(tau)
        
        
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
        """
        Returns the borders of stretch zones realisation
        
        Parameters:
        -----------
        tau: float or float array
        
        xtol: float
            Accuracy for the bisection solver

        maxiter: int
            Maximal number of iterations for the bisection solver
            
        Returns:
        --------
        left_borders : float array shape (N, *tau.shape)
        
        right_borders : float array shape (N, *tau.shape)
        """
        tau = np.array(tau)
        
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
        """
        Returns the uplift rate (nu-value) for each moment tau.
        Just returns the uplift rate at chi=0, nu(tau, chi=0)
        
        Parameters:
        -----------
        tau : float or float array
            The time argument of the uplift rate function
        
        chi : float or float array
            The spatial argument of the uplift rate function
            The result value does not depend on this argument, 
            but sometimes this function should take 2 paraeters.
        
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
        tau = np.array(tau)
        chi = np.array(chi)
        if tau.ndim > chi.ndim:
            use_shape = tau.shape
        else:
            use_shape = chi.shape
            
        nu = step_function(tau, borders=self.patch_starts, values=np.append(rate_before, self.uplift_rates))
        nu = np.ones(use_shape)*nu
        return nu