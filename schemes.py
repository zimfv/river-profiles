import numpy as np
import scipy as sp

def approx1(nu, initial, border, n=1.0, dtau=1e-3, dchi=1e-3, ntau=200, nchi=1000):
    """
    Returns the solution of equation
    $$
    \cfrac{\partial\lambda}{\partial\tau} = \nu(\chi, \tau) - (\cfrac{\partial\lambda}{\partial\chi})^n
    $$
    using difference scheme:
    $$
    L_i^{j+1} = L_i^j + d\tau\cdot\nu(i\cdot d\chi, j\cdot d\tau) + d\tau\cdot(\cfrac{L_i^j - L_{i-1}^j}{d\chi})^n
    $$
    
    Parameters:
    -----------
    nu: function of 2 arguments
    
    initial: function of 1 argument
        lambda(chi) for tau = 0
        
    border: function of 1 argument
        lambda(tau), for chi = 0
        
    n : float
    
    dtau : float
    
    dchi : float
    
    ntau : int
        Number of tau grid lines
    
    nchi : int
        Number of chi grid lines
    
    Returns:
    --------
    sols: np.array shape (ntau, nchi)
        Solutions
    
    taus: np.array shape (ntau)
        tau values
        
    chis: np.array shape (nchis)
        chi values
    """
    taus = dtau*np.arange(ntau)
    chis = dchi*np.arange(nchi)
    
    sols = np.zeros([ntau, nchi])
    try:
        sols[0, :] = initial(chis)
    except TypeError:
        sols[0, :] = [initial(chi) for chi in chis]
    try:
        sols[:, 0] = border(taus)
    except TypeError:
        sols[:, 0] = [border(tau) for tau in taus]
    
    for j in range(ntau - 1):
        jnus = np.array([nu(j*dtau, chi) for chi in chis[1:]]) # should add a condition to nu taking array arguments
        sols[j+1, 1:] = sols[j, 1:] + dtau*jnus - dtau*((sols[j, 1:] - sols[j, :-1])/dchi)**n
            
    return sols, taus, chis




def approx1nonlinear(nu, initial, border, n=1.0, dtau=1e-3, dchi=1e-3, ntau=200, nchi=1000, 
                     method=sp.optimize.fsolve, use_fprime=True, bar=None):
    """
    Returns the solution of equation
    $$
    \cfrac{\partial\lambda}{\partial\tau} = \nu(\chi, \tau) - (\cfrac{\partial\lambda}{\partial\chi})^n
    $$
    using difference scheme:
    $$
    \cfrac{\lambda^{j+1}_{k+1} - \lambda^{j}_{k+1}}{s} = \nu^{j+1}_{k+1} - (
    \cfrac{\lambda^{j+1}_{k+1} - \lambda^{j+1}_{k}}{h})^n
    $$
    
    Parameters:
    -----------
    nu: function of 2 arguments
    
    initial: function of 1 argument
        lambda(chi) for tau = 0
        
    border: function of 1 argument
        lambda(tau), for chi = 0
        
    n : float
    
    dtau : float
    
    dchi : float
    
    ntau : int
        Number of tau grid lines
    
    nchi : int
        Number of chi grid lines
    
    method : function
        Method to solve nonlinear equation
        
    use_fprime : bool
        Use an analytical derivative as fprime parameter in method
        
    bar: tqdm bar or None:
        bar to update each iteration.
        Not draw bar if it's None
        
    Returns:
    --------
    sols: np.array shape (ntau, nchi)
        Solutions
    
    taus: np.array shape (ntau)
        tau values
        
    chis: np.array shape (nchis)
        chi values
    """
    taus = dtau*np.arange(ntau)
    chis = dchi*np.arange(nchi)
    
    sols = np.zeros([ntau, nchi])
    try:
        sols[0, :] = initial(chis)
    except TypeError:
        sols[0, :] = [initial(chi) for chi in chis]
    try:
        sols[:, 0] = border(taus)
    except TypeError:
        sols[:, 0] = [border(tau) for tau in taus]
    
    for j in range(ntau - 1):
        for k in range(nchi - 1):
            f = lambda sol: sol - sols[j, k+1] + dtau*((sol - sols[j+1, k])/dchi)**n - dtau*nu((j+1)*dtau, (k+1)*dchi)
            if use_fprime:
                fprime = lambda sol: 1 - (dtau/dchi)*((sol - sols[j+1, k])/dchi)**(n - 1)
                sols[j+1, k+1] = method(f, x0=sols[j, k], fprime=fprime)
            else:
                sols[j+1, k+1] = method(f, x0=sols[j, k]) 
            if bar is not None:
                bar.update()
    return sols, taus, chis




def approx2nonlinear(nu, initial, border, n=1.0, dtau=1e-3, dchi=1e-3, ntau=200, nchi=1000, 
                     method=sp.optimize.fsolve, use_fprime=True, bar=None):
    """
    Returns the solution of equation
    $$
    \cfrac{\partial\lambda}{\partial\tau} = \nu(\chi, \tau) - (\cfrac{\partial\lambda}{\partial\chi})^n
    $$
    using difference scheme:
    $$
    \cfrac{1}{2s}(\lambda_{k+1}^{j+1} + \lambda_{k}^{j+1} - \lambda_{k}^{j+1} - \lambda_k^j) = 
    \nu((j+0.5)s, (k+0.5)h) - 
    (\cfrac{1}{2h}(\lambda_{k+1}^{j+1} - \lambda_{k}^{j+1} + \lambda_{k}^{j+1} - \lambda_k^j))^n + 
    O(s + h^n)
    $$
    
    Parameters:
    -----------
    nu: function of 2 arguments
    
    initial: function of 1 argument
        lambda(chi) for tau = 0
        
    border: function of 1 argument
        lambda(tau), for chi = 0
        
    n : float
    
    dtau : float
    
    dchi : float
    
    ntau : int
        Number of tau grid lines
    
    nchi : int
        Number of chi grid lines
    
    method : function
        Method to solve nonlinear equation
        
    use_fprime : bool
        Use an analytical derivative as fprime parameter in method
        
    bar: tqdm bar or None:
        bar to update each iteration.
        Not draw bar if it's None
        
    Returns:
    --------
    sols: np.array shape (ntau, nchi)
        Solutions
    
    taus: np.array shape (ntau)
        tau values
        
    chis: np.array shape (nchis)
        chi values
    """
    taus = dtau*np.arange(ntau)
    chis = dchi*np.arange(nchi)
    
    sols = np.zeros([ntau, nchi])
    try:
        sols[0, :] = initial(chis)
    except TypeError:
        sols[0, :] = [initial(chi) for chi in chis]
    try:
        sols[:, 0] = border(taus)
    except TypeError:
        sols[:, 0] = [border(tau) for tau in taus]
    
    for j in range(ntau - 1):
        for k in range(nchi - 1):
            f = lambda sol: 0.5*(sol + sols[j+1, k] - sols[j, k+1] - sols[j, k])/dtau + (0.5*(sol - sols[j+1, k] + sols[j, k+1] - sols[j, k])/dchi)**n - nu((j + 0.5)*dtau, (k + 0.5)*dchi)
            if use_fprime:
                fprime = lambda sol: 0.5/dtau + (0.5/dchi)*(0.5*(sol - sols[j+1, k] + sols[j, k+1] - sols[j, k])/dchi)**(n - 1)
                fprime2 = lambda sol: (0.5/dchi)**2 * (0.5*(sol - sols[j+1, k] + sols[j, k+1] - sols[j, k])/dchi)**(n - 2)
                try:
                    sols[j+1, k+1] = method(f, x0=sols[j, k], fprime=fprime, fprime2=fprime2)
                except TypeError:
                    sols[j+1, k+1] = method(f, x0=sols[j, k], fprime=fprime)
            else:
                sols[j+1, k+1] = method(f, x0=sols[j, k]) 
            if bar is not None:
                bar.update()
    return sols, taus, chis