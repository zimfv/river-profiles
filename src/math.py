import numpy as np


def lengths_check(borders, values, raise_error=True):
    """
    Returns True, if lengths are correct: borders array should contain 1 element less, than values array
    If this is not correct, it will raise error if raise_error or return Flase in other case.
    """
    r = len(values) - len(borders) == 1
    if raise_error and not r:
        msg = f'Wrong arrays lengths: borders array should contain 1 element less, than values array. But their length are {len(borders)} and {len(values)}.'
        raise ValueError(msg)
        
    borders = np.array(borders)
    if (borders[1:] - borders[:-1] <= 0).any():
        warnings.warn("The borders array is not strictly increasing")
    return r


def step_function(x, borders=[], values=[0]):
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


def step_integral(x0, x1, borders=[], values=[0], negative_backward=True):
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
        
    negative_backward: bool
        Set I(x0, x1) = -I(x1, x0) for cases, when x0 > x1 
        
    Returns:
    --------
    r : float or float array shape x.shape
    """
    lengths_check(borders, values, raise_error=True)
    
    x0 = np.array(x0)
    x1 = np.array(x1)*np.ones(x0.shape)
    x0 = np.array(x0)*np.ones(x1.shape)
    
    borders_ = np.concatenate([[-np.inf], borders, [+np.inf]])
    
    ones = np.ones(np.append(len(borders_), x1.shape).astype(int))
    res = ones*borders_.reshape(np.append(len(borders_), np.ones(len(x1.shape))).astype(int))
    res[res < x0] = (x0*ones)[res < x0]
    res[res > x1] = (x1*ones)[res > x1]
    res = res[1:] - res[:-1]
    res = res*np.array(values).reshape(np.append(len(values), np.ones(len(res.shape) - 1)).astype(int))
    res = res.sum(axis=0)
    
    if negative_backward and (x0 > x1).any():
        res = np.array(res)
        res[x0 > x1] = -step_integral(x0=x1[x0 > x1], x1=x0[x0 > x1], borders=borders, values=values)
        
    return res


def solve_bisect(f, x0, x1, xtol=1e-8, maxiter=200, iteration=0):
    """
    Bisection recursive solve method
    https://en.wikipedia.org/wiki/Bisection_method
    
    Parameters:
    -----------
    f : function, takes float or np.array as argument
    
    x0, x1: float or np.array same shape, as f-argument should be
    
    xtol: float
        The calculation will terminate if the relative error between two consecutive iterates is at most xtol.
    
    maxiter: int
        Maximal number of iterations
    
    iteration : int
        Iteration number
    
    Returns:
    --------
    float or np.array same shape, as x0 and x1
    """
    x0 = np.array(x0)
    x1 = np.array(x1)
    if (np.sign(f(x0))*np.sign(f(x1)) > 0).any():
        msg = 'f(x0) and f(x1) should have different signs.'
        raise ValueError(msg)
    xm = 0.5*(x0 + x1)
    change0 = np.sign(f(x0)) == np.sign(f(xm))
    change1 = np.sign(f(x1)) == np.sign(f(xm))
    x0[change0] = xm[change0]
    x1[change1] = xm[change1]
    if (abs(x0 - x1)[np.invert(np.isnan(x0 - x1))] < xtol).all() or (iteration == maxiter):
        return xm
    return solve_bisect(f, x0, x1, xtol=xtol, maxiter=maxiter, iteration=iteration + 1)