from bguriver.slope_patches import SlopePatches
slp = SlopePatches([0, 1, 3], [3, 2, 1], 1)
print(slp.get_elevation(4, 4))

from bguriver.approximation import approximate
n = slp.n
nu = slp.get_nu_value
initial = lambda chi: chi*slp.get_slopes()[0] # we sopouse, that nu(0) define the initial solution
border = lambda tau: tau*0

dtau = 0.05
ntau = 61
dchi = 0.05
nchi = 61

from tqdm import tqdm
with tqdm(total=(ntau-1)*(nchi - 1)) as bar:
    # approximate the soulution
    sol, tau, chi = approximate(nu, initial, border, n, 
                                dtau=dtau, ntau=ntau, 
                                dchi=dchi, nchi=nchi, 
                                bar=bar)

print(f'sol.shape = {sol.shape}')
print(f'tau.shape = {tau.shape}')
print(f'chi.shape = {chi.shape}')