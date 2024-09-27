from bguriver.slope_patches import SlopePatches

n = 1.33
patch_starts = [0, 1, 3]
uplift_rates = [1, 2, 1]
slp = SlopePatches(patch_starts, uplift_rates, n)
print(f'SlopePatches: n={n:.2f}')
for i in range(slp.count()):
    print(f'tau_{i}={slp.patch_starts[i]:.2f} -> nu_{i}={slp.uplift_rates[i]:.2f}')
print()

tau, chi = 2.5, 2.5
lam = slp.get_elevation(tau, chi)
print(f'lam(tau={tau:.2f}, chi={chi:.2f}) = {lam:.2f}.\n')

tau = 7
rpl, rpr = slp.get_patches_relisation_borders(tau)
zpl, zpr = slp.get_stretch_zones_relisation_borders(tau)
print(f'Realisation borders at moment tau={tau:.2f}:')
for i in range(slp.count() - 1):
    print(f'Patch {i}:        [{rpl[i]:.2f}, {rpr[i]:.2f}]')
    print(f'Stretch Zone {i}: [{zpl[i]:.2f}, {zpr[i]:.2f}]')
i = slp.count() - 1
print(f'Patch {i}:        [{rpl[i]:.2f}, {rpr[i]:.2f}]\n')




from bguriver.approximation import approximate
n = slp.n
nu = slp.get_nu_value
initial = lambda chi: chi*slp.get_slopes()[0] # we sopouse, that nu(0) define the initial solution
border = lambda tau: tau*0

dtau = 0.1
ntau = 50
dchi = 0.1
nchi = 50

sol, tau, chi = approximate(nu, initial, border, n, 
                            dtau=dtau, ntau=ntau, 
                            dchi=dchi, nchi=nchi)
lam = slp.get_elevation(tau, chi)

print(f'sol.shape = {sol.shape}')
print(f'tau.shape = {tau.shape}')
print(f'chi.shape = {chi.shape}')
print(f'lam.shape = {lam.shape}')

mae = abs(lam - sol).mean()
print(f'MAE: {mae:.4f}')