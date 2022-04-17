# %% IMPORTS

import numpy as np
import matplotlib.pyplot as plt
from extra_packages import OutputInterface as oi
import SFALinearPulse as sfa_lin
from matplotlib.colors import LogNorm
from itertools import repeat

# %% ENABLE MULTIPROCESSING

import multiprocessing

try:
    cpus = multiprocessing.cpu_count()
except NotImplementedError:
    cpus = 8   # arbitrary default

pool = multiprocessing.Pool(processes=cpus)

# %% GET GTO COEFFICIENTS

output = oi.OutputInterface('output_files/R-CHBrClF.out')
Ip = abs(output.saved_orbitals[output.HOMO][0])
kappa = np.sqrt(2*Ip)
omega = 0.057     # Frequency of light
Up = 0.22         # Ponderomotive potential
CEP = np.pi/2
N = 2

gto_coeffs = np.array(output.output_GTOs())

# %% CALCULATE THE PMD

dp = 4. / 150
px_lst = np.arange(-1, 1, dp)
pz_lst = np.arange(-1.5, 1.5, dp)
py = 0.
px_grd, pz_grd = np.meshgrid(px_lst, pz_lst)

# GTO
sfa_gto = sfa_lin.SFALinearPulse(Ip, Up, omega, N, CEP, 'GTO_MO_SPA')
M_gto_grd = np.array(pool.starmap(sfa_gto.Mxz_List, zip(px_grd, repeat(py), pz_grd, repeat(gto_coeffs))))
M_gto_sqr = np.abs(np.flip(M_gto_grd, 0))**2
print('GTO done!')

# %% PLOT RESULT

plt.imshow(M_gto_sqr, interpolation='bicubic', cmap='inferno',
                 norm=LogNorm(vmax=np.max(M_gto_sqr), vmin=np.max(M_gto_sqr)*1e-3),
                 extent=(-1, 1, -1.5, 1.5))
plt.xlabel(r'$p_\perp$ (a.u.)')
#ax2.set_ylabel(r'$p_\parallel$ (a.u.)')
plt.colorbar()

plt.show()

# %% EOF