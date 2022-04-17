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

output = oi.OutputInterface('output_files/N2.out')
Ip = abs(output.saved_orbitals[output.HOMO][0])
kappa = np.sqrt(2*Ip)
omega = 0.057     # Frequency of light
Up = 0.22         # Ponderomotive potential
CEP = np.pi/2
N = 2

gto_coeffs = np.array(output.output_GTOs())
#gto_coeffs[:, -2] += 5.0

# %% CALCULATE THE PMD

dp = 4. / 150
px_lst = np.arange(-1, 1, dp)
pz_lst = np.arange(-1.5, 1.5, dp)
py = 0.
px_grd, pz_grd = np.meshgrid(px_lst, pz_lst)

# GTO
sfa_gto = sfa_lin.SFALinearPulse(Ip, Up, omega, N, CEP, 'GTO')
M_gto_grd = np.array(pool.starmap(sfa_gto.Mxz_List, zip(px_grd, repeat(py), pz_grd, repeat(gto_coeffs))))
M_gto_sqr = np.abs(np.flip(M_gto_grd, 0))**2
print('GTO done!')

# %% PLOT RESULT

plt.imshow(M_gto_sqr, interpolation='bicubic', cmap='Spectral_r',
                 norm=LogNorm(vmax=np.max(M_gto_sqr), vmin=np.max(M_gto_sqr)*1e-4),
                 extent=(-3, 3, -4, 4))
plt.xlabel(r'$p_\perp$ (a.u.)')
#ax2.set_ylabel(r'$p_\parallel$ (a.u.)')
plt.colorbar()

plt.show()

# %% EOF

'''
                   ----ENERGY BASED RESULTS----

 DIPOLE #        X                  Y                  Z   (A.U.)
 ###################################################################
        #  -1.6592090E-01      5.5949848E-01     -1.1915733E-01

  ALPHA #        X                  Y                  Z   (A.U.)
 ###################################################################
   X    #   3.5776779E+01     -4.5957440E-02     -6.8949239E+00
   Y    #  -4.5957440E-02      2.5523784E+01      8.2148174E-01
   Z    #  -6.8949239E+00      8.2148174E-01      4.4949015E+01

  BETA  #        X                  Y                  Z   (A.U.)
 ###################################################################
  XX    #   2.2239419E+00     -1.1161319E+01     -8.5547072E+00
  YY    #   7.0422175E+00     -2.2190306E+01     -3.7002792E+00
  ZZ    #   7.1217983E+00     -9.3673407E+00     -5.3854592E+01

  GAMMA #        XX                 YY                 ZZ   (A.U.)
 ###################################################################
  XX    #   5.2932592E+02      1.4915713E+02      3.6197889E+02
  YY    #   1.4915713E+02      3.4015102E+02      1.4370016E+02
  ZZ    #   3.6197889E+02      1.4370016E+02      1.2460077E+03




                   ----DIPOLE BASED RESULTS----

 DIPOLE #        X                  Y                  Z   (A.U.)
 ###################################################################
        #  -1.6592090E-01      5.5949847E-01     -1.1915733E-01

  ALPHA #        X                  Y                  Z   (A.U.)
 ###################################################################
   X    #   3.5776780E+01     -4.5960821E-02     -6.8949531E+00
   Y    #  -4.5934040E-02      2.5523810E+01      8.2147580E-01
   Z    #  -6.8949582E+00      8.2147680E-01      4.4949026E+01

  BETA  #        X                  Y                  Z   (A.U.)
 ###################################################################
  XX    #   2.2290237E+00     -1.1152321E+01     -8.5522619E+00
  YY    #   7.0482816E+00     -2.2181237E+01     -3.7000092E+00
  ZZ    #   7.1303242E+00     -9.3639087E+00     -5.3854892E+01

  GAMMA #        XX                 YY                 ZZ   (A.U.)
 ###################################################################
  XX    #   5.2465179E+02      9.6016599E+01      3.3324422E+02
  YY    #   9.6016599E+01      2.8137826E+02      8.0701856E+01
  ZZ    #   3.3324422E+02      8.0701856E+01      1.2123215E+03

'''