# distutils: language = c++
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Pulse version of the SFA for circular polarization
"""
Created on Tuesday Nov 24 10:29:00 2020

@author: asmaxwell

SFA pulse class
"""
#import scipy.optimize as op
import scipy.integrate as it
import scipy.special as sp
cimport scipy.special.cython_special as csp
import mpmath as mp
import functools
import multiprocessing
import scipy.integrate as it
from scipy.optimize import root, minimize
from scipy.interpolate import CubicSpline

import numpy as np
cimport numpy as np
from numpy cimport ndarray

cimport cython

from libcpp cimport bool

from libc.math cimport sin as sin_re
#from libc.math cimport sinc as sinc_re
from libc.math cimport cos as cos_re
#from libc.math cimport cot as cot_re
from libc.math cimport acos as acos_re
from libc.math cimport exp as exp_re
from libc.math cimport sqrt as sqrt_re
from libc.math cimport abs as abs_re
from libc.math cimport pow as pow_re
from libc.math cimport atan2 as atan2_re
from libc.math cimport tgamma as gamma

#from libc.math cimport  cyl_bessel_j as bessel_j_re

cdef extern from "sph_harm.h":
    double complex sph_harm(double complex x, double complex y, double complex z, double complex r, int l, int m)

#fused types
ctypedef fused dbl_or_cmplx:
    double
    double complex

ctypedef fused c_dbl_int:
    int
    double
    double complex

cdef extern from "<complex.h>" namespace "std" nogil:
    double complex exp(double complex z)
    double complex sin(double complex z)
    double complex cos(double complex z)
    #double complex cot(double complex z)
    double complex sqrt(double complex z)
    double complex acos(double complex z)
    double complex log(double complex z)
    double complex atan2(double complex z)
    #double complex pow(double complex z, c_dbl_int z)
    double real(double complex z)
    double imag(double complex z)


cdef double complex I1 = 1j
cdef double Pi = np.pi
cdef double rtPi = np.sqrt(Pi)
cdef double rt2 = np.sqrt(2.)
cdef int cacheSize = 2 ** 20

#Code to import old CQSFA prefactor code
# cdef extern from "HpgForms.h":
#     double complex calculateHMatrixElement(int target, double Ip, double complex pz, double px,  double complex ts, double complex Eft, double complex Aft, double theta)

### shorcut functions to efficntly switch trig between real and complex varients
cdef sin_c(dbl_or_cmplx t):
    if (dbl_or_cmplx is double):
        return sin_re(t)
    else:
        return sin(t)
cdef cos_c(dbl_or_cmplx t):
    if (dbl_or_cmplx is double):
        return cos_re(t)
    else:
        return cos(t)

cdef cot_c(dbl_or_cmplx t):
    return cos_c(t) / sin_c(t)

cdef class SFALinearPulse:
    '''
        Class to compute the transition amplitude M(p) and its dervative M_g(p) using the SFA and saddle point approximation
    '''
    #memeber variables like in C++!
    cdef readonly double Ip, Up, rtUp, omega, CEP, AlignmentAngle, kappa, Z
    cdef readonly int N, num_int
    cdef public int OAM
    cdef readonly str target
    cdef public list radial_splines
    #cdef object __weakref__ # enable weak referencing support

    def __init__(self, Ip_=0.5, Up_=0.44, omega_=0.057, N_=6, CEP_=0., target_="None", OAM_=1000, Z_=1):
        """
            Initialise field and target parameters defaults correspond to 800nm wl and 2 10^14 W/cm^2 intensity
            for the target 0=He, HeTheta=1, Ne=2, Ar=3, ArEx_4S=4, Xe=5, N2=6, N2Theta=7, O2Theta=8, H = 9 (and default case e.g. any other number)
            The ionization prefactor must be set independently
        """
        #Set pulse and targetvparameters
        self.Ip = Ip_
        self.Up = Up_
        self.rtUp = np.sqrt(Up_)  #must change this if Up is changed! Fixed by making Up readonly
        self.omega = omega_
        self.kappa = np.sqrt(2 * Ip_)
        self.Z = Z_
        self.OAM = OAM_  # Standard value of 1000 means that no OAM is not used

        self.N = N_
        self.CEP = CEP_

        self.AlignmentAngle = 0.
        self.target = target_
        self.test_target()
        self.radial_splines = []

    def test_target(self):
        """
        Give a warning if wrong target is selected
        """
        target_list = ['hyd1s_analytic', 'GTO', 'GTO_dress', 'asymp', 'asymp_martiny', 'dipole2', 'GTO_MO_SPA',
                       'dipole', 'dress_dip', 'cutoff', 'test', 'asymp_stark', 'single_clm']
        if not self.target == 'None' and self.target not in target_list:
            print('Warning: The chosen target is not known! Will calculate with prefractor set to 1.')

    #@functools.lru_cache(maxsize=cacheSize)
    cdef dbl_or_cmplx F(s, dbl_or_cmplx t):
        '''
        envelope for Sin^2 laser pulse
        '''
        #need to be fast can get evalutated millions of times
        if (real(t) < 0 or real(t) > 2 * s.N * Pi / s.omega):
            return 0
        return 2 * s.rtUp * sin_c(s.omega * t / (2 * s.N)) ** 2

    #@functools.lru_cache(maxsize=cacheSize)
    cpdef dbl_or_cmplx Af(s, dbl_or_cmplx t):
        '''
        Vector potential for Gaussian laser pulse
        '''
        if (real(t) < 0 or real(t) > 2 * s.N * Pi / s.omega):
            return 0
        cdef dbl_or_cmplx envelope, carrier
        envelope = s.F(t)
        if (dbl_or_cmplx is double):
            carrier = cos_re(s.omega * t + s.CEP)
        else:
            carrier = cos(s.omega * t + s.CEP)
        return envelope * carrier

    #@functools.lru_cache(maxsize=cacheSize)
    cpdef dbl_or_cmplx Ef(s, dbl_or_cmplx t):
        '''
        Electric Field for Sin^2 laser pulse
        '''
        if (real(t) < 0 or real(t) > 2 * s.N * Pi / s.omega):
            return 0
        cdef dbl_or_cmplx cos1, sin1, cos2, sin2, env
        env = s.omega * s.F(t)
        cos1 = cos_c(s.omega * t + s.CEP)
        sin1 = sin_c(s.omega * t + s.CEP)

        cos2 = cos_c(s.omega * t / (2 * s.N))
        sin2 = sin_c(s.omega * t / (2 * s.N))
        return env * sin1 - 2 * s.rtUp * s.omega * cos1 * cos2 * sin2 / s.N

    #@functools.lru_cache(maxsize=cacheSize)
    cpdef dbl_or_cmplx AfI(s, dbl_or_cmplx t):
        '''
            Integral of vector potential
        '''
        cdef double factor, a1, a2, a3, AI0 = 0
        cdef dbl_or_cmplx s1, s2, s3
        #Case if N==1 as general expression will diverge
        if (s.N == 1):
            A10 = (3 * s.rtUp * sin_re(s.CEP)) / (4. * s.omega)
            return -(s.rtUp * (2 * t * s.omega * cos_re(s.CEP) - 4 * sin_c(t * s.omega + s.CEP) + sin_c(
                2 * t * s.omega + s.CEP))) / (4. * rt2 * s.omega) - AI0
        else:
            AI0 = (s.rtUp * sin_re(s.CEP)) / (s.omega - s.N ** 2 * s.omega)  # AfI(0) to ensure limits are correct.
            a1 = s.N * s.N - 1
            a2 = (s.N + 1.) / s.N
            a3 = (s.N - 1.) / s.N
            factor = s.rtUp / (2 * s.omega * a1)

            s1 = sin_c(s.omega * t + s.CEP)
            s2 = sin_c(s.CEP + a2 * s.omega * t)
            s3 = sin_c(s.CEP + a3 * s.omega * t)

            return factor * (2 * a1 * s1 - s.N * s.N * (a3 * s2 + a2 * s3)) - AI0

    #@functools.lru_cache(maxsize=cacheSize)
    cpdef dbl_or_cmplx Af2I(s, dbl_or_cmplx t):
        '''
            Integral of vector potential squared
        '''
        if (s.N == 1):
            AI20 = (-25 * s.Up * sin_re(2 * s.CEP)) / (96. * s.omega)
            return (3 * s.Up * (4 * (6 * t * s.omega + t * s.omega * cos_re(2 * s.CEP) - 8 * sin_c(t * s.omega) + sin_c(
                2 * t * s.omega) + 3 * sin_c(2 * (s.CEP + t * s.omega)) - 4 * sin_c(2 * s.CEP + t * s.omega)) + sin_c(
                2 * (s.CEP + 2 * t * s.omega))) - 16 * s.Up0 * sin_c(2 * s.CEP + 3 * t * s.omega)) / (96. * s.omega)

        AI20 = (3 * s.Up * cos_re(s.CEP) * sin_re(s.CEP)) / (
                    4 * s.omega - 20 * s.N ** 2 * s.omega + 16 * s.N ** 4 * s.omega)
        cdef dbl_or_cmplx s1, s2, s3, s4, s5, s6, s7, c1
        cdef double a1 = s.omega * (s.N + 1) / (s.N), a2 = s.omega * (s.N - 1) / (s.N), a3 = s.omega * (2 * s.N + 1) / (
            s.N), a4 = s.omega * (2 * s.N - 1) / (s.N)

        s1 = sin_c(2 * s.omega * t)
        s2 = sin_c(s.omega * t / s.N)
        s3 = sin_c(2 * s.omega * t / s.N)
        s4 = sin_c(2 * s.CEP + a3 * t)
        s5 = sin_c(2 * s.CEP + 2 * a2 * t)
        s6 = sin_c(2 * s.CEP + 2 * a1 * t)
        s7 = sin_c(2 * s.CEP + a4 * t)
        c1 = cos_c(2 * s.omega * t)

        return (s.Up / 16) * (12 * t + 6 * c1 * sin_re(2 * s.CEP) / s.omega + 6 * s1 * cos_re(2 * s.CEP) / s.omega
                              - (16 * s.N / s.omega) * s2 + (
                                          2 * s.N / s.omega) * s3 - 8 * s4 / a3 + s5 / a2 + s6 / a1 - 8 * s7 / a4)

    #@functools.lru_cache(maxsize=cacheSize)
    cpdef dbl_or_cmplx S(s, double p, double theta, double phi, dbl_or_cmplx t):
        '''
        Action as given by the SFA for a Pulse
        '''
        cdef dbl_or_cmplx tTerms = (s.Ip + 0.5 * p * p) * t
        cdef dbl_or_cmplx linAI = p * cos_re(theta) * s.AfI(t)
        cdef dbl_or_cmplx quadAI = 0.5 * s.Af2I(t)
        return tTerms + linAI + quadAI

    cpdef dbl_or_cmplx DS(s, double p, double theta, double phi, dbl_or_cmplx t):
        cdef pz = p * cos_re(theta)
        return pz * s.Af(t) + 0.5 * s.Af(t) * s.Af(t) + 0.5 * p * p + s.Ip

    #@functools.lru_cache(maxsize=cacheSize)
    #REDO for linear field and with CEP!!!!!
    cpdef DSZ(s, double p, double theta, double phi):
        '''
            Derivative of the action tranformed by t->i N/omega Log[z] for esay solving
            This creates an 2(N+1) polynomial which can be efficeintly solved and the solutions easily 
            transformed back. This function passes the roots of the polynomial as numpy array so they can be solved using np.roots.          
            It is clear there will be N+1 solutions and their complex 
            conjugates.
        '''
        #Note for the linear case input varible phi is redundant, keep for genrality though

        cdef double complex exp_CEP = exp(I1 * s.CEP)

        #costruct polynomial in z of order 4*(N+1)
        poly_coeffs = np.zeros(4 * s.N + 4 + 1) + 0. * I1

        #0 order terms
        cdef c0 = s.Up * (exp_CEP * exp_CEP)
        poly_coeffs[0:5] = [c0 / 16, -c0 / 4, 3 * c0 / 8, -c0 / 4, c0 / 16]

        #N order terms (+= accounts for cases where coefficients combine)
        cdef c1 = p * s.rtUp * cos_re(theta) * exp_CEP
        poly_coeffs[s.N + 1:s.N + 4] += [-c1 / 2, c1, -c1 / 2]

        #2N order terms
        poly_coeffs[2 * s.N:2 * s.N + 5] += [s.Up / 8, -s.Up / 2, 2 * s.Ip + p * p + 3 * s.Up / 4, -s.Up / 2, s.Up / 8]

        #3N order terms
        cdef c3 = p * s.rtUp * cos_re(theta) / exp_CEP
        poly_coeffs[3 * s.N + 1:3 * s.N + 4] += [-c3 / 2, c3, -c3 / 2]

        #4N order terms
        cdef c4 = s.Up / (exp_CEP * exp_CEP)
        poly_coeffs[4 * s.N:] += [c4 / 16, -c4 / 4, 3 * c4 / 8, -c4 / 4, c4 / 16]

        return poly_coeffs

    cpdef double complex DSZ_val(s, double p, double theta, double phi, dbl_or_cmplx z):
        poly_coeffs = s.DSZ(p, theta, phi)
        cdef double complex sum_val = 0
        for n in range(0, len(poly_coeffs)):
            sum_val += poly_coeffs[n] * z ** n
        return sum_val

    cdef double complex addIfRealNeg(s, double complex ts):
        if (real(ts) < 0):
            return ts + 2 * Pi * s.N / s.omega
        else:
            return ts

    #@functools.lru_cache(maxsize=cacheSize)
    #Ensure correction transformation is done
    cpdef TimesGen(s, double p, double theta, double phi):
        '''
            Solution for times found by transforming the derivative of the action into
            a 4(N+1) polynomial and solving using np.roots. This should be a very effiecint way to get 
            all roots
        '''
        poly_coeffs = s.DSZ(p, theta, phi)
        z_roots = np.polynomial.polynomial.polyroots(poly_coeffs)
        #now we must transform back using t=I N Log[z]/omega
        ts_roots = I1 * s.N * np.log(z_roots) / s.omega
        ts_roots = [ts for ts in ts_roots if imag(ts) > 0]  #remove (divergent) solutions with negative imag
        #make sure all t values are in the domain [0, 2*pi*N/omega]
        ts_roots = [s.addIfRealNeg(ts) for ts in ts_roots]
        #sort real parts to easily select specific solutions
        return sorted(ts_roots, key=np.real)


    cpdef dbl_or_cmplx EI(self, dbl_or_cmplx t):
        """
        Integral of the electric field 
        """
        cdef double A0 = 2 * self.rtUp
        if self.N == 1:
            return A0 / 4 * (np.cos(2 * self.omega * t + self.CEP) - 2 * np.cos(self.omega * t + self.CEP))
        else:
            return A0 / 2 * np.cos(self.omega * t + self.CEP) * (np.cos(self.omega * t / self.N) - 1)


    cpdef dbl_or_cmplx EI2(self, dbl_or_cmplx t):
        """
        Integral of electric field squared 
        """
        cdef A0 = 2 * self.rtUp
        if self.N == 1:
            return A0**2 * self.omega * (t * self.omega / 4. - np.sin(self.omega * t) / 4. - np.sin(2. * self.omega * t + 2. * self.CEP)/16. + np.sin(3. * self.omega * t + 2. * self.CEP) / 12. - np.sin(4. * self.omega * t + 2. * self.CEP) / 32.)
        else:
            sin1 = np.sin(self.omega * (2. * self.N - 1.) * t/self.N + 2. * self.CEP)
            return A0**2 * self.omega * (-self.N * np.sin(self.omega*t/self.N) / 4. + self.omega * (3. * self.N**2 + 1.) * t
                            / (16. * self.N ** 2) + (self.N - 1.) * sin1 / (8. * (2. * self.N - 1.)) + (self.N + 1.) * sin1
                            / (8. * (2. * self.N + 1.)) + (self.N ** 2 - 1.) * np.sin(2. * self.omega * t/self.N)
                            / (32. * self.N) - (3. * self.N ** 2 - 1.) * np.sin(2. * self.omega * t + 2. * self.CEP)
                            / (32. * self.N ** 2) - (self.N ** 2 - 2.*self.N + 1.) * sin1
                            /(64. * self.N * (self.N - 1.)) - (self.N ** 2 + 2. * self.N + 1.) * sin1 / (64. * self.N * (self.N + 1.)))


    cpdef dbl_or_cmplx DEf(self, dbl_or_cmplx t):
        """
        Derivative of the electric field 
        """
        cdef double A0 = 2*self.rtUp
        cdef dbl_or_cmplx sin1 = np.sin(self.omega*t/(2*self.N))
        cdef dbl_or_cmplx cos1 = np.cos(self.omega * t / (2 * self.N))
        cdef dbl_or_cmplx cos2 = np.cos(self.omega*t + self.CEP)
        return A0*self.omega**2 * (-cos1**2 * cos2/(2*self.N**2) +
                                   2*sin1*np.sin(self.omega*t+self.CEP)*cos1/self.N +
                                   sin1**2*cos2/(2*self.N**2) + sin1**2 * cos2)


    cpdef double stark_minimize_action(self, double[:] t_vec, double p, double theta, double phi, pol_list):
        """
        The derivative of the stark action in a form ready to minimize. 
        pol_list = [mu_z, alpha_zz]
        """
        cdef double complex t = t_vec[0] + I1 * t_vec[1]
        return np.abs(self.DS(p, theta, phi, t) + pol_list[0] * self.Ef(t) + 0.5 * (self.Ef(t)) ** 2 * pol_list[1])


    cpdef double complex[:] stark_times_gen(self, double p, double theta, double phi, pol_list):
        """
        Saddle point times for the stark shifted action.
        """
        # First get the normal SP times
        cdef double complex[:] SP_times = np.array(self.TimesGen(p, theta, phi), dtype=complex)
        cdef double complex[:] res_list = np.array([], dtype=complex)
        cdef double max_t = 2*Pi*self.N / self.omega

        for time in SP_times:
            res = minimize(self.stark_minimize_action, np.array([time.real, time.imag]), args=(p, theta, phi, pol_list), bounds=((0, max_t), (0, max_t)))
            res_list = np.append(res_list, res.x[0] + 1.j * res.x[1])
        return np.array(res_list)


    cpdef dbl_or_cmplx S_stark(self, double p, double theta, double phi, dbl_or_cmplx ts, pol_list):
        """
        Action modified with the Stark shift. 
        pol_list is the dipole moment and polarizability [mu_z, alpha_zz]
        """
        return pol_list[0] * self.EI(ts) + pol_list[1] / 2 * self.EI2(ts) + self.S(p, theta, phi, ts)


    cpdef double complex DDS_stark(self, double p, double theta, double phi, double complex ts, pol_list):
        """
        Double derivative of action with Stark modification 
        """
        return self.DDS(p, theta, phi, ts) + pol_list[0] * self.DEf(ts) + pol_list[1] * self.Ef(ts) * self.DEf(ts)


    cdef double complex mo_spe(self, double complex t, double p, double theta, double phi, double Rz):
        # px = p * sin_re(theta) * cos_re(phi) # py = p * sin_re(theta) * sin_re(phi) # pz = p * cos_re(theta)
        cdef double complex E_dot_R = self.Ef(t) * Rz
        cdef pz = p * cos_re(theta)
        return pz * self.Af(t) + 0.5 * self.Af(t) * self.Af(t) + 0.5 * p * p + self.Ip + E_dot_R

    cpdef mo_spe_real(self, Y, double p, double theta, double phi, double Rz):
        cdef double t_re, t_im
        t_re, t_im = Y
        cdef double complex ts = t_re + 1j * t_im
        val = self.mo_spe(ts, p, theta, phi, Rz)
        return np.array([val.real, val.imag])

    cpdef mo_times_gen(self, double p, double theta, double phi, double Rz):
        guess_times = self.TimesGen(p, theta, phi)
        times = []
        for ts in guess_times:
            sol = root(self.mo_spe_real, np.array([ts.real, ts.imag]), args=(p, theta, phi, Rz))
            times.append(complex(*sol.x))
            if not sol.success:
                print('Root finder failed to converge!')
        return np.array(times, dtype=complex)

    cpdef double complex DDEf(self, double complex t):
        '''
        Double time derivative of the electric field
        '''
        cdef double prefactor = - 2.0 * sqrt_re(self.Up) * self.omega**2 / self.N**2
        cdef double complex term1 = cos(self.omega * t + self.CEP) * (self.N**2 + 1.0) \
                                    * cos(self.omega * t / (2.0 * self.N))**2
        cdef double complex term2 = -2.0 * self.N * sin(self.omega * t / (2.0 * self.N)) \
                                    * cos(self.omega * t / (2.0 * self.N)) * sin(self.omega * t + self.CEP)
        cdef double complex term3 = -cos(self.omega * t + self.CEP) * (self.N**2 + 0.5)
        return prefactor * (term1 + term2 + term3)

    cpdef double complex DDPhi(self, double p, double theta, double phi, double complex t, double Rz):
        '''
        Double time derivative of the modified action for the MO-SPA
        '''
        cdef pz = p * cos_re(theta)
        cdef double complex term1 = -(pz + self.Af(t)) * self.Ef(t)
        cdef double complex term2 = Rz * self.DDEf(t)
        return term1 + term2


    #1 varible determinant for saddle point approximation
    cpdef double complex DDS(s, double p, double theta, double phi, double complex t):
        '''Second order derivative of action wrt t'''
        return -(p * cos_re(theta) + s.Af(t)) * s.Ef(t)

    #### PREFACTORS and support functions ####
    cpdef double complex sph_harm_OAM(self, double complex x, double complex y, double complex z,
                                      double complex p, int l, int m, double phi):
        """
        Spherical harmonic function in Cartesian coordinates that also allows for OAM selection
        """
        if self.OAM == 1000:  # No OAM selection - return the whole thing!
            return sph_harm(x, y, z, p, l, m)
        else:  # OAM selection - we must kill the exponential from the spherical harmonic
            return sph_harm(x, y, z, p, l, m) / np.exp(1j * m * phi)


    cpdef double I2_factor(self, int l):
        """
        Factor needed to calculate the transition amplitude for an asymptotic WF in the LG 
        """
        cdef double nu = self.Z / self.kappa
        cdef double factor = self.kappa ** (nu + 2.) * gamma(l + nu + 3.) / gamma(l + 3./2.) * \
                             (self.kappa**(l+1.) * 2.**(l+0.5))**(-1) * 2**(-nu-2.)

        # Determine the hypergeometric function in the saddle points (z=1 always):
        cdef double a = 0.5 * (l - nu - 1.)
        cdef double b = 0.5 * (l - nu)
        cdef double c = l + 3. / 2.
        cdef double hyper_geo = gamma(c) * gamma(c - a - b) / (gamma(c - a) * gamma(c - b))

        return hyper_geo * factor


    cpdef double complex d_asymp_Er(self, double p, double theta, double phi, double complex ts, clm_array):
        """
        Prefactor for E*r for the asymptotic wave function
        clm_array is the expansion coeffs. for the asymptotic wave function. 
        """
        cdef double I2_p, I2_m, alpha_p, alpha_m, sn, theta_t, px, py, pz
        cdef double complex factor1, factor2, pz_t, p_t

        # Values needed
        cdef double complex d_res = 0.
        cdef double complex ddS = self.DDS(p, theta, 0., ts)
        cdef double nu = self.Z / self.kappa
        cdef int max_l = clm_array.shape[1]
        cdef double complex E = self.Ef(ts)

        # Find the coordinates
        px = p * sin_re(theta) * cos_re(phi)
        py = p * sin_re(theta) * sin_re(phi)
        pz = p * cos_re(theta)

        sn = 1. if self.Af(ts).imag > 0 else -1.
        pz_t = 1j * sn * sqrt_re(2 * self.Ip + px ** 2 + py ** 2)  # pz+s.Af(ts) #tilde{pz}
        p_t = 1j * sqrt_re(2 * self.Ip)  # sqrt(px**2+py**2+pz_2**2) # =modulus{tilde{p}}

        # Find terms only dependent on ts
        factor1 = 1.j**(nu/2.+1) * E * self.kappa**nu * gamma(nu/2.+1)/(2*np.sqrt(2*np.pi)) \
                                     * ddS**(-nu-2.) * (2*ddS)**(nu/2.+1)

        # Now loop over the contributions from the Clm's
        for l in range(0, max_l):
            # Find stuff not dependent on m:
            factor2 = (p_t/(1.j * self.kappa))**l

            for m in range(-l, l + 1):
                # Get clm
                sign = 0 if m >= 0 else 1
                clm = clm_array[sign, l, abs(m)]

                if np.abs(clm) == 0:  # Don't calculate if it's all 0...
                    continue
                if self.OAM != 1000:  # Possibility for OAM selection - only takes the ones with matching m
                    if m != self.OAM:
                        continue

                # Factors from recursion of spherical harmonics:
                alpha_p = np.sqrt((l - m + 1.) * (l + m + 1.) / ((2. * l + 1.) * (2. * l + 3.)))
                alpha_m = np.sqrt((l - m) * (l + m) / ((2. * l - 1.) * (2. * l + 1.)))

                # Now add to the result
                d_res += clm * factor2 * (alpha_m * self.sph_harm_OAM(px, py, pz_t, p_t, l-1, m, phi) * self.kappa**2/p_t
                                          - alpha_p * self.sph_harm_OAM(px, py, pz_t, p_t, l+1, m, phi)*p_t)

        if self.OAM == 1000:  # OAM selection is not activated
            return d_res * factor1
        else:  # OAM selection is activated - remember the i**OAM!
            return d_res * 1j ** self.OAM * factor1


    cpdef double complex d_asymp_martiny(self, double p, double theta, double phi, double complex ts, clm_array):
        """
        This is the prefactor calculated in Martiny's phd (eq. 4.29), simplified using an identity for the 
        gamma functions.
        """
        cdef double sn, theta_t, px, py, pz
        cdef double complex factor1, factor2, pz_t, p_t
        cdef double complex d_res = 0.

        # Values needed
        cdef double complex ddS = self.DDS(p, theta, 0., ts)
        cdef double nu = self.Z / self.kappa
        cdef int max_l = clm_array.shape[1]

        # Find the coordinates
        px = p * sin_re(theta) * cos_re(phi)
        py = p * sin_re(theta) * sin_re(phi)
        pz = p * cos_re(theta)
        sn = 1. if self.Af(ts).imag > 0 else -1.
        pz_t = 1j * sn * sqrt_re(2 * self.Ip + px ** 2 + py ** 2)  # pz+s.Af(ts) : tilde{pz}
        p_t = 1j * sqrt_re(2 * self.Ip)  # sqrt(px**2+py**2+pz_2**2) : modulus{tilde{p}}

        #cos_theta_t = np.imag(pz_t) / np.imag(p_t)  # pz_t and p_t are both imaginary in saddle points

        # Start calculating the prefactor itself! Everything not depending on l,m:
        factor1 = gamma(nu / 2. + 1.) / sqrt_re(2 * np.pi) * (1j * self.kappa) ** nu * (2. / 1j * ddS) ** (
                    nu / 2.) / ddS ** nu

        # Everything depending on l,m:
        for l in range(0, max_l):
            factor2 = (p_t / (1.j * self.kappa)) ** l

            for m in range(-l, l + 1):
                # Get clm
                sign = 0 if m >= 0 else 1
                clm = clm_array[sign, l, abs(m)]

                if abs(clm) == 0:  # Don't calculate if it's zero anyway...
                    continue
                if self.OAM != 1000:  # Possibility for OAM selection - only calculate m values matching OAM
                    if m != self.OAM:
                        continue

                d_res += factor2 * clm * self.sph_harm_OAM(px, py, pz_t, p_t, l, m, phi)

        if self.OAM == 1000:  # OAM selection is not activated
            return d_res * factor1 * 1j
        else:  # OAM selection is activated - remember the i**OAM!
            return d_res * factor1 * 1j * 1j ** self.OAM

    cpdef double complex d_asymp_lm(self, double p, double theta, double phi, double complex ts, lm_array):
        """
        E*r prefactor for a single l and m, for use in analytical orientation averaging.
        Equivalent to the Martiny prefactor, but with
        """
        l, m = lm_array
        cdef double sn, theta_t, px, py, pz
        cdef double complex factor1, factor2, pz_t, p_t
        cdef double complex d_res = 0.

        # Values needed
        cdef double complex ddS = self.DDS(p, theta, 0., ts)
        cdef double nu = self.Z / self.kappa

        # Find the coordinates
        px = p * sin_re(theta) * cos_re(phi)
        py = p * sin_re(theta) * sin_re(phi)
        pz = p * cos_re(theta)
        sn = 1. if self.Af(ts).imag > 0 else -1.
        pz_t = 1j * sn * sqrt_re(2 * self.Ip + px ** 2 + py ** 2)  # pz+s.Af(ts) : tilde{pz}
        p_t = 1j * sqrt_re(2 * self.Ip)  # sqrt(px**2+py**2+pz_2**2) : modulus{tilde{p}}

        # Start calculating the prefactor itself! Everything not depending on l,m:
        factor1 = gamma(nu / 2. + 1.) / sqrt_re(2 * np.pi) * (1j * self.kappa) ** nu * (2. / 1j * ddS) ** (nu / 2.) \
                  / ddS ** nu

        factor2 = (p_t / (1.j * self.kappa)) ** l

        d_res += factor2 * self.sph_harm_OAM(px, py, pz_t, p_t, l, m, phi)

        return d_res * factor1 * 1j


    cpdef double complex d_dipole(self, double p, double theta, double phi, double complex ts, alpha_list, clm_array):
        """
        Prefactor for dipole matrix term...
        """
        cdef double sn, px, py, pz, alpha_p, alpha_m, beta
        cdef double complex pz_t, p_t, factor1, res_prime
        cdef double complex d_res = 0.

        # Values needed
        cdef double complex ddS = self.DDS(p, theta, 0., ts)
        cdef double nu = self.Z / self.kappa
        cdef int max_l = clm_array.shape[1]
        cdef double complex E_field = self.Ef(ts)

        # Find dipole through polarizability
        cdef double complex mu_x = alpha_list[0] * E_field
        cdef double complex mu_y = alpha_list[1] * E_field
        cdef double complex mu_z = alpha_list[2] * E_field

        # Find the coordinates
        px = p * sin_re(theta) * cos_re(phi)
        py = p * sin_re(theta) * sin_re(phi)
        pz = p * cos_re(theta)
        sn = 1. if self.Af(ts).imag > 0 else -1.
        pz_t = 1j * sn * sqrt_re(2 * self.Ip + px ** 2 + py ** 2)  # pz+s.Af(ts) : tilde{pz}
        p_t = 1j * sqrt_re(2 * self.Ip)  # sqrt(px**2+py**2+pz_2**2) : modulus{tilde{p}}

        # Calculate everything not dependend on l and m:
        factor1 = self.kappa**(nu-2.) * 2.**(-3./2.) * 1.j**(nu-1.) * gamma((nu-1.)/2.) / np.sqrt(np.pi) \
                  * (2. / 1.j * ddS) ** ((nu - 1.) / 2.) / (ddS ** (nu - 1.))

        # Loops for the rest (l,m and l',m')
        for l in range(0, max_l):
            for m in range(-l, l + 1):
                # Get clm
                sign = 0 if m >= 0 else 1
                clm = clm_array[sign, l, abs(m)]
                if abs(clm) == 0:  # Don't calculate if it's zero anyway...
                    continue

                # Find l' and m' terms
                res_prime = 0.
                for li in [-1, 1]:
                    if li + l < 0:  # These are 0
                        continue
                    for mi in [-1, 0, 1]:
                        if abs(mi + m) > l:  # These are 0
                            continue

                        if self.OAM != 1000:  # Possibility for OAM selection - only calculate m values matching OAM
                            if m + mi != self.OAM:
                                continue

                        alpha_p = 0.  # Reset values
                        alpha_m = 0.
                        beta = 0.
                        if mi == 1:  # Find the alpha+ terms
                            alpha_p = -1. * (1. if li == 1 else 0.) * np.sqrt((l+m+1.)*(l+m+2.)/((2.*l+1.)*(2.*l+3.))) \
                                      + (1. if li == -1 else 0.) * np.sqrt((l-m)*(l-m-1.)/((2.*l-1.)*(2.*l+1.)))
                        elif mi == 0:  # Find the beta term
                            beta = (1. if li == 1 else 0.) * np.sqrt((l-m+1.)*(l+m+1.)/((2.*l+1.)*(2.*l+3.))) \
                                   + (1. if li == -1 else 0.) * np.sqrt((l-m)*(l+m)/((2.*l-1.)*(2.*l+1.)))
                        else:  # Find the alpha- terms
                            alpha_m = 1. * (1. if li == 1 else 0.) * np.sqrt((l-m+1.)*(l-m+2.)/((2.*l+1.)*(2.*l+3.))) \
                                      - (1. if li == -1 else 0.) * np.sqrt((l+m)*(l+m-1.)/((2.*l-1.)*(2.*l+1.)))

                        # Now calculate the contribution from the specific l', m'
                        res_prime += (-1.j)**(l+li) * self.sph_harm_OAM(px, py, pz_t, p_t, l+li, m+mi, phi) * p_t**(l+li) / self.kappa**(l+li) \
                                    * (mu_x/2.*(alpha_p + alpha_m) + mu_y/(2.*1.j)*(alpha_p - alpha_m) + mu_z*beta)

                d_res += clm * res_prime

        if self.OAM == 1000:  # OAM selection is not activated
            return 1.j * factor1 * d_res
        else:  # OAM selection is activated - remember the i**OAM!
            return 1.j * factor1 * d_res * 1j ** self.OAM



    cpdef dbl_or_cmplx hermite_poly(self, int n, dbl_or_cmplx z):
        if n == 0:
            return 1.
        elif n == 1:
            return 2. * z
        elif n == 2:
            return 4. * z ** 2 - 2.
        elif n == 3:
            return 8 * z ** 3 - 12. * z
        elif n == 4:
            return 16. * z ** 4 - 48. * z ** 2 + 12.
        elif n == 5:
            return 32. * z ** 5 - 160. * z ** 3 + 120. * z
        else:
            print('Too high n: not implemented')
            return 0.0


    cpdef double complex di_gto(self, double p, double theta, double phi, dbl_or_cmplx ts,
                                double front_factor, double alpha, int i, int j, int k,
                                double xa, double ya, double za):
        """ Bound state prefactor d(p,t)=<p+A|H|0> in the LG using LCAO and GTOs from GAMESS coefficients """
        cdef double px = p * sin_re(theta) * cos_re(phi)
        cdef double py = p * sin_re(theta) * sin_re(phi)
        cdef double complex pz = p * cos_re(theta) + self.Af(ts)

        cdef double complex result = 2. ** (-(i + j + k) - 3. / 2) * alpha ** (-(i + j + k) / 2. - 3. / 2) \
                                     * np.exp(-1j * (px * xa + py * ya + pz * za)) \
                                     * np.exp(-(px ** 2 + py ** 2 + pz ** 2) / (4. * alpha)) \
                                     * np.exp(-1.j * np.pi * (i + j + k) / 2.) * self.Ef(ts) \
                                     * self.hermite_poly(i, px / (2. * sqrt_re(alpha))) \
                                     * self.hermite_poly(j, py / (2. * sqrt_re(alpha))) \
                                     * (za * self.hermite_poly(k, pz / (2. * sqrt_re(alpha))) -
                                        1.j / (2. * sqrt_re(alpha)) * self.hermite_poly(k + 1,
                                                                                        pz / (2. * sqrt_re(alpha))))
        return front_factor * result


    cpdef double complex d_gto(self, double p, double theta, double phi, dbl_or_cmplx ts, coefficients):
        """ Total GTO prefactor from GAMESS coefficients """
        cdef double complex result = 0.
        cdef double x_a, y_a, z_a, alpha, front_factor
        cdef int i, j, k

        for row in coefficients:
            front_factor, alpha, i, j, k, x_a, y_a, z_a = row
            result += self.di_gto(p, theta, phi, ts, front_factor, alpha, i, j, k, x_a, y_a, z_a)
        return result


    cpdef double complex di_gto_dress(self, double p, double theta, double phi, dbl_or_cmplx ts,
                                      double front_factor, double alpha, int i, int j, int k,
                                      double xa, double ya, double za):
        """ Dressed bound state prefactor d(p,t)=<p+A|H|0> in the LG using LCAO and GTOs from GAMESS coefficients """
        cdef double px = p * sin_re(theta) * cos_re(phi)
        cdef double py = p * sin_re(theta) * sin_re(phi)
        cdef double complex pz = p * cos_re(theta) + self.Af(ts)

        cdef double complex result = 2. ** (-(i + j + k) - 3. / 2) * alpha ** (-(i + j + k) / 2. - 3. / 2) \
                                     * np.exp(-1j * (px * xa + py * ya + p * cos_re(theta) * za)) \
                                     * np.exp(-(px ** 2 + py ** 2 + pz ** 2) / (4. * alpha)) \
                                     * np.exp(-1.j * np.pi * (i + j + k) / 2.) * self.Ef(ts) \
                                     * self.hermite_poly(i, px / (2. * sqrt_re(alpha))) \
                                     * self.hermite_poly(j, py / (2. * sqrt_re(alpha))) \
                                     * (-1.j / (2. * sqrt_re(alpha)) * self.hermite_poly(k + 1,
                                                                                         pz / (2. * sqrt_re(alpha))))
        return front_factor * result


    cpdef double complex d_gto_dress(self, double p, double theta, double phi, dbl_or_cmplx ts, coefficients):
        """ Total dressed GTO prefactor from GAMESS coefficients """
        cdef double complex result = 0.
        cdef double x_a, y_a, z_a, alpha, front_factor
        cdef int i, j, k

        for row in coefficients:
            front_factor, alpha, i, j, k, x_a, y_a, z_a = row
            result += self.di_gto_dress(p, theta, phi, ts, front_factor, alpha, i, j, k, x_a, y_a, z_a)
        return result

    cdef double complex legendre_poly(self, int l, int m, double complex x):
        if abs(m) > l:
            return 0.0
        cdef double complex pref = gamma(l + m + 1)/gamma(l - m + 1) * (1.0 - x*x)**(m/2.0) / (2.0**m * gamma(m + 1))
        return pref * sp.hyp2f1(m - l, m + l + 1, m + 1, (1.0 - x) / 2.0)


    cdef sphericalY(self, int l, int m, double complex cos_theta, double phi):
        cdef double complex pref = sqrt_re((2.0 * l + 1.0) / 2.0 * gamma(l - m) / gamma(l + m))
        return pref * self.legendre_poly(l, m, cos_theta) * exp(-1j * m * phi)

    cdef double complex radial_int(self, double p, double theta, double phi, int l, int m, double complex ts):
        if abs(m) > l:
            return 0.0

        cdef double nu = 1.0/self.kappa
        return 1.0j**(nu - 1.0) * p**l * gamma((l + nu)/2.0) / (self.kappa**(l + 3.0*nu + 2.0) * (l + nu)) \
               * (-2.0j*self.DDS(p, theta, phi, ts))**((nu - 1.0)/2) / self.DDS(p, theta, phi, ts)**(nu - 1.0)
        #return 1.0j**(nu - 1.0)/4.0 * self.kappa**(nu - l - 2.0) * gamma((nu - 1.0)/2.0) * p**l \
        #       * (2.0j*self.DDS(p, theta, phi, ts))**((nu - 1.0)/2.0) / (self.DDS(p, theta, phi, ts))**(nu - 1.0)


    cpdef double complex di_dip(self, double p, double theta, double phi, double complex ts,
                               double complex c_lm, int l, int m):
        """ Asymptotic dipole interaction prefactor """
        cdef double complex result = 0.0 + 0.0j
        cdef double px = p * sin_re(theta) * cos_re(phi)
        cdef double py = p * sin_re(theta) * sin_re(phi)
        cdef double pz = p * cos_re(theta)
        cdef double sn = 1. if self.Af(ts).imag > 0 else -1.
        cdef double complex pz_t = 1j * sn * sqrt_re(2 * self.Ip + px ** 2 + py ** 2)  # pz+s.Af(ts) : tilde{pz}
        cdef double complex p_t = 1j * sqrt_re(2 * self.Ip)  # sqrt(px**2+py**2+pz_2**2) : modulus{tilde{p}}
        cdef double cos_theta_t = np.imag(pz_t) / np.imag(p_t)  # pz_t and p_t are both imaginary in saddle points

        mu = np.array([0.0, 0.0, 0.0]) # np.array([-1.6592090E-01, 5.5949848E-01, -1.1915733E-01])
        alpha = np.array([[3.5776779E+01, -4.5957440E-02, -6.8949239E+00],
                          [-4.5957440E-02, 2.5523784E+01, 8.2148174E-01],
                          [-6.8949239E+00, 8.2148174E-01, 4.4949015E+01]])

        cdef double complex mu_norm = np.sqrt(mu[0]**2 + mu[1]**2 + (mu[2] + np.sum(alpha[:,2])*self.Ef(ts))**2)

        cdef double complex x_term_1 = -1.0j**(l + 1.0) * (-0.5) \
                                       * sqrt_re((l + m + 1.0)*(l + m + 2.0)/((2.0*l + 1.0)*(2.0*l + 3.0))) \
                                       * self.sph_harm(m + 1, l + 1, cos_theta_t, phi) \
                                       * self.radial_int(p, theta, phi, l + 1, m + 1, ts)

        cdef double complex x_term_2 = -1.0j**(l - 1.0) * 0.5 \
                                       * sqrt_re((l - m)*(l - m - 1.0)/((2.0*l - 1.0)*(2.0*l + 1.0))) \
                                       * self.sph_harm(m + 1, l - 1, cos_theta_t, phi) \
                                       * self.radial_int(p, theta, phi, l - 1, m + 1, ts)

        cdef double complex x_term_3 = -1.0j**(l + 1.0) * 0.5 \
                                       * sqrt_re((l - m + 1.0)*(l - m + 2.0)/((2.0*l + 1.0)*(2.0*l + 3.0))) \
                                       * self.sph_harm(m - 1, l + 1, cos_theta_t, phi) \
                                       * self.radial_int(p, theta, phi, l + 1, m - 1, ts)

        cdef double complex x_term_4 = -1.0j**(l - 1.0) * (-0.5) \
                                       * sqrt_re((l + m)*(l + m - 1)/((2.0*l + 1.0)*(2.0*l - 1.0))) \
                                       * self.sph_harm(m - 1, l - 1, cos_theta_t, phi) \
                                       * self.radial_int(p, theta, phi, l - 1, m - 1, ts)

        cdef double complex y_term_1 = -1.0j**(l + 1.0) * (-1.0/2.0j) \
                                       * sqrt_re((l + m + 1.0)*(l + m + 2.0)/((2.0*l + 1)*(2.0*l + 3.0))) \
                                       * self.sph_harm(m + 1, l + 1, cos_theta_t, phi) \
                                       * self.radial_int(p, theta, phi, l + 1, m + 1, ts)
        cdef double complex y_term_2 = -1.0j**(l - 1.0) * 1.0/2.0j \
                                       * sqrt_re((l - m)*(l - m - 1.0)/((2.0*l - 1.0)*(2.0*l + 1.0))) \
                                       * self.sph_harm(m + 1, l - 1, cos_theta_t, phi) \
                                       * self.radial_int(p, theta, phi, l - 1, m + 1, ts)

        cdef double complex y_term_3 = -1.0j**(l + 1) * (-1.0/2.0j) \
                                       * sqrt_re((l - m + 1.0)*(l - m + 2.0)/((2.0*l + 1.0)*(2.0*l + 3.0))) \
                                       * self.sph_harm(m - 1, l + 1, cos_theta_t, phi) \
                                       * self.radial_int(p, theta, phi, l + 1, m - 1, ts)

        cdef double complex y_term_4 = -1.0j**(l - 1.0) * 1.0/2.0j \
                                       * sqrt_re((l + m)*(l + m - 1.0)/((2.0*l - 1.0)*(2.0*l + 1.0))) \
                                       * self.sph_harm(m - 1, l - 1, cos_theta_t, phi) \
                                       * self.radial_int(p, theta, phi, l - 1, m - 1, ts)

        cdef double complex z_term_1 = -1.0j**(l + 1.0) \
                                       * sqrt_re((l - m + 1.0)*(l + m + 1.0)/((2.0*l + 1.0)*(2.0*l + 3.0))) \
                                       * self.sph_harm(m, l + 1, cos_theta_t, phi) \
                                       * self.radial_int(p, theta, phi, l + 1, m, ts)

        cdef double complex z_term_2 = -1.0j**(l - 1.0) \
                                       * sqrt_re((l - m)*(l + m)/((2.0*l - 1.0)*(2.0*l + 1.0))) \
                                       * self.sph_harm(m, l - 1, cos_theta_t, phi) \
                                       * self.radial_int(p, theta, phi, l - 1, m, ts)

        result += -sqrt_re(np.pi/2.0) * c_lm * mu_norm * (x_term_1 + x_term_2 + x_term_3 + x_term_4 + y_term_1
                                                          + y_term_2 + y_term_3 + y_term_4 + z_term_1 + z_term_2)
        terms = [x_term_1, x_term_2, x_term_3, x_term_4, y_term_1, y_term_2, y_term_3, y_term_4, z_term_1, z_term_2]
        for i, term in enumerate(terms):
            if np.isnan(term):
                print(f'Term {i} is NaN: {term}')
        return result


    cpdef double complex d_dip(self, double p, double theta, double phi, double complex ts, coeffs):
        cdef double complex result = 0.0 + 0.0j
        cdef l_max = coeffs.shape[1]
        cdef double complex clm
        cdef int sgn
        for l in range(0, l_max):
            for m in range(-l, l + 1):
                sgn = 0 if m >= 0 else 1
                clm = coeffs[sgn, l, abs(m)]
                if abs(clm) == 0:
                    continue
                result += self.di_dip(p, theta, phi, ts, clm, l, m)
        if np.isnan(result) or np.isinf(result):
            print(f'Result is {np.inf if np.isinf(result) else np.nan}')
        return result


    def load_splines(self, r_cutoff, max_p=10, max_l=10):
        """
        Test for loading splines!
        """
        nu = self.Z / self.kappa
        r_list = np.linspace(0, r_cutoff, 1000)
        p_list = np.linspace(0, max_p, 200)

        integrand = lambda r, p, l: r ** (nu - 1) * np.exp(-self.kappa * r) * sp.spherical_jn(l, p * r)
        integral = lambda p, l: np.trapz(integrand(r_list, p, l), r_list)

        for l in range(max_l + 1):
            int_list = [integral(p, l) for p in p_list]
            self.radial_splines.append(CubicSpline(p_list, int_list))


    cpdef double complex d_test_cutoff(self, double p, double theta, double phi, double complex t, clm_array, alpha_list):
        """
        Test for the radial cutoff integral 
        """
        cdef double sn, px, py, pz, alpha_p, alpha_m, beta
        cdef double complex pz_t, p_t, res_prime
        cdef double complex d_res = 0.

        # Values needed
        cdef double nu = self.Z / self.kappa
        cdef int max_l = clm_array.shape[1]
        cdef double complex E_field = self.Ef(t)  # 2*np.sqrt(self.Up)
        cdef double complex r_integral

        # Find dipole through polarizability
        cdef double complex mu_x = alpha_list[0] * E_field
        cdef double complex mu_y = alpha_list[1] * E_field
        cdef double complex mu_z = alpha_list[2] * E_field

        # Find the coordinates
        px = p * sin_re(theta) * cos_re(phi)
        py = p * sin_re(theta) * sin_re(phi)
        pz = p * cos_re(theta)
        sn = 1. if self.Af(t).imag > 0 else -1.
        pz_t = 1j * sn * sqrt_re(2 * self.Ip + px ** 2 + py ** 2)  # pz+s.Af(ts) : tilde{pz}
        p_t = 1j * sqrt_re(2 * self.Ip)  # sqrt(px**2+py**2+pz_2**2) : modulus{tilde{p}}

        # Loops for the rest (l,m and l',m')
        for l in range(0, max_l):
            for m in range(-l, l + 1):
                # Get clm
                sign = 0 if m >= 0 else 1
                clm = clm_array[sign, l, abs(m)]
                if abs(clm) == 0:  # Don't calculate if it's zero anyway...
                    continue

                # Find l' and m' terms
                res_prime = 0.
                for li in [-1, 1]:
                    if li + l < 0:  # These are 0
                        continue
                    for mi in [-1, 0, 1]:
                        if abs(mi + m) > l:  # These are 0
                            continue

                        alpha_p = 0.  # Reset values
                        alpha_m = 0.
                        beta = 0.
                        if mi == 1:  # Find the alpha+ terms
                            alpha_p = -1. * (1. if li == 1 else 0.) * np.sqrt(
                                (l + m + 1.) * (l + m + 2.) / ((2. * l + 1.) * (2. * l + 3.))) \
                                      + (1. if li == -1 else 0.) * np.sqrt(
                                (l - m) * (l - m - 1.) / ((2. * l - 1.) * (2. * l + 1.)))
                        elif mi == 0:  # Find the beta term
                            beta = (1. if li == 1 else 0.) * np.sqrt(
                                (l - m + 1.) * (l + m + 1.) / ((2. * l + 1.) * (2. * l + 3.))) \
                                   + (1. if li == -1 else 0.) * np.sqrt(
                                (l - m) * (l + m) / ((2. * l - 1.) * (2. * l + 1.)))
                        else:  # Find the alpha- terms
                            alpha_m = 1. * (1. if li == 1 else 0.) * np.sqrt(
                                (l - m + 1.) * (l - m + 2.) / ((2. * l + 1.) * (2. * l + 3.))) \
                                      - (1. if li == -1 else 0.) * np.sqrt(
                                (l + m) * (l + m - 1.) / ((2. * l - 1.) * (2. * l + 1.)))

                        # Numerically calculate the radial integral
                        r_num_list = np.linspace(0, 2.996291714390912, 100)
                        r_integral = np.trapz(r_num_list ** (nu - 1.) * np.exp(-self.kappa * r_num_list) * sp.spherical_jn(l + li, p_t * r_num_list), r_num_list)

                        #r_integral = self.radial_splines[l + li](p_t)

                        # Now calculate the contribution from the specific l', m'
                        res_prime += r_integral * sph_harm(px, py, pz_t, p_t, l + li, m + mi) * (-1.j) ** (l + li) \
                                     * (mu_x / 2. * (alpha_p + alpha_m) + mu_y / (2. * 1.j) * (
                                    alpha_p - alpha_m) + mu_z * beta)

                d_res += clm * res_prime

        return d_res


    cpdef double complex d0_test(self, double p, double theta, double phi, double complex ts):
        """
        Function for test input 
        """
        return 27. * self.Ef(ts)


    cpdef double complex d_asymp_stark(self, double p, double theta, double phi, double complex ts,
                                         clm_array):
        """
        Asymptotic prefactor using Stark shift correction. This is essentialy the prefactor without expanding using 
        identity for times from the saddle-point.
        """
        cdef double sn, theta_t, px, py, pz
        cdef double complex factor1, factor2, pz_t, p_t
        cdef double complex d_res = 0.

        # Values needed
        cdef double nu = self.Z / self.kappa
        cdef int max_l = clm_array.shape[1]

        # Find the coordinates
        px = p * sin_re(theta) * cos_re(phi)
        py = p * sin_re(theta) * sin_re(phi)
        pz = p * cos_re(theta)
        #sn = 1. if self.Af(ts).imag > 0. else -1.
        #pz_t = 1j * sn * sqrt_re(2 * self.Ip + px ** 2 + py ** 2)  # pz+s.Af(ts) : tilde{pz}
        #p_t = 1j * sqrt_re(2 * self.Ip)  # sqrt(px**2+py**2+pz_2**2) : modulus{tilde{p}}
        pz_t = pz + self.Af(ts)
        p_t = (pz_t**2 + px**2 + py**2)**(1./2.)

        # Start calculating the prefactor itself! Everything not depending on l,m:
        factor1 = 2.**(-nu-1.) * self.DS(p, theta, phi, ts)**(-nu) * self.kappa**nu

        #print('factor1 : ', factor1)
        #print('DS : ', self.DS(p, theta, phi, ts)**(-nu))
        #print('kappa : ', self.kappa**nu)

        # Everything depending on l,m:
        for l in range(0, max_l):
            if np.sum(np.abs(clm_array[0,l,:])) + np.sum(np.abs(clm_array[1,l,:])) == 0:
                continue  #Skip this row if it's all 0

            factor2 = (p_t / (1.j * self.kappa))**l * 2**(-l-0.5) * gamma(l+nu+2.) / gamma(l+3./2.) \
                      * sp.hyp2f1(0.5*(l-nu), 0.5*(l-nu+1.), l+3./2., -p_t**2/self.kappa**2)

            for m in range(-l, l + 1):
                if self.OAM != 1000:  # Possibility for OAM selection - only calculate m values matching OAM
                    if m != self.OAM:
                        continue

                # Get clm
                sign = 0 if m >= 0 else 1
                clm = clm_array[sign, l, abs(m)]
                if abs(clm) == 0:  # Don't calculate if it's zero anyway...
                    continue

                d_res += factor2 * clm * self.sph_harm_OAM(px, py, pz_t, p_t, l, m, phi)

        if self.OAM == 1000:  # OAM selection is not activated
            return d_res * factor1 * 1.j
        else:  # OAM selection is activated - remember the i**OAM!
            return d_res * factor1 * 1.j * 1.j ** self.OAM


    cpdef double complex asymp_transform(self, double p, double theta, double phi, double t, clm_array):
        """
        Fourier transform of the asymptotic wave function 
        """
        cdef double px, py, pz, p_t, pz_t
        cdef double complex factor1, factor2
        cdef double nu = self.Z / self.kappa
        cdef double complex res = 0

        # Coordinates
        px = p * sin_re(theta) * cos_re(phi)
        py = p * sin_re(theta) * sin_re(phi)
        pz = p * cos_re(theta)
        pz_t = pz + self.Af(t)
        p_t = sqrt_re(px**2. + py**2. + pz_t**2.)

        # Factor without m,l
        factor1 = self.kappa ** nu * (self.kappa ** 2 + p_t ** 2) ** (-nu - 1.)

        # Factors with m,l
        for l in range(0, clm_array.shape[1]):
            if np.sum(np.abs(clm_array[0, l, :])) + np.sum(np.abs(clm_array[1, l, :])) == 0:
                continue

            factor2 = gamma(l + nu + 2.) / gamma(l + 3. / 2.) * (p_t / self.kappa) ** l \
                      * sp.hyp2f1(0.5 * (l - nu), 0.5 * (l - nu + 1.), l + 3. / 2.,
                                  -p_t ** 2 / self.kappa ** 2) * 2. ** (-l - 0.5) * (-1.j) ** l

            for m in range(-l, l + 1):
                if self.OAM != 1000:  # Possibility for OAM selection - only calculate m values matching OAM
                    if m != self.OAM:
                        continue

                sign = 0 if m >= 0 else 1
                clm = clm_array[sign, l, abs(m)]
                res += clm * self.sph_harm_OAM(px, py, pz_t, p_t, l, m, phi) * factor2

        if self.OAM == 1000:
            return res * factor1
        else:
            return res * factor1 * 1.j**self.OAM


    cpdef double complex d0(self, double p, double theta, double phi, double complex ts, double complex[:,:,:] state_array, double[:] alpha_list):
        """
        Function to select the right prefactor based on self.target 
        """
        if self.target == "GTO":
            return self.d_gto(p, theta, phi, ts, state_array)
        elif self.target == 'GTO_dress':
            return self.d_gto_dress(p, theta, phi, ts, state_array)
        elif self.target == 'asymp':
            return self.d_asymp_Er(p, theta, phi, ts, state_array)
        elif self.target == 'asymp_martiny':
            return self.d_asymp_martiny(p, theta, phi, ts, state_array)
        elif self.target == 'dipole':
            return self.d_dipole(p, theta, phi, ts, alpha_list, state_array)
        elif self.target == 'cutoff':
            return self.d_test_cutoff(p, theta, phi, ts, state_array, alpha_list)
        elif self.target == 'dipole2':
            return self.d_dip(p, theta, phi, ts, state_array)
        elif self.target == 'dress_dip':
            return self.d_gto_dress(p, theta, phi, ts, state_array[0]) + self.d_dip(p, theta, phi, ts, state_array[1])
        elif self.target == 'test':
            return self.d0_test(p, theta, phi, ts)
        elif self.target == 'single_clm':
            return self.d_asymp_lm(p, theta, phi, ts, state_array)
        else:
            return 1.


    @cython.boundscheck(False)  # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    #@functools.lru_cache(maxsize=cacheSize)
    cpdef double complex M(s, double p, double theta, double phi, state_array, alpha_list=None,
                           double tf = np.inf):  # double pz, double px, double t, int N, int eLim):
        '''
        Final transition amplitude
        Constructed as sum 
        '''
        cdef double complex MSum = 0.
        #cdef double complex d = 0.
        #cdef double x_a, y_a, z_a, alpha, front_factor
        #cdef int i, j, k
        #cdef double complex ts = 0.
        #cdef double px = p * sin_re(theta) * cos_re(phi)
        #cdef double py = p * sin_re(theta) * sin_re(phi)
        #cdef double complex pz = p * cos_re(theta)

        #if s.target == 'GTO_MO_SPA':
        #    for row in state_array:
        #        front_factor, alpha, i, j, k, x_a, y_a, z_a = row
        #        times = s.mo_times_gen(p, theta, phi, z_a)
        #        for ts in times:
        #            if real(ts) < tf:
        #                d = s.di_gto(p, theta, phi, ts, front_factor, alpha, i, j, k, x_a, y_a, z_a)
        #                det = sqrt(2. * Pi * I1 / s.DDPhi(p, theta, phi, ts, z_a))
        #                exp_phi = np.exp(I1 * s.S(p, theta, phi, ts) - I1 * (px * x_a + py * y_a + (pz + s.Af(ts)) * z_a))
        #                MSum += det * d * exp_phi
        #    return MSum

        if s.target == 'asymp_stark':
            times = s.stark_times_gen(p, theta, phi, alpha_list)

            for ts in times:
                det = sqrt(2. * Pi * I1 / s.DDS_stark(p, theta, phi, ts, alpha_list))
                expS = exp(I1 * s.S_stark(p, theta, phi, ts, alpha_list))
                d0 = s.d_asymp_stark(p, theta, phi, ts, state_array)

                #print(det)
                #print(expS)
                #print(d0)
                MSum += d0 * det * expS

            # Remember boundary terms!
            return MSum - s.asymp_transform(p, theta, phi, 0., state_array) \
                   * (np.exp(1.j * s.S_stark(p, theta, phi, 2. * s.N * np.pi / s.omega, alpha_list)) -
                      np.exp(1.j * s.S_stark(p, theta, phi, 0., alpha_list)))

        else:  # Normal stuff
            times = s.TimesGen(p, theta, phi)
            for ts in times:
                if (real(ts) < tf):
                    det = sqrt(2. * Pi * I1 / s.DDS(p, theta, phi, ts))
                    expS = exp(I1 * s.S(p, theta, phi, ts))
                    d0 = s.d0(p, theta, phi, ts, state_array, alpha_list)
                    MSum += d0 * det * expS
            return MSum


    # Transition amplitude in cartesian co-ordinates
    cpdef double complex Mxy(s, px, py, pz, state_array, alpha_list=None, tf = np.inf):
        cdef double p = sqrt_re(px * px + py * py + pz * pz)
        cdef double theta = acos_re(pz / p)
        cdef double phi = atan2_re(py, px)
        return s.M(p, theta, phi, state_array, alpha_list, tf)


    # List comprehension over cartesian transition amplitude
    @cython.boundscheck(False)  # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def Mxy_List(s, pxList, pyList, double pz, tf = np.inf, state_array=None, alpha_list=None):
        return np.array([s.Mxy(px, py, pz, state_array, alpha_list, tf) for px, py in zip(pxList, pyList)])


    def Mxz_List(s, pxList, double py, pzList, state_array, alpha_list = None, tf = np.inf):
        return np.array([s.Mxy(px, py, pz, state_array, alpha_list, tf) for px, pz in zip(pxList, pzList)])


    #### Code for exact integration of the prefactor! ####
    cpdef double complex d0_analytic_hyd1s(self, double p, double theta, double phi, double t):
        """
        Analytic transition amplitude for ground state of hydrogen, see eq. 75 + 76 in Milosevic review. 
        ONLY works for numerical calculations, as divergent in saddle points! 
        """
        return 1. / (sqrt_re(2.) * np.pi) / self.DS(p, theta, phi, t)


    cpdef double complex d0_asymp_martiny(self, double p, double theta, double phi, double t, clm_array):
        """
        Numerical version of prefactor from Martiny's thesis, found above
        """
        return self.asymp_transform(p, theta, phi, t, clm_array) \
               * self.DS(p, theta, phi, t) #* np.exp(1.j * self.S(p, theta, phi, t))


    cpdef double complex d_pol_num(self, double p, double theta, double phi, double t, clm_array, alpha_list):
        """
        Numerical calculation of the cutoff for the polarized prefactor
        """
        cdef double sn, px, py, pz, alpha_p, alpha_m, beta, pz_t, p_t,
        cdef double complex res_prime
        cdef double complex d_res = 0.

        # Values needed
        cdef double nu = self.Z / self.kappa
        cdef int max_l = clm_array.shape[1]
        cdef double E_field = self.Ef(t)
        cdef double r_integral

        # Find dipole through polarizability
        cdef double mu_x = alpha_list[0] * E_field
        cdef double mu_y = alpha_list[1] * E_field
        cdef double mu_z = alpha_list[2] * E_field

        # Find the coordinates
        px = p * sin_re(theta) * cos_re(phi)
        py = p * sin_re(theta) * sin_re(phi)
        pz = p * cos_re(theta)
        pz_t = pz + E_field
        p_t = np.sqrt(px**2 + py**2 + pz_t**2)


        # Loops for the rest (l,m and l',m')
        for l in range(0, max_l):
            for m in range(-l, l + 1):
                # Get clm
                sign = 0 if m >= 0 else 1
                clm = clm_array[sign, l, abs(m)]
                if abs(clm) == 0:  # Don't calculate if it's zero anyway...
                    continue

                # Find l' and m' terms
                res_prime = 0.
                for li in [-1, 1]:
                    if li + l < 0:  # These are 0
                        continue
                    for mi in [-1, 0, 1]:
                        if abs(mi + m) > l:  # These are 0
                            continue

                        alpha_p = 0.  # Reset values
                        alpha_m = 0.
                        beta = 0.
                        if mi == 1:  # Find the alpha+ terms
                            alpha_p = -1. * (1. if li == 1 else 0.) * np.sqrt(
                                (l + m + 1.) * (l + m + 2.) / ((2. * l + 1.) * (2. * l + 3.))) \
                                      + (1. if li == -1 else 0.) * np.sqrt(
                                (l - m) * (l - m - 1.) / ((2. * l - 1.) * (2. * l + 1.)))
                        elif mi == 0:  # Find the beta term
                            beta = (1. if li == 1 else 0.) * np.sqrt(
                                (l - m + 1.) * (l + m + 1.) / ((2. * l + 1.) * (2. * l + 3.))) \
                                   + (1. if li == -1 else 0.) * np.sqrt(
                                (l - m) * (l + m) / ((2. * l - 1.) * (2. * l + 1.)))
                        else:  # Find the alpha- terms
                            alpha_m = 1. * (1. if li == 1 else 0.) * np.sqrt(
                                (l - m + 1.) * (l - m + 2.) / ((2. * l + 1.) * (2. * l + 3.))) \
                                      - (1. if li == -1 else 0.) * np.sqrt(
                                (l + m) * (l + m - 1.) / ((2. * l - 1.) * (2. * l + 1.)))

                        # Numerically calculate the radial integral
                        #r_num_list = np.linspace(0, 2.996291714390912, 100)
                        #r_integral = np.trapz(r_num_list**(nu-1.)*np.exp(-self.kappa*r_num_list)*sp.spherical_jn(l+li, p_t * r_num_list), r_num_list)

                        r_integral = self.radial_splines[l+li](p_t)

                        # Now calculate the contribution from the specific l', m'
                        res_prime += r_integral * sph_harm(px, py, pz_t, p_t, l+li, m+mi) * (-1.j)**(l+li) \
                                     * (mu_x/2.*(alpha_p + alpha_m) + mu_y/(2.*1.j)*(alpha_p - alpha_m) + mu_z*beta)

                d_res += clm * res_prime

        return d_res


    cpdef double complex d0_test_num(self, double p, double theta, double phi, double t):
        """
        Function for test input 
        """
        return 27. * self.Ef(t)


    cpdef double complex d0_num(self, double p, double theta, double phi, double t, state_array=None, alpha_list=None):
        """
        Function to select the right prefactor based on self.target 
        """
        if self.target == "GTO":
            return self.d_gto(p, theta, phi, t, state_array)
        elif self.target == "GTO_dress":
            return self.d_gto_dress(p, theta, phi, t, state_array)
        elif self.target == 'asymp_martiny':
            return self.d0_asymp_martiny(p, theta, phi, t, state_array)
        elif self.target == 'cutoff':
            return self.d_pol_num(p, theta, phi, t, state_array, alpha_list)
        elif self.target == 'hyd1s_analytic':
            return self.d0_analytic_hyd1s(p, theta, phi, t)
        elif self.target == 'asymp_stark':
            return self.d_asymp_stark(p, theta, phi, t, state_array)
        elif self.target == 'test':
            return self.d0_test_num(p, theta, phi, t)
        else:
            return 1.


    cpdef double M_integrand(self, double t, double p, double theta, double phi, int cmplx, state_array=None, alpha_list=None):
        """
        The integrand of the time integral for numerical integration.
        Currently only the exponentiated action.
        """
        cdef double complex prefactor = self.d0_num(p, theta, phi, t, state_array, alpha_list=alpha_list)
        cdef double complex exp_S
        if self.target != 'asymp_stark':
            exp_S = np.exp(I1 * self.S(p, theta, phi, t))  # Normal exp of action
        else:
            exp_S = np.exp(I1 * self.S_stark(p, theta, phi, t, alpha_list))  # Modified stark action

        cdef double complex Mi = exp_S * prefactor

        if cmplx == 1:
            return np.imag(Mi)
        else:
            return np.real(Mi)
        #return Mi

    cpdef double complex M_num(s, double p, double theta, double phi, double tf = np.inf,
                               state_array=None, alpha_list=None):
        '''
        Final transition amplitude computed numerically 
        '''
        cdef double M_im, M_re, er1, er2
        cdef double complex bound1, bound2
        cdef double complex M

        #bound1 = 1. / (1j * (s.Ip + 0.5 * p * p))
        #bound2 = np.exp(1j * s.N * Pi * (8 * s.Ip + 4 * p * p + 3 * s.Up) / (4 * s.omega)) / (
        #            1j * (s.Ip + 0.5 * p * p))

        M_im = it.quad(s.M_integrand, 0., 2 * s.N * Pi / s.omega, args=(p, theta, phi, 1, state_array, alpha_list),
                       limit=1000, epsrel=1.e-8)[0]
        M_re = it.quad(s.M_integrand, 0., 2 * s.N * Pi / s.omega, args=(p, theta, phi, 0, state_array, alpha_list),
                       limit=1000,
                       epsrel=1.e-8)[0]

        #cdef double dt = 0.1
        #cdef double[:] t_list = np.arange(0., 2*s.N*Pi/s.omega, dt)
        #M_list = [s.M_integrand(t, p, theta, phi, 0, state_array, alpha_list, r_c) for t in t_list]
        #cdef double complex[:] M_list = np.zeros_like(t_list, dtype=complex)

        #cdef int list_len = t_list.shape[0]
        #cdef int i
        #for i in range(list_len):
        #    M_list[i] = s.M_integrand(t_list[i], p, theta, phi, 0, state_array, alpha_list, r_c)

        #M = np.trapz(M_list, t_list)

        if s.target == 'hyd1s_analytic':
            return 1.j * (M_re + 1.j * M_im) - 1 / (sqrt_re(2.) * np.pi) * s.DS(p, theta, phi, 0.) ** (-2.) \
                   * (np.exp(1.j * s.S(p, theta, phi, 2. * s.N * np.pi / s.omega)) - np.exp(
                1.j * s.S(p, theta, phi, 0.)))

        elif s.target == 'asymp_martiny':
            return 1.j * (M_re + 1.j * M_im) - s.asymp_transform(p, theta, phi, 0., state_array) \
                   * (np.exp(1.j * s.S(p, theta, phi, 2. * s.N * np.pi / s.omega)) - np.exp(
                1.j * s.S(p, theta, phi, 0.)))

        elif s.target == 'asymp_stark':
            return (M_re + 1.j * M_im) - s.asymp_transform(p, theta, phi, 0., state_array) \
                   * (np.exp(1.j * s.S_stark(p, theta, phi, 2. * s.N * np.pi / s.omega, alpha_list)) - np.exp(
                1.j * s.S_stark(p, theta, phi, 0., alpha_list)))

        else:
            return M_re + 1j * M_im  #+ bound1 - bound2
            #return M

    # Transition amplitude in cartesian co-ordinates
    cpdef double complex Mxy_num(s, px, py, pz, tf = np.inf, state_array=None, alpha_list=None):
        cdef double p = sqrt_re(px * px + py * py + pz * pz)
        cdef double theta = acos_re(pz / p)
        cdef double phi = atan2_re(py, px)
        return s.M_num(p, theta, phi, tf, state_array, alpha_list)

    def Mxz_List_num(s, pxList, double py, pzList, state_array=None, alpha_list=None, tf = np.inf):
        return np.array([s.Mxy_num(px, py, pz, tf, state_array, alpha_list) for px, pz in zip(pxList, pzList)])


    ####   ---   OAM Functions   ---   ####
    cpdef Ml(s, double p, double theta, state_array=None, int Nphi = 250):
        """
        This is the fourier series coeiffint of M to get the OAM distribusion.
        It is computed taking advantage of the FFT
        """
        phiList = np.linspace(-Pi, Pi, Nphi)
        MphiList = [s.M(p, theta, phi, state_array=state_array, alpha_list=np.array([0.])) for phi in phiList]
        return np.fft.fft(MphiList) / Nphi

    @cython.boundscheck(False)  # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def Ml_List(s, pList, theta, Nphi = 250):
        return np.array([[abs(M) ** 2 for M in s.Ml(p, theta, Nphi)] for p in pList]).T

    cpdef Mlxz(s, px, pz, state_array=None, int Nphi = 250):
        """
        convert Ml to cartesian coordinates, note px is the perpendicular coordinate st. px^2 = px^2+py^2
        """
        cdef double p = sqrt_re(px * px + pz * pz)
        cdef double theta = acos_re(pz / p)
        return s.Ml(p, theta, state_array, Nphi)

    @cython.boundscheck(False)  # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def Mlxz_List(s, pxList, pzList, state_array, Nphi = 250):
        return np.array([[abs(M) ** 2 for M in s.Mlxz(px, pz, state_array, Nphi)] for px, pz in zip(pxList, pzList)])



    #####   ---   code for spectra

    cpdef double spec_integrand(self, double phi, double theta, double p, state_array=None, err=1e-4):
        return abs(self.M(p, theta, phi, state_array=state_array))**2*sin_re(theta)

    cpdef double spec_phi(self, double theta, double p, state_array=None, err=1e-4):
        return it.quad(self.spec_integrand, 0.0, 2.0*Pi, epsabs=err, epsrel=err, limit=500,
                      args=(theta, p, state_array, err))[0]

    cpdef double hemisphere_integrated(self, double E, state_array=None, err=1e-4):
        p = np.sqrt(2.0*E)
        return it.quad(self.spec_phi, 0.0, Pi/2.0, epsabs=err, epsrel=err, limit=500,
                       args=(p, state_array, err))[0]

    cpdef double Spectra(s, double E, state_array=None, double phi = 0., double t = np.inf,
                         double err = 1.0e-4, int limit = 500):
        """Function for the spectra"""
        Norm_val, Norm_error = it.quad(s.Spec_Norm, 0, Pi, args=(phi, E, t, state_array),
                                       epsabs=err, epsrel=err, limit=limit)
        return Norm_val

    #Spectra Norm integrand
    cpdef double Spec_Norm(s, double theta, double phi, double E, double t = np.inf, state_array=None):
        """Function to compute the integrand of the theta integral for the 'spectra"""
        #phrased in spherical coordinates
        #cdef double complex M, Mg
        cdef double px, pr
        pr = sqrt_re(2 * E)
        px = pr * sin_re(theta)

        return abs(s.M(pr, theta, phi, state_array=state_array)) ** 2 * px * 2 * Pi

