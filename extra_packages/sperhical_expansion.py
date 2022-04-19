# %%
import numpy as np
import matplotlib.pyplot as plt
import pyshtools as pysh
from scipy.special import sph_harm
from scipy.special import lpmv as assoc_legendre
from . import OutputInterface
from scipy.signal import find_peaks
import scipy.special as sp
from scipy.optimize import curve_fit
# %%


def eval_GTOs(x, y, z, param_list):
    """
    Evaluate gaussian type orbitals from parameter list
    """
    res = 0
    for params in param_list:
        MO, alpha, i, j, k, x0, y0, z0 = params
        xi = x-x0
        yi = y-y0
        zi = z-z0
        res += MO * xi**i * yi**j * zi**k * np.exp(-alpha * (xi**2 + yi**2 + zi**2))
    return res


def spherical_expansion(func, N, plot_coeff=False):
    """
    Expands a given function of theta and phi in spherical harmonics (f(theta, phi)).
    Output is on the form C_(l,m) = cilm[0,l,m] and C_(l,-m) = cilm[1,l,m].
    """
    if N % 2 != 0:
        print('N should be an even number! Incrementing by one')
        N += 1

    phi_list = np.arange(0, 360, 360/N)
    theta_list = np.arange(0, 180, 180/N)
    func_grid = np.zeros((N, N), dtype=complex)

    theta_rad = np.deg2rad(theta_list)
    phi_rad = np.deg2rad(phi_list)

    for i, theta in enumerate(theta_rad):
        for j, phi in enumerate(phi_rad):
            func_grid[i, j] = func(theta, phi)
    grid = pysh.SHGrid.from_array(func_grid,  copy=True)
    sh_grid = grid.expand(normalization='ortho', csphase=-1)

    if plot_coeff:
        fig, ax = sh_grid.plot_spectrum2d()
        plt.show()

    return sh_grid.to_array()


def eval_sph_from_coeff(theta, phi, coeff_array):
    """
    Evaluate a linear combination of spherical harmonics from the array of coefficients
    at some given direction
    """
    max_l = coeff_array.shape[1]
    res = 0 + 0j
    for l in range(max_l):
        for m in range(-l, l + 1, 1):
            if m >= 0:
                sign = 0
            else:
                sign = 1
            res += coeff_array[sign, l, abs(m)] * sph_harm(m, l, phi, theta)
    return res


def get_asymp_from_sph_coeff(GTO_sph_coeff, r, Ip, Z=1):
    """
    Gets the coeffciencts for the asymptotic expansion, given the coefficients for the sperhical expansion.
    The matching is made at a given value of r.
    """
    kappa = np.sqrt(2*Ip)

    flm_list = np.zeros_like(GTO_sph_coeff)
    max_l = flm_list.shape[1]
    radial_part = np.exp(kappa*r) / r**(Z/kappa-1)

    for l in range(max_l):
        if l == 0:
            flm_list[0, 0, 0] = GTO_sph_coeff[0, 0, 0] * radial_part
            continue
        for m in range(-l, l+1, 1):
            if m >= 0:
                sign = 0
            else:
                sign = 1
            flm_list[sign, l, abs(m)] = GTO_sph_coeff[sign, l, abs(m)] * radial_part
    return flm_list


def get_as_coeffs(func, r, n_samp, Ip, Z=1, abs_thresh=1e-6, normalized=False):
    """
    Get the asymptotic coefficients for a given value of r in a.u.

    :param func: Function to expand, func(r, theta, phi).
    :param r: The r to determine the asymptotic coefficients at.
    :param n_samp: Number of points used in the spherical expansion. Determines the accuracy.
    :param Ip: Ionization potential.
    :param Z: Charge of the leftover core.
    :param abs_thresh: Threshold value of coeffs. in the expansion. All below this is set to 0.
    :param normalized: If True, the coefficients will be normalized.
    """
    flm_lst = spherical_expansion(lambda theta, phi: func(r, theta, phi), n_samp, plot_coeff=False)

    clm_lst = np.zeros_like(flm_lst, dtype=complex)
    l_max = flm_lst.shape[1]
    kappa = np.sqrt(2 * abs(Ip))
    radial = r ** (Z / kappa - 1) * np.exp(-kappa * r)
    for l in range(l_max):
        for m in range(-l, l + 1):
            sgn = 0 if m >= 0 else 1
            clm = flm_lst[sgn, l, abs(m)] / radial
            clm_lst[sgn, l, abs(m)] = clm if abs(clm) > abs_thresh else 0

    if normalized:
        return clm_lst / np.sum(np.abs(clm_lst)**2)
    else:
        return clm_lst


def get_asymptotic_coeffs(func, n_r, n_samp, Ip, Z=1, interval=None, plot=False, normalized=False):
    """
    Gets the coefficients for the asymptotic wave function from the function
    """
    if interval is None:
        interval = [2, 17.5]

    # First find flms as a function of r
    r_lst = np.linspace(interval[0], interval[-1], n_r)
    flm_lst = []
    for i, r in enumerate(r_lst):
        print(f'Evaluating at r={r:.4f} \t Nr. {i + 1}/{n_r}')
        flm_lst.append(spherical_expansion(lambda theta, phi: func(r, theta, phi), n_samp, plot_coeff=False))
    flm_lst = np.array(flm_lst)

    # Loop through all of the coeffiicients and find the constant clm.
    # If it is lower than a threshold, it is set to zero
    ABS_THRESH = 1e-2
    clm_lst = np.zeros_like(flm_lst[0])
    l_max = flm_lst.shape[2]
    kappa = np.sqrt(2*np.abs(Ip))
    radial = lambda r, k: r**(Z/k - 1) * np.exp(-k*r)

    for l in range(l_max):
        for m in range(-l, l + 1):
            sgn = 0 if m >= 0 else 1
            clm = flm_lst[:, sgn, l, abs(m)] / radial(r_lst, kappa)
            idx = find_peaks(np.abs(clm))[0]
            if idx.size > 0:
                print(f'l={l} \t m={m} \t {r_lst[idx]} \t {clm[idx]}')
                val = clm[idx[-1]]
                clm_lst[sgn, l, abs(m)] = (val if np.abs(val) > ABS_THRESH else 0)
                if plot:
                    FUSK_FACTOR = 70
                    plt.figure(facecolor='white')
                    plt.plot(r_lst, np.abs(clm), label=r'$c_{\ell m}$')
                    plt.plot(r_lst, FUSK_FACTOR * np.abs(flm_lst[:, sgn, l, abs(m)]), label=r'$f_{\ell m}$')
                    plt.xlabel(r'Radius $r$')
                    plt.ylabel(r'Absolute amplitude')
                    plt.plot(r_lst[idx], np.abs(clm[idx]), 'o')
                    plt.legend(frameon=False)
                    plt.minorticks_on()
                    plt.show()

    if normalized:
        return clm_lst / np.sum(np.abs(clm_lst)**2)
    else:
        return clm_lst


def convert_list_to_clm_array(coeff_list, fill_to=None):
    """
    Converts a flat array of Clm coeffs into the form of the Clm array used in this script: array[sign, l, abs(m)],
    with sign = 1 means -m and sign = 1 is +m. Note that the list must have the correct number of entries to match a
    l value. The fill_to parameter enables filling up to a given l with dummy variable -1 when past given values.
    """
    if fill_to is None:  # Here we match the array size to the given coeff_list
        max_l = np.sqrt(len(coeff_list)) - 1
    else:  # Here we calculate for a fixed size array
        max_l = int(fill_to) - 1

    if max_l % 1 != 0:
        print('Invalid length of list!')
        return None

    if type(coeff_list[0]) == np.complex128 or type(coeff_list[0]) == complex:
        clm_array = np.zeros((2, int(max_l) + 1, int(max_l) + 1), dtype=complex)
    else:
        clm_array = np.zeros((2, int(max_l) + 1, int(max_l) + 1))
    counter = 0

    for l in range(int(max_l + 1)):
        for m in range(-l, l+1):
            sign = 0 if m >= 0 else 1

            if counter <= len(coeff_list)-1:  # Use the value given
                clm_array[sign, l, abs(m)] = coeff_list[counter]
            else:  # Place dummy variable -1
                clm_array[sign, l, abs(m)] = -1
            counter += 1

    return clm_array


def get_as_from_r_array(func, r_array, n_samp, Ip, Z=1, normalized=False):
    """
    Calculates the asymptotic coeffs of func(r, theta, phi) in the given r values.
    Input is a array of r values in same form as the Clm arrays used in this script.
    Use convert_list_to_clm_array() first, if you have a flat list of values.
    """
    max_l = r_array.shape[1]
    clm_array = np.zeros_like(r_array, dtype=complex)
    kappa = np.sqrt(2 * abs(Ip))

    old_r = -1  # To make sure we don't keep calculating the Laplace expansion if r is the same
    for l in range(max_l):
        print(f'Calculating for l = {l} and m =', end='')
        for m in range(-l, l+1):
            print(f' {m},', end='')
            sign = 0 if m >= 0 else 1
            r = r_array[sign, l, abs(m)]

            if r != old_r:  # r is different - find new Laplace expansion!
                flm_list = spherical_expansion(lambda theta, phi: func(r, theta, phi), n_samp, plot_coeff=False)
                radial = r ** (Z / kappa - 1) * np.exp(-kappa * r)
                old_r = r

            clm_array[sign, l, abs(m)] = flm_list[sign, l, abs(m)] / radial  # Get the asymptotic coefficient

        print('\n')
    if normalized:
        return clm_array / np.sum(np.abs(clm_array) ** 2)
    else:
        return clm_array


def convert_clm_array_to_list(clm_array):
    """
    Converts an array of the clm_array shape to a list (sorted by l and then m)
    """
    max_l = clm_array.shape[1]
    clm_list = []

    for l in range(max_l):
        for m in range(-l, l+1):
            sign = 0 if m >= 0 else 1
            clm_list.append(clm_array[sign, l, abs(m)])

    return clm_list


def replace_dummy_variables(func, clm_array, r_value, Ip, Z=1):
    """
    Replace the dummy values ((-1) from convert_list_to_clm_array) in the clm_array with values matched in r
    """
    max_l = clm_array.shape[1]
    res_array = np.zeros_like(clm_array, dtype=complex)
    flm_lst = spherical_expansion(lambda theta, phi: func(r_value, theta, phi), max_l*2, plot_coeff=False)
    kappa = np.sqrt(2 * abs(Ip))

    radial = r_value ** (Z / kappa - 1) * np.exp(-kappa * r_value)
    for l in range(max_l):
        for m in range(-l, l + 1):
            sign = 0 if m >= 0 else 1
            clm = clm_array[sign, l, abs(m)]

            if clm != -1:
                res_array[sign, l, abs(m)] = clm
            else:
                res_array[sign, l, abs(m)] = flm_lst[sign, l, abs(m)] / radial

    return res_array


def extend_to_higher_l(func, clm_array, new_l, r_value, Ip, Z=1):
    """
    Extends the clm_array by matching the new coefficients at the given r value
    """
    if new_l <= clm_array.shape[1] - 1:
        print('The new l value is less or equal to the original! Returning original array...')
        return clm_array

    clm_list = convert_clm_array_to_list(clm_array)
    array_extended = convert_list_to_clm_array(clm_list, fill_to=new_l)
    clm_array_extended = replace_dummy_variables(func, array_extended, r_value, Ip, Z)
    return clm_array_extended


def get_asymp_fit(func, r_list, n_samp, Ip, orbital_nr=None, Z=1, return_flm=False, threshold=None):
    """
    Gets the asymptotic coefficients from a function func(r, theta, phi), matched to the asymptotic WF on the interval
    r_list.
    """
    kappa = np.sqrt(2*Ip)

    # First get the flm's
    f_lms = []
    for i, r in enumerate(r_list):
        print(f'\rEvaluating at r={r:.4f} \t Nr. {i + 1}/{len(r_list)}', end='')
        if orbital_nr is None:
            f_lms.append(spherical_expansion(lambda theta, phi: func(r, theta, phi), n_samp))
        else:
            f_lms.append(spherical_expansion(lambda theta, phi: func(r, theta, phi, orbital_nr), n_samp))
    f_lms = np.array(f_lms)
    print('')  # Just to break the line

    # Then do the fitting to find the asymptotic coefficients
    def asymp_ting(r, clm_real, clm_imag):
        # This is just a function to split up the real and imaginary parts...
        r_cal = r[:len(r)//2]
        res_real = clm_real * r_cal**(Z/kappa - 1) * np.exp(-kappa*r_cal)
        res_imag = clm_imag * r_cal**(Z/kappa - 1) * np.exp(-kappa*r_cal)
        return np.hstack([res_real, res_imag])

    clm_array = np.zeros_like(f_lms[0], dtype=complex)

    print('Now fitting!')
    for l in range(clm_array.shape[1]):
        for m in range(-l, l+1):
            sign = 0 if m >= 0 else 1

            flm_real = np.real(f_lms[:, sign, l, abs(m)])
            flm_imag = np.imag(f_lms[:, sign, l, abs(m)])
            popt, _ = curve_fit(asymp_ting, np.hstack([r_list, r_list]), np.hstack([flm_real, flm_imag]), p0=[1, 1])

            if threshold is None:
                clm_array[sign, l, abs(m)] = popt[0] + 1j*popt[1]
            else:
                if abs(popt[0] + 1j*popt[1]) > threshold:
                    clm_array[sign, l, abs(m)] = popt[0] + 1j*popt[1]
                else:
                    clm_array[sign, l, abs(m)] = 0.

    print('Done!')
    if return_flm:
        return clm_array, f_lms
    else:
        return clm_array


def laplace_several_r(func, r_list, n_samp, orbital_nr=None):
    """Gets the Laplace expansion coeffs. for several r values"""
    f_lms = []
    for i, r in enumerate(r_list):
        print(f'\rEvaluating at r={r:.4f} \t Nr. {i + 1}/{len(r_list)}', end='')
        if orbital_nr is None:
            f_lms.append(spherical_expansion(lambda theta, phi: func(r, theta, phi), n_samp))
        else:
            f_lms.append(spherical_expansion(lambda theta, phi: func(r, theta, phi, orbital_nr), n_samp))
    return np.array(f_lms)


def eval_asymptotic_cart(x, y, z, coeffs, Ip, Z=1):
    """
    Evaluates the asymptotic wave function in cartesian coordinates
    """
    l_max = coeffs.shape[1]
    kappa = np.sqrt(2 * abs(Ip))
    eta = 2 * Z / kappa + 5
    radial_norm = np.sqrt((2 * kappa) ** eta / sp.gamma(eta))

    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(np.sqrt(x**2 + y**2), z)
    phi = np.arctan2(y, x)

    if type(x) is np.ndarray:
        radial_part, angular_sum = np.zeros_like(x, dtype=complex), np.zeros_like(x, dtype=complex)
    else:
        radial_part, angular_sum = 0, 0

    radial_part += radial_norm * r**(Z / kappa - 1) * np.exp(-kappa * r)
    for l in range(l_max):
        for m in range(-l, l + 1):
            sgn = 1 if m >= 0 else 0  # This should be the other way around?
            angular_sum += coeffs[sgn, l, m]*sp.sph_harm(m, l, phi, theta)

    return radial_part * angular_sum


def eval_asymptotic(r, theta, phi, coeff_array, Ip, Z=1):
    """
    Evaluates the asymptotic wave function in spherical coordinates from the c_lm array of coefficients
    """
    max_l = coeff_array.shape[1]
    res = 0 + 0j
    kappa = np.sqrt(2*Ip)
    radial_part = np.exp(-kappa*r) * r**(Z/kappa-1)
    for l in range(max_l):
        if l == 0:
            res += coeff_array[0, 0, 0] * sph_harm(0, 0, phi, theta) * radial_part
            continue
        for m in range(-l, l + 1, 1):
            if m >= 0:
                sign = 0
            else:
                sign = 1
            res += coeff_array[sign, l, abs(m)] * sph_harm(m, l, phi, theta) * radial_part
    return res


def cylindrical_from_spherical(r_par, r_perp, sph_coeffs):
    """
    Calculates the Fourier coefficients of a function given as a spherical expansion.
    Note that if sph_coeffs depend on r, this is valid for this value of r only (r = sqrt(r_par**2 + r_perp**2))!
    """
    r = np.sqrt(r_perp**2 + r_par**2)
    cos_theta = r_par / r
    theta = np.arccos(cos_theta)
    max_m = sph_coeffs.shape[1]-1

    fm_list = []
    for m in range(-max_m, max_m+1):
        sign = 0 if m >= 0 else 1
        fm = 0
        for l, flm in enumerate(sph_coeffs[sign, :, abs(m)]):
            if abs(m) > l or flm == 0 or flm == 0j:
                continue
            #N_lm = np.sqrt((2*l+1)/(4*np.pi) * np.math.factorial(l-m)/np.math.factorial(l+1))  # Should have these precalculated?
            #fm += flm * N_lm * assoc_legendre(m, l, cos_theta)
            fm += flm * sph_harm(m, l, 0, theta)  # This might be better?
        fm_list.append(fm)
    return fm_list


def eval_cylindrical(phi, coeff_list):
    """
    Evaluates a function given as a Fourier series with coefficients in coeff_list
    """
    res = 0
    max_m = int((len(coeff_list) - 1)/2)
    for m, coeff in zip(range(-max_m, max_m+1), coeff_list):
        res += np.exp(1j * m * phi) * coeff
    return res

