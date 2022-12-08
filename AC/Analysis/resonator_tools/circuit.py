import warnings
import logging
import numpy as np
import time
import scipy.optimize as spopt
from scipy.constants import hbar
from scipy.interpolate import splrep, splev
from scipy.ndimage.filters import gaussian_filter1d

# from Geresdi_lab_code.AC.Analysis import resonator_tools

from Geresdi_lab_code.AC.Analysis.resonator_tools.utilities import plotting, save_load, Watt2dBm, dBm2Watt
from Geresdi_lab_code.AC.Analysis.resonator_tools.circlefit import circlefit
from Geresdi_lab_code.AC.Analysis.resonator_tools.calibration import calibration
import importlib

##
## z_data_raw denotes the raw data
## z_data denotes the normalized data
##		  

class reflection_port_original(circlefit, save_load, plotting, calibration):
    '''
    normal direct port probed in reflection
    '''

    def __init__(self, f_data=None, z_data_raw=None):
        self.porttype = 'direct'
        self.fitresults = {}
        self.z_data = None
        if f_data is not None:
            self.f_data = np.array(f_data)
        else:
            self.f_data = None
        if z_data_raw is not None:
            self.z_data_raw = np.array(z_data_raw)
        else:
            self.z_data = None
        self.phasefitsmooth = 3


    def _S11(self, f, fr, k_c, k_i):
        '''
        use either frequency or angular frequency units
        for all quantities
        k_l=k_c+k_i: total (loaded) coupling rate
        k_c: coupling rate
        k_i: internal loss rate
        '''
        return ((k_c - k_i) + 2j * (f - fr)) / ((k_c + k_i) - 2j * (f - fr))


    def get_delay(self, f_data, z_data, delay=None, ignoreslope=True, guess=True):
        '''
        ignoreslope option not used here
        retrieves the cable delay assuming the ideal resonance has a circular shape
        modifies the cable delay until the shape Im(S21) vs Re(S21) is circular
        see "do_calibration"
        '''
        maxval = np.max(np.absolute(z_data))
        z_data = z_data / maxval
        A1, A2, A3, A4, fr, Ql = self._fit_skewed_lorentzian(f_data, z_data)
        if self.df_error / fr > 0.0001 or self.dQl_error / Ql > 0.1:
            # print("WARNING: Calibration using Lorentz fit failed, trying phase fit...")
            A1 = np.mean(np.absolute(z_data))
            A2 = 0.
            A3 = 0.
            A4 = 0.
            # fr = np.mean(f_data)
            f = splrep(f_data, np.unwrap(np.angle(z_data)), k=5, s=self.phasefitsmooth)
            fr = f_data[np.argmax(np.absolute(splev(f_data, f, der=1)))]
            Ql = 1e4
        if ignoreslope == True:
            A2 = 0.
        else:
            A2 = 0.
            print(
                "WARNING: The ignoreslope option is ignored! Corrections to the baseline should be done manually prior to fitting.")
            print("see also: resonator_tools.calibration.fit_baseline_amp() etc. for help on fitting the baseline.")
            print("There is also an example ipython notebook for using this function.")
            print(
                "However, make sure to understand the impact of the baseline (parasitic coupled resonances etc.) on your system.")
        # z_data = (np.absolute(z_data)-A2*(f_data-fr)) * np.exp(np.angle(z_data)*1j)  #usually not necessary
        if delay is None:
            if guess == True:
                delay = self._guess_delay(f_data, z_data)
            else:
                delay = 0.
            delay = self._fit_delay(f_data, z_data, delay, maxiter=300)
        params = [A1, A2, A3, A4, fr, Ql]
        return delay, params


    def do_calibration(self, f_data, z_data, ignoreslope=True, guessdelay=True, fixed_delay=None):
        '''
        calculating parameters for normalization
        '''
        # if highQi:
        #     delay, params = self.manual_fit_phase(f_data, z_data, delay=fixed_delay, ignoreslope=ignoreslope, guesses=None)
        #else:
        delay, params = self.get_delay(f_data, z_data, ignoreslope=ignoreslope, guess=guessdelay, delay=fixed_delay)
        z_data = (z_data - params[1] * (f_data - params[4])) * np.exp(2. * 1j * np.pi * delay * f_data)
        xc, yc, r0 = self._fit_circle(z_data)
        zc = np.complex(xc, yc)

        fitparams = self._phase_fit(f_data, self._center(z_data, zc), 0., np.absolute(params[5]), params[4])
        theta, Ql, fr = fitparams
        beta = self._periodic_boundary(theta + np.pi, np.pi)  ###
        offrespoint = np.complex((xc + r0 * np.cos(beta)), (yc + r0 * np.sin(beta)))

        alpha = self._periodic_boundary(np.angle(offrespoint) + np.pi, np.pi)
        # a = np.absolute(offrespoint)
        # alpha = np.angle(zc)
        a = r0 + np.absolute(zc)
        return delay, a, alpha, fr, Ql, params[1], params[4]


    def do_normalization(self, f_data, z_data, delay, amp_norm, alpha, A2, frcal):
        '''
        transforming resonator into canonical position
        '''
        return (z_data - A2 * (f_data - frcal)) / amp_norm * np.exp(1j * (-alpha + 2. * np.pi * delay * f_data))


    def circlefit(self, f_data, z_data, fr=None, Ql=None, refine_results=False, calc_errors=True):
        '''
        S11 version of the circlefit
        '''

        if fr is None: fr = f_data[np.argmin(np.absolute(z_data))]
        if Ql is None: Ql = 1e6
        xc, yc, r0 = self._fit_circle(z_data, refine_results=refine_results)
        phi0 = -np.arcsin(yc / r0)
        theta0 = self._periodic_boundary(phi0 + np.pi, np.pi)
        z_data_corr = self._center(z_data, np.complex(xc, yc))
        theta0, Ql, fr = self._phase_fit(f_data, z_data_corr, theta0, Ql, fr)
        # print("Ql from phasefit is: " + str(Ql))
        Qi = Ql / (1. - r0)
        Qc = 1. / (1. / Ql - 1. / Qi)

        results = {"Qi": Qi, "Qc": Qc, "Ql": Ql, "fr": fr, "theta0": theta0}

        # calculation of the error
        p = [fr, Qc, Ql]
        # chi_square, errors = rt.get_errors(rt.residuals_notch_ideal,f_data,z_data,p)
        if calc_errors == True:
            chi_square, cov = self._get_cov_fast_directrefl(f_data, z_data, p)
            # chi_square, cov = rt.get_cov(rt.residuals_notch_ideal,f_data,z_data,p)

            if cov is not None:
                errors = np.sqrt(np.diagonal(cov))
                fr_err, Qc_err, Ql_err = errors
                # calc Qi with error prop (sum the squares of the variances and covariaces)
                dQl = 1. / ((1. / Ql - 1. / Qc) ** 2 * Ql ** 2)
                dQc = - 1. / ((1. / Ql - 1. / Qc) ** 2 * Qc ** 2)
                Qi_err = np.sqrt(
                    (dQl ** 2 * cov[2][2]) + (dQc ** 2 * cov[1][1]) + (2 * dQl * dQc * cov[2][1]))  # with correlations
                errors = {"Ql_err": Ql_err, "Qc_err": Qc_err, "fr_err": fr_err, "chi_square": chi_square, "Qi_err": Qi_err}
                results.update(errors)
            else:
                print("WARNING: Error calculation failed!")
        else:
            # just calc chisquared:
            fun2 = lambda x: self._residuals_notch_ideal(x, f_data, z_data) ** 2
            chi_square = 1. / float(len(f_data) - len(p)) * (fun2(p)).sum()
            errors = {"chi_square": chi_square}
            results.update(errors)

        return results


    def autofit(self, electric_delay=None, fcrop=None):
        '''
        automatic calibration and fitting
        electric_delay: set the electric delay manually
        fcrop = (f1,f2) : crop the frequency range used for fitting

        '''
        if fcrop is None:
            self._fid = np.ones(self.f_data.size, dtype=bool)
        else:
            f1, f2 = fcrop
            self._fid = np.logical_and(self.f_data >= f1, self.f_data <= f2)

        delay, amp_norm, alpha, fr, Ql, A2, frcal = \
            self.do_calibration(self.f_data[self._fid],
                                self.z_data_raw[self._fid],
                                ignoreslope=True,
                                guessdelay=False,
                                fixed_delay=electric_delay,
                                )
        self.z_data = self.do_normalization(self.f_data, self.z_data_raw, delay, amp_norm, alpha, A2, frcal)
        self.fitresults = self.circlefit(self.f_data[self._fid], self.z_data[self._fid], fr, Ql, refine_results=False,
                                         calc_errors=True)
        self.fitresults['delay'] = delay
        self.fitresults['a'] = amp_norm
        self.fitresults['alpha'] = alpha
        self.z_data_sim = A2 * (self.f_data - frcal) + self._S11_directrefl(self.f_data, fr=self.fitresults["fr"],
                                                                            Ql=self.fitresults["Ql"],
                                                                            Qc=self.fitresults["Qc"], a=amp_norm,
                                                                            alpha=alpha, delay=delay)
        self.z_data_sim_norm = self._S11_directrefl(self.f_data, fr=self.fitresults["fr"], Ql=self.fitresults["Ql"],
                                                    Qc=self.fitresults["Qc"], a=1., alpha=0., delay=0.)
        self._delay = delay



    def _S11_directrefl(self, f, fr=10e9, Ql=900, Qc=1000., a=1., alpha=0., delay=.0):
        '''
        full model for notch type resonances
        '''
        return a * np.exp(np.complex(0, alpha)) * np.exp(-2j * np.pi * f * delay) * (
                    2. * Ql / Qc - 1. + 2j * Ql * (fr - f) / fr) / (1. - 2j * Ql * (fr - f) / fr)


    def get_single_photon_limit(self, unit='dBm'):
        '''
        returns the amout of power in units of W necessary
        to maintain one photon on average in the cavity
        unit can be 'dbm' or 'watt'
        '''
        if self.fitresults != {}:
            fr = self.fitresults['fr']
            k_c = 2 * np.pi * fr / self.fitresults['Qc']
            k_i = 2 * np.pi * fr / self.fitresults['Qi']
            if unit == 'dBm':
                return Watt2dBm(1. / (4. * k_c / (2. * np.pi * hbar * fr * (k_c + k_i) ** 2)))
            elif unit == 'watt':
                return 1. / (4. * k_c / (2. * np.pi * hbar * fr * (k_c + k_i) ** 2))

        else:
            warnings.warn('Please perform the fit first', UserWarning)
            return None


    def get_photons_in_resonator(self, power, unit='dBm'):
        '''
        returns the average number of photons
        for a given power (defaul unit is 'dbm')
        unit can be 'dBm' or 'watt'
        '''
        if self.fitresults != {}:
            if unit == 'dBm':
                power = dBm2Watt(power)
            fr = self.fitresults['fr']
            k_c = 2 * np.pi * fr / self.fitresults['Qc']
            k_i = 2 * np.pi * fr / self.fitresults['Qi']
            return 4. * k_c / (2. * np.pi * hbar * fr * (k_c + k_i) ** 2) * power
        else:
            warnings.warn('Please perform the fit first', UserWarning)
            return None

    def _phase_dist(self, angle):
        """
        Maps angle [-2pi, +2pi] to phase distance on circle [0, pi]
        """
        return np.pi - np.abs(np.pi - np.abs(angle))

    @classmethod
    def phase_centered(cls, f, fr, Ql, theta, delay=0.):
        """
        Yields the phase response of a strongly overcoupled (Qi >> Qc) resonator
        in reflection which corresponds to a circle centered around the origin.
        Additionally, a linear background slope is accounted for if needed.

        inputs:
        - fr: Resonance frequency
        - Ql: Loaded quality factor (and since Qi >> Qc also Ql = Qc)
        - theta: Offset phase
        - delay (opt.): Time delay between output and input signal leading to
                        linearly frequency dependent phase shift
        """
        return theta - 2 * np.pi * delay * (f - fr) + 2. * np.arctan(2. * Ql * (1. - f / fr))
    ###############################
    def manual_fit_phase(self, f_data, z_data, delay=None, ignoreslope=True, guesses=None):
        """
        Fits the phase response of a strongly overcoupled (Qi >> Qc) resonator
        in reflection which corresponds to a circle centered around the origin
        (cf‌. phase_centered()).

        inputs:
        - z_data: Scattering data of which the phase should be fit. Data must be
                  distributed around origin ("circle-like").
        - guesses (opt.): If not given, initial guesses for the fit parameters
                          will be determined. If given, should contain useful
                          guesses for fit parameters as a tuple (fr, Ql, delay)

        outputs:
        - fr: Resonance frequency
        - Ql: Loaded quality factor
        - theta: Offset phase
        - delay: Time delay between output and input signal leading to linearly
                 frequency dependent phase shift
        """
        phase = np.unwrap(np.angle(z_data))

        # For centered circle roll-off should be close to 2pi. If not warn user.
        if np.max(phase) - np.min(phase) <= 0.8 * 2 * np.pi:
            logging.warning(
                "Data does not cover a full circle (only {:.1f}".format(
                    np.max(phase) - np.min(phase)
                )
                + " rad). Increase the frequency span around the resonance?"
            )
            roll_off = np.max(phase) - np.min(phase)
        else:
            roll_off = 2 * np.pi
        #roll_off = np.max(phase) - np.min(phase)

        # Set useful starting parameters
        if guesses is None:
            # Use maximum of derivative of phase as guess for fr
            # phase_smooth = splrep(self.f_data, phase, k=5, s=100)
            # phase_derivative = splev(self.f_data, phase_smooth, der=1)
            phase_smooth = gaussian_filter1d(phase, 30)
            phase_derivative = np.gradient(phase_smooth)
            fr_guess = self.f_data[np.argmax(np.abs(phase_derivative))]
            Ql_guess = 2 * fr_guess / (self.f_data[-1] - self.f_data[0])
            # Estimate delay from background slope of phase (substract roll-off)
            slope = phase[-1] - phase[0] + roll_off
            delay_guess = -slope / (2 * np.pi * (self.f_data[-1] - self.f_data[0]))
        else:
            fr_guess, Ql_guess, delay_guess = guesses
        # This one seems stable and we do not need a manual guess for it
        theta_guess = 0.5 * (np.mean(phase[:5]) + np.mean(phase[-5:]))

        # Fit model with less parameters first to improve stability of fit
        def residuals_Ql(params):
            Ql, = params
            return residuals_full((fr_guess, Ql, theta_guess, delay_guess))

        def residuals_fr_theta(params):
            fr, theta = params
            return residuals_full((fr, Ql_guess, theta, delay_guess))

        # def residuals_Ql_delay(params):
        # Ql, delay = params
        # return residuals_full((fr_guess, Ql, theta_guess, delay))
        def residuals_delay(params):
            delay, = params
            return residuals_full((fr_guess, Ql_guess, theta_guess, delay))

        def residuals_fr_Ql(params):
            fr, Ql = params
            return residuals_full((fr, Ql, theta_guess, delay_guess))

        # def residuals_fr(params):
        # fr, = params
        # return residuals_full((fr, Ql_guess, theta_guess, delay_guess))
        def residuals_full(params):
            return self._phase_dist(
                phase - self.phase_centered(self.f_data, *params)
            )

        p_final = spopt.leastsq(residuals_Ql, [Ql_guess])
        Ql_guess, = p_final[0]
        p_final = spopt.leastsq(residuals_fr_theta, [fr_guess, theta_guess])
        fr_guess, theta_guess = p_final[0]
        p_final = spopt.leastsq(residuals_delay, [delay_guess])
        delay_guess, = p_final[0]
        p_final = spopt.leastsq(residuals_fr_Ql, [fr_guess, Ql_guess])
        fr_guess, Ql_guess = p_final[0]
        # p_final = spopt.leastsq(residuals_fr, [fr_guess])
        # fr_guess, = p_final[0]
        # p_final = spopt.leastsq(residuals_Ql, [Ql_guess])
        # Ql_guess, = p_final[0]
        p_final = spopt.leastsq(residuals_full, [
            fr_guess, Ql_guess, theta_guess, delay_guess
        ])
        A1, A2, A3, A4, fr, Ql = self._fit_skewed_lorentzian(f_data, z_data)
        if self.df_error / fr > 0.0001 or self.dQl_error / Ql > 0.1:
            print("WARNING: Calibration using Lorentz fit failed, trying phase fit...")
            A1 = np.mean(np.absolute(z_data))
            A2 = 0.
            A3 = 0.
            A4 = 0.
            # fr = np.mean(f_data)
            f = splrep(f_data, np.unwrap(np.angle(z_data)), k=5, s=self.phasefitsmooth)
            fr = f_data[np.argmax(np.absolute(splev(f_data, f, der=1)))]
            Ql = 1e4
        if ignoreslope == True:
            A2 = 0.
        else:
            A2 = 0.
            print(
                "WARNING: The ignoreslope option is ignored! Corrections to the baseline should be done manually prior to fitting.")
            print("see also: resonator_tools.calibration.fit_baseline_amp() etc. for help on fitting the baseline.")
            print("There is also an example ipython notebook for using this function.")
            print(
                "However, make sure to understand the impact of the baseline (parasitic coupled resonances etc.) on your system.")
            # z_data = (np.absolute(z_data)-A2*(f_data-fr)) * np.exp(np.angle(z_data)*1j)  #usually not necessary
        params = [A1, A2, A3, A4, fr, Ql]

        return p_final[0][3], params

class reflection_port_phi(circlefit, save_load, plotting, calibration):
    '''
    normal direct port probed in reflection
    '''

    def __init__(self, f_data=None, z_data_raw=None):
        self.porttype = 'direct'
        self.fitresults = {}
        self.z_data = None
        if f_data is not None:
            self.f_data = np.array(f_data)
        else:
            self.f_data = None
        if z_data_raw is not None:
            self.z_data_raw = np.array(z_data_raw)
        else:
            self.z_data = None
        self.phasefitsmooth = 3

    def _S11(self, f, fr, k_c, k_i):
        '''
        use either frequency or angular frequency units
        for all quantities
        k_l=k_c+k_i: total (loaded) coupling rate
        k_c: coupling rate
        k_i: internal loss rate
        '''
        return ((k_c - k_i) + 2j * (f - fr)) / ((k_c + k_i) - 2j * (f - fr))

    def get_delay(self, f_data, z_data, delay=None, ignoreslope=True, guess=True):
        '''
        ignoreslope option not used here
        retrieves the cable delay assuming the ideal resonance has a circular shape
        modifies the cable delay until the shape Im(S21) vs Re(S21) is circular
        see "do_calibration"
        '''
        maxval = np.max(np.absolute(z_data))
        z_data = z_data / maxval
        A1, A2, A3, A4, fr, Ql = self._fit_skewed_lorentzian(f_data, z_data)

        # print(A1, A2, A3, A4, fr, Ql)
        if self.df_error / fr > 0.0001 or self.dQl_error / Ql > 0.1:
            print("WARNING: Calibration using Lorentz fit failed, trying phase fit...")
            A1 = np.mean(np.absolute(z_data))
            A2 = 0.
            A3 = 0.
            A4 = 0.
            # fr = np.mean(f_data)
            f = splrep(f_data, np.unwrap(np.angle(z_data)), k=5, s=self.phasefitsmooth)
            fr = f_data[np.argmax(np.absolute(splev(f_data, f, der=1)))]
            Ql = 1e4
        if ignoreslope == True:
            A2 = 0.
        else:
            A2 = 0.
            print(
                "WARNING: The ignoreslope option is ignored! Corrections to the baseline should be done manually prior to fitting.")
            print("see also: resonator_tools.calibration.fit_baseline_amp() etc. for help on fitting the baseline.")
            print("There is also an example ipython notebook for using this function.")
            print(
                "However, make sure to understand the impact of the baseline (parasitic coupled resonances etc.) on your system.")
            # z_data = (np.absolute(z_data)-A2*(f_data-fr)) * np.exp(np.angle(z_data)*1j)  #usually not necessary
        if delay is None:
            if guess == True:
                delay = self._guess_delay(f_data, z_data)
            else:
                delay = 0.
            delay = self._fit_delay(f_data, z_data, delay, maxiter=800)
        params = [A1, A2, A3, A4, fr, Ql]
        return delay, params

    def do_calibration(self, f_data, z_data,
                       fixed_delay=None,
                       ignoreslope=True,
                       guesses=None, guessdelay = True,
                       ):
        '''
        calculating parameters for normalization
        '''

        delay, params = self.get_delay(f_data, z_data, ignoreslope=ignoreslope, guess=guessdelay, delay=fixed_delay)

        z_data = (z_data - params[1] * (f_data - params[4])) * np.exp(2. * 1j * np.pi * delay * f_data)
        xc, yc, r0 = self._fit_circle(z_data)
        zc = np.complex(xc, yc)

        fitparams = self._phase_fit(f_data, self._center(z_data, zc), 0., np.absolute(params[5]), params[4])
        theta, Ql, fr = fitparams
        beta = self._periodic_boundary(theta + np.pi, np.pi)  ###
        #print(f'beta2 = {beta}, xc = {xc}, yc = {yc}')
        offrespoint = np.complex((xc + r0 * np.cos(beta)), (yc + r0 * np.sin(beta)))
        alpha = self._periodic_boundary(np.angle(offrespoint) + np.pi, np.pi)
        #print(f"offresspoint2 = {offrespoint}")
        # print(f"alpha = {alpha}")
        a = np.absolute(offrespoint)
        #print(a)
        # alpha = np.angle(zc)
        #a = r0 + np.absolute(zc)
        #a = np.absolute(zc)
        # print(f"do calbiration: r0 = {self.r0}, zc = {zc}")
        return delay, a, alpha, fr, Ql, params[1], params[4]

    def do_alpha_calibration(self, f_data, z_data,
                       fixed_delay=None,
                       ignoreslope=True,
                       guesses=None, guessdelay = True,
                       ):
        '''
        calculating parameters for normalization
        '''
        #if highQi:
        #    delay, params = self.manual_fit_phase(f_data, z_data, delay=fixed_delay, ignoreslope=ignoreslope,
        #                                          guesses=None)
        #else:
        delay, params = self.get_delay(f_data, z_data, ignoreslope=ignoreslope, guess=guessdelay, delay=fixed_delay)
        #print(f'delay = {delay}')
        z_data = (z_data - params[1] * (f_data - params[4])) * np.exp(2. * 1j * np.pi * delay * f_data)
        xc, yc, r0 = self._fit_circle(z_data)
        zc = np.complex(xc, yc)

        fitparams = self._phase_fit(f_data, self._center(z_data, zc), 0., np.absolute(params[5]), params[4])
        theta = fitparams[0]
        beta = self._periodic_boundary(theta + np.pi, np.pi)  ###
        offrespoint = np.complex((xc + r0 * np.cos(beta)), (yc + r0 * np.sin(beta)))
        #print(offrespoint)
        alpha = self._periodic_boundary(np.angle(offrespoint) + np.pi, np.pi)
        return delay, alpha

    def do_normalization(self, f_data, z_data, delay, amp_norm, alpha, A2, frcal):
        '''
        transforming resonator into canonical position
        '''
        # print(-A2*(f_data-frcal))
        return (z_data - A2 * (f_data - frcal)) / amp_norm * np.exp(1j * (-alpha + 2. * np.pi * delay * f_data))

    def circlefit(self, f_data, z_data, fr=None, Ql=None,
                  refine_results=False,
                  calc_errors=True,
                  fit_type = 'DCM'):
        '''
        S11 version of the circlefit
        '''

        if fr is None: fr = f_data[np.argmin(np.absolute(z_data))]
        if Ql is None: Ql = 1e4

        xc, yc, r0 = self._fit_circle(z_data, refine_results=refine_results)
        #phi0 = -np.arcsin(yc / r0) Original
        # phi0 = np.arcsin(yc / r0) Correct (?) adjustment
        phi0 = np.arctan(yc / (xc + 1))
        # print(f"yc = {yc}, xc = {xc},  r0 = {r0}")
        if xc < -1:
            self.phi = phi0 + np.pi
        else:
            self.phi = phi0
        #theta0 = self._periodic_boundary(-phi0 + np.pi, np.pi)
        theta0 = self._periodic_boundary(-phi0 + np.pi, np.pi)
        z_data_corr = self._center(z_data, np.complex(xc, yc))
        theta0, Ql, fr = self._phase_fit(f_data, z_data_corr, theta0, Ql, fr)
        # print("Ql from phasefit is: " + str(Ql))
        if fit_type == 'DCM':
            #print("Doing DCM fit")
            #Qi = Ql / (1. - r0)
            #Qc = 1. / (1. / Ql - 1. / Qi)
            Qc = Ql / (r0)
            Qc_complex = Qc * np.exp(np.complex(0, phi0))
            #print(1/Qc)
            #print(1/Qc_complex)
            Qi = 1. / (1. / Ql - 1. / Qc)
            #Qi = 1. / (1. / Ql - np.real(1/ Qc_complex))

        else:
            #print("Doing phi CM fit")
            Qi = Ql / (1. - r0)
            Qc = 1. / (1. / Ql - 1. / Qi)
            #Qc = Ql / (r0)
            #Qi = 1. / (1. / Ql - 1/ Qc)

        results = {"Qi": Qi, "Qc": Qc, "Ql": Ql, "fr": fr, "theta0": theta0}

        # calculation of the error
        p = [fr, Qc, Ql]
        # chi_square, errors = rt.get_errors(rt.residuals_notch_ideal,f_data,z_data,p)
        if calc_errors == True:
            chi_square, cov = self._get_cov_fast_directrefl(f_data, z_data, p)
            # chi_square, cov = rt.get_cov(rt.residuals_notch_ideal,f_data,z_data,p)

            if cov is not None:
                errors = np.sqrt(np.diagonal(cov))
                fr_err, Qc_err, Ql_err = errors
                # calc Qi with error prop (sum the squares of the variances and covariaces)
                dQl = 1. / ((1. / Ql - 1. / Qc) ** 2 * Ql ** 2)
                dQc = - 1. / ((1. / Ql - 1. / Qc) ** 2 * Qc ** 2)
                Qi_err = np.sqrt(
                    (dQl ** 2 * cov[2][2]) + (dQc ** 2 * cov[1][1]) + (2 * dQl * dQc * cov[2][1]))  # with correlations
                errors = {"Ql_err": Ql_err, "Qc_err": Qc_err, "fr_err": fr_err, "chi_square": chi_square,
                          "Qi_err": Qi_err}
                results.update(errors)
            else:
                print("WARNING: Error calculation failed!")
        else:
            # just calc chisquared:
            fun2 = lambda x: self._residuals_notch_ideal(x, f_data, z_data) ** 2
            chi_square = 1. / float(len(f_data) - len(p)) * (fun2(p)).sum()
            errors = {"chi_square": chi_square}
            results.update(errors)

        return results

    def autofit(self, electric_delay=None, fcrop=None, manual_calibrate = False,  fit_type = 'DCM'):
        '''
        automatic calibration and fitting
        Manual calibrate works best with a lot of delay AND when Qi >> Qc. Then, set electric delay DOES NOT DO ANYTHING.
        fit_type has to get some extra love still
        electric_delay: set the electric delay manually if man_calib = False
        fcrop = (f1,f2) : crop the frequency range used for fitting
        '''

        if fcrop is None:
            self._fid = np.ones(self.f_data.size, dtype=bool)
        else:
            f1, f2 = fcrop
            self._fid = np.logical_and(self.f_data >= f1, self.f_data <= f2)

        if manual_calibrate:
            delay, amp_norm, alpha, fr, Ql, A2, frcal = \
                self._manual_calibrate(self.f_data[self._fid], self.z_data_raw[self._fid], ignoreslope=True,
                                       guessdelay=False, fixed_delay=electric_delay)

            #print(amp_norm, alpha, delay, fr, Ql, A2,)
#            delay, alpha = self.do_alpha_calibration(self.f_data[self._fid],
#                                    self.z_data_raw[self._fid],
#                                    ignoreslope=True,
#                                    guessdelay=False,
#                                    )
        else:
            delay, amp_norm, alpha, fr, Ql, A2, frcal = \
                self.do_calibration(self.f_data[self._fid],
                                    self.z_data_raw[self._fid],
                                    ignoreslope=True,
                                    guessdelay=False,
                                    fixed_delay=electric_delay,
                                    )

        #print( amp_norm, alpha, delay, fr, Ql, A2 )

        self.z_data = self.do_normalization(self.f_data, self.z_data_raw, delay, amp_norm, alpha, A2, frcal)
        self.fitresults = self.circlefit(self.f_data[self._fid], self.z_data[self._fid], fr, Ql, refine_results=True,
                                         calc_errors=True, fit_type = fit_type)
        self.fitresults['delay'] = delay
        self.fitresults['a'] = amp_norm
        self.fitresults['alpha'] = alpha
        self.fitresults['phi'] = self.phi
        self.z_data_sim = A2 * (self.f_data - frcal) + self._S11_directrefl(self.f_data, fr=self.fitresults["fr"],
                                                                            Ql=self.fitresults["Ql"],
                                                                            Qc=self.fitresults["Qc"], a=-amp_norm,
                                                                            alpha=alpha, delay=delay)
        self.z_data_sim_norm = self._S11_directrefl(self.f_data, fr=self.fitresults["fr"], Ql=self.fitresults["Ql"],
                                                    Qc=self.fitresults["Qc"], a=1., alpha=0., delay=0.)
        self._delay = delay

    def _S11_directrefl(self, f, fr=10e9, Ql=900, Qc=1000., a=1., alpha=0., delay=.0):
        '''
        full model for notch type resonances
        '''
        phi = self.phi
        #phi = np.pi-1
        complexQc = Qc * np.exp(-1j * phi)
        return a * np.exp(1j * (alpha - 2 * np.pi * f * delay)) * (
                1. - 2. * Ql / (complexQc * 1 * (1. + 2j * Ql * (f / fr - 1.)))
        )
        # return a*np.exp(np.complex(0,alpha))*np.exp(-2j*np.pi*f*delay) * ( 2.*Ql/complexQc - 1. + 2j*Ql*(fr-f)/fr ) / ( 1. - 2j*Ql*(fr-f)/fr )

    def get_single_photon_limit(self, unit='dBm'):
        '''
        returns the amout of power in units of W necessary
        to maintain one photon on average in the cavity
        unit can be 'dbm' or 'watt'
        '''
        if self.fitresults != {}:
            fr = self.fitresults['fr']
            k_c = 2 * np.pi * fr / self.fitresults['Qc']
            k_i = 2 * np.pi * fr / self.fitresults['Qi']
            if unit == 'dBm':
                return Watt2dBm(1. / (4. * k_c / (2. * np.pi * hbar * fr * (k_c + k_i) ** 2)))
            elif unit == 'watt':
                return 1. / (4. * k_c / (2. * np.pi * hbar * fr * (k_c + k_i) ** 2))

        else:
            warnings.warn('Please perform the fit first', UserWarning)
            return None

    def get_photons_in_resonator(self, power, unit='dBm'):
        '''
        returns the average number of photons
        for a given power (defaul unit is 'dbm')
        unit can be 'dBm' or 'watt'
        '''
        if self.fitresults != {}:
            if unit == 'dBm':
                power = dBm2Watt(power)
            fr = self.fitresults['fr']
            k_c = 2 * np.pi * fr / self.fitresults['Qc']
            k_i = 2 * np.pi * fr / self.fitresults['Qi']
            return 4. * k_c / (2. * np.pi * hbar * fr * (k_c + k_i) ** 2) * power
        else:
            warnings.warn('Please perform the fit first', UserWarning)
            return None

    ##############################
    # Added some stuff of myself #
    ##############################

    def coupling_C(self, Z_r, Z_0=50):
        Q_e = self.fitresults.get('Qc')
        Q_inter = self.fitresults.get('Qi')
        w_r = 2 * np.pi * self.fitresults.get('fr')

        C_couple = np.sqrt(np.pi / (Z_0 * Z_r * Q_e)) * 1 / (w_r)
        # Critical coupling is when Q_i = Q_c
        C_crit = np.sqrt(np.pi / (Z_0 * Z_r * Q_inter)) * 1 / (w_r)
        self.fitresults.update({
            "C_crit": C_crit,
            "C_couple": C_couple
        })

    def manual_fit_phase(self, f_data, z_data, delay=None, ignoreslope=True, guesses=None):
        """
        Fits the phase response of a strongly overcoupled (Qi >> Qc) resonator
        in reflection which corresponds to a circle centered around the origin
        (cf‌. phase_centered()).

        inputs:
        - z_data: Scattering data of which the phase should be fit. Data must be
                  distributed around origin ("circle-like").
        - guesses (opt.): If not given, initial guesses for the fit parameters
                          will be determined. If given, should contain useful
                          guesses for fit parameters as a tuple (fr, Ql, delay)

        outputs:
        - fr: Resonance frequency
        - Ql: Loaded quality factor
        - theta: Offset phase
        - delay: Time delay between output and input signal leading to linearly
                 frequency dependent phase shift
        """
        phase = np.unwrap(np.angle(z_data))

        # For centered circle roll-off should be close to 2pi. If not warn user.
        if np.max(phase) - np.min(phase) <= 0.8 * 2 * np.pi:
            logging.warning(
                "Data does not cover a full circle (only {:.1f}".format(
                    np.max(phase) - np.min(phase)
                )
                + " rad). Increase the frequency span around the resonance?"
            )
            roll_off = np.max(phase) - np.min(phase)
        else:
            roll_off = 2 * np.pi

        # Set useful starting parameters
        if guesses is None:
            # Use maximum of derivative of phase as guess for fr
            # phase_smooth = splrep(self.f_data, phase, k=5, s=100)
            # phase_derivative = splev(self.f_data, phase_smooth, der=1)
            phase_smooth = gaussian_filter1d(phase, 30)
            phase_derivative = np.gradient(phase_smooth)
            fr_guess = self.f_data[np.argmax(np.abs(phase_derivative))]
            Ql_guess = 2 * fr_guess / (self.f_data[-1] - self.f_data[0])
            # Estimate delay from background slope of phase (substract roll-off)
            slope = phase[-1] - phase[0] + roll_off
            delay_guess = -slope / (2 * np.pi * (self.f_data[-1] - self.f_data[0]))
        else:
            fr_guess, Ql_guess, delay_guess = guesses
        # This one seems stable and we do not need a manual guess for it
        theta_guess = 0.5 * (np.mean(phase[:5]) + np.mean(phase[-5:]))

        # Fit model with less parameters first to improve stability of fit

        def residuals_Ql(params):
            Ql, = params
            return residuals_full((fr_guess, Ql, theta_guess, delay_guess))

        def residuals_fr_theta(params):
            fr, theta = params
            return residuals_full((fr, Ql_guess, theta, delay_guess))

        # def residuals_Ql_delay(params):
        # Ql, delay = params
        # return residuals_full((fr_guess, Ql, theta_guess, delay))
        def residuals_delay(params):
            delay, = params
            return residuals_full((fr_guess, Ql_guess, theta_guess, delay))

        def residuals_fr_Ql(params):
            fr, Ql = params
            return residuals_full((fr, Ql, theta_guess, delay_guess))

        # def residuals_fr(params):
        # fr, = params
        # return residuals_full((fr, Ql_guess, theta_guess, delay_guess))
        def residuals_full(params):
            return self._phase_dist(
                phase - self.phase_centered(self.f_data, *params)
            )

        p_final = spopt.leastsq(residuals_Ql, [Ql_guess])
        Ql_guess, = p_final[0]
        p_final = spopt.leastsq(residuals_fr_theta, [fr_guess, theta_guess])
        fr_guess, theta_guess = p_final[0]
        p_final = spopt.leastsq(residuals_delay, [delay_guess])
        delay_guess, = p_final[0]
        p_final = spopt.leastsq(residuals_fr_Ql, [fr_guess, Ql_guess])
        fr_guess, Ql_guess = p_final[0]
        # p_final = spopt.leastsq(residuals_fr, [fr_guess])
        # fr_guess, = p_final[0]
        # p_final = spopt.leastsq(residuals_Ql, [Ql_guess])
        # Ql_guess, = p_final[0]
        p_final = spopt.leastsq(residuals_full, [
            fr_guess, Ql_guess, theta_guess, delay_guess
        ])

        A1, A2, A3, A4, fr, Ql = self._fit_skewed_lorentzian(f_data, z_data)
        if self.df_error / fr > 0.0001 or self.dQl_error / Ql > 0.1:
            print("WARNING: Calibration using Lorentz fit failed, trying phase fit...")
            A1 = np.mean(np.absolute(z_data))
            A2 = 0.
            A3 = 0.
            A4 = 0.
            # fr = np.mean(f_data)
            f = splrep(f_data, np.unwrap(np.angle(z_data)), k=5, s=self.phasefitsmooth)
            fr = f_data[np.argmax(np.absolute(splev(f_data, f, der=1)))]
            Ql = 1e4
        if ignoreslope == True:
            A2 = 0.
        else:
            A2 = 0.
            print(
                "WARNING: The ignoreslope option is ignored! Corrections to the baseline should be done manually prior to fitting.")
            print("see also: resonator_tools.calibration.fit_baseline_amp() etc. for help on fitting the baseline.")
            print("There is also an example ipython notebook for using this function.")
            print(
                "However, make sure to understand the impact of the baseline (parasitic coupled resonances etc.) on your system.")
            # z_data = (np.absolute(z_data)-A2*(f_data-fr)) * np.exp(np.angle(z_data)*1j)  #usually not necessary
        params = [A1, A2, A3, A4, fr, Ql]

        return p_final[0][3], params

    def _manual_fit_delay(self):
        """
        Finds the cable delay by repeatedly centering the "circle" and fitting
        the slope of the phase response.
        """

        # Translate data to origin
        xc, yc, r0 = self._fit_circle(self.z_data_raw)
        z_data = self.z_data_raw - complex(xc, yc)

        # Find first estimate of parameters
        fr, Ql, theta, self.delay = self._fit_phase(f_data, z_data)

        # Do not overreact (see end of for loop)
        self.delay *= 0.05

        # Iterate to improve result for delay
        for i in range(self.fit_delay_max_iterations):
            # Translate new best fit data to origin
            z_data = self.z_data_raw * np.exp(2j * np.pi * self.delay * self.f_data)
            xc, yc, r0 = self._fit_circle(z_data)
            z_data -= complex(xc, yc)

            # Find correction to current delay
            guesses = (fr, Ql, 5e-11)
            fr, Ql, theta, delay_corr = self._fit_phase(z_data, guesses)

            # Stop if correction would be smaller than "measurable"
            phase_fit = self.phase_centered(self.f_data, fr, Ql, theta, delay_corr)
            residuals = np.unwrap(np.angle(z_data)) - phase_fit
            if 2 * np.pi * (self.f_data[-1] - self.f_data[0]) * delay_corr <= np.std(residuals):
                break

            # Avoid overcorrection that makes procedure switch between positive
            # and negative delays
            if delay_corr * self.delay < 0:  # different sign -> be careful
                if abs(delay_corr) > abs(self.delay):
                    self.delay *= 0.5
                else:
                    # delay += 0.1*delay_corr
                    self.delay += 0.1 * np.sign(delay_corr) * 5e-11
            else:  # same direction -> can converge faster
                if abs(delay_corr) >= 1e-8:
                    self.delay += min(delay_corr, self.delay)
                elif abs(delay_corr) >= 1e-9:
                    self.delay *= 1.1
                else:
                    self.delay += delay_corr

        if 2 * np.pi * (self.f_data[-1] - self.f_data[0]) * delay_corr > np.std(residuals):
            logging.warning(
                "Delay could not be fit properly!"
            )

        # Store result in dictionary (also for backwards-compatibility)
        # self.fitresults["delay"] = self.delay

    def _manual_calibrate(self, f_data, z_data, ignoreslope=True, guessdelay=True, fixed_delay=None):
        """
        Finds the parameters for normalization of the scattering data. See
        Sij of port classes for explanation of parameters.
        """

        # Correct for delay and translate circle to origin
        self.delay, params = self.manual_fit_phase(f_data, z_data)
        #print(self.delay)
        z_data = self.z_data_raw * np.exp(2j * np.pi * self.delay * self.f_data)
        xc, yc, self.r0 = self._fit_circle(z_data)
        zc = complex(xc, yc)
        #print(zc)
        z_data -= zc

        #print(xc, yc)
        # Find off-resonant point by fitting offset phase
        # (centered circle corresponds to lossless resonator in reflection)

        self.fr, self.Ql, theta, self.delay_remaining = self._fit_phase_x(z_data)
        self.theta = self._periodic_boundary(theta, np.pi)

        beta = self._periodic_boundary(theta + np.pi, np.pi)
         #print(f'beta = {beta}, xc = {xc}, yc = {yc}')
        offrespoint = zc + self.r0 * np.cos(beta) + 1j * self.r0 * np.sin(beta)
        self.offrespoint = offrespoint

        #print(f"offresspoint = {offrespoint}")
        self.a = np.absolute(offrespoint)
        #print(self.a)
        #self.alpha = np.angle(offrespoint) + np.pi
        self.alpha = self._periodic_boundary(np.angle(offrespoint) + np.pi, np.pi)
        #print(f"alpha = {np.angle(offrespoint)}?")
        self.phi = self._periodic_boundary(beta - self.alpha, np.pi)
        #print(f"phi = {self.phi}")
        # Store radius for later calculation
        self.r0 /= self.a
        # print(f"man calib: r0 = {self.r0}, zc = {zc}")
        return self.delay, self.a, self.alpha, self.fr, self.Ql, params[1], params[4]

    def _fit_phase_x(self, z_data, guesses=None):
        """
        Fits the phase response of a strongly overcoupled (Qi >> Qc) resonator
        in reflection which corresponds to a circle centered around the origin
        (cf‌. phase_centered()).

        inputs:
        - z_data: Scattering data of which the phase should be fit. Data must be
                  distributed around origin ("circle-like").
        - guesses (opt.): If not given, initial guesses for the fit parameters
                          will be determined. If given, should contain useful
                          guesses for fit parameters as a tuple (fr, Ql, delay)

        outputs:
        - fr: Resonance frequency
        - Ql: Loaded quality factor
        - theta: Offset phase
        - delay: Time delay between output and input signal leading to linearly
                 frequency dependent phase shift
        """
        phase = np.unwrap(np.angle(z_data))

        # For centered circle roll-off should be close to 2pi. If not warn user.
        if np.max(phase) - np.min(phase) <= 0.8 * 2 * np.pi:
            logging.warning(
                "Data does not cover a full circle (only {:.1f}".format(
                    np.max(phase) - np.min(phase)
                )
                + " rad). Increase the frequency span around the resonance?"
            )
            roll_off = np.max(phase) - np.min(phase)
        else:
            roll_off = 2 * np.pi

        # Set useful starting parameters
        if guesses is None:
            # Use maximum of derivative of phase as guess for fr
            # phase_smooth = splrep(self.f_data, phase, k=5, s=100)
            # phase_derivative = splev(self.f_data, phase_smooth, der=1)
            phase_smooth = gaussian_filter1d(phase, 30)
            phase_derivative = np.gradient(phase_smooth)
            fr_guess = self.f_data[np.argmax(np.abs(phase_derivative))]
            Ql_guess = 2 * fr_guess / (self.f_data[-1] - self.f_data[0])
            # Estimate delay from background slope of phase (substract roll-off)
            slope = phase[-1] - phase[0] + roll_off
            delay_guess = -slope / (2 * np.pi * (self.f_data[-1] - self.f_data[0]))
        else:
            fr_guess, Ql_guess, delay_guess = guesses
        # This one seems stable and we do not need ual guess for it
        theta_guess = 0.5 * (np.mean(phase[:5]) + np.mean(phase[-5:]))

        # Fit model with less parameters first to improve stability of fit

        def residuals_Ql(params):
            Ql, = params
            return residuals_full((fr_guess, Ql, theta_guess, delay_guess))

        def residuals_fr_theta(params):
            fr, theta = params
            return residuals_full((fr, Ql_guess, theta, delay_guess))

        # def residuals_Ql_delay(params):
        # Ql, delay = params
        # return residuals_full((fr_guess, Ql, theta_guess, delay))
        def residuals_delay(params):
            delay, = params
            return residuals_full((fr_guess, Ql_guess, theta_guess, delay))

        def residuals_fr_Ql(params):
            fr, Ql = params
            return residuals_full((fr, Ql, theta_guess, delay_guess))

        # def residuals_fr(params):
        # fr, = params
        # return residuals_full((fr, Ql_guess, theta_guess, delay_guess))
        def residuals_full(params):
            return self._phase_dist(
                phase - self.phase_centered(self.f_data, *params)
            )

        p_final = spopt.leastsq(residuals_Ql, [Ql_guess])
        Ql_guess, = p_final[0]
        p_final = spopt.leastsq(residuals_fr_theta, [fr_guess, theta_guess])
        fr_guess, theta_guess = p_final[0]
        p_final = spopt.leastsq(residuals_delay, [delay_guess])
        delay_guess, = p_final[0]
        p_final = spopt.leastsq(residuals_fr_Ql, [fr_guess, Ql_guess])
        fr_guess, Ql_guess = p_final[0]
        # p_final = spopt.leastsq(residuals_fr, [fr_guess])
        # fr_guess, = p_final[0]
        # p_final = spopt.leastsq(residuals_Ql, [Ql_guess])
        # Ql_guess, = p_final[0]
        p_final = spopt.leastsq(residuals_full, [
            fr_guess, Ql_guess, theta_guess, delay_guess
        ])

        return p_final[0]

    def _phase_dist(self, angle):
        """
        Maps angle [-2pi, +2pi] to phase distance on circle [0, pi]
        """
        return np.pi - np.abs(np.pi - np.abs(angle))

    @classmethod
    def phase_centered(cls, f, fr, Ql, theta, delay=0.):
        """
        Yields the phase response of a strongly overcoupled (Qi >> Qc) resonator
        in reflection which corresponds to a circle centered around the origin.
        Additionally, a linear background slope is accounted for if needed.
        
        inputs:
        - fr: Resonance frequency
        - Ql: Loaded quality factor (and since Qi >> Qc also Ql = Qc)
        - theta: Offset phase
        - delay (opt.): Time delay between output and input signal leading to
                        linearly frequency dependent phase shift
        """
        return theta - 2 * np.pi * delay * (f - fr) + 2. * np.arctan(2. * Ql * (1. - f / fr))


########################################################################################################################
########################################################################################################################


class resonator(object):
    '''
    Universal resonator analysis class
    It can handle different kinds of ports and assymetric resonators.
    '''

    def __init__(self, ports={}, comment=None):
        '''
        initializes the resonator class object
        ports (dictionary {key:value}): specify the name and properties of the coupling ports
            e.g. ports = {'1':'direct', '2':'notch'}
        comment: add a comment
        '''
        self.comment = comment
        self.port = {}
        self.transm = {}
        if len(ports) > 0:
            for key, pname in iter(ports.items()):
                if pname == 'direct':
                    self.port.update({key: reflection_port()})
                elif pname == 'notch':
                    self.port.update({key: notch_port()})
                else:
                    warnings.warn("Undefined input type! Use 'direct' or 'notch'.", SyntaxWarning)
        if len(self.port) == 0: warnings.warn("Resonator has no coupling ports!", UserWarning)

    def add_port(self, key, pname):
        if pname == 'direct':
            self.port.update({key: reflection_port()})
        elif pname == 'notch':
            self.port.update({key: notch_port()})
        else:
            warnings.warn("Undefined input type! Use 'direct' or 'notch'.", SyntaxWarning)
        if len(self.port) == 0: warnings.warn("Resonator has no coupling ports!", UserWarning)

    def delete_port(self, key):
        del self.port[key]
        if len(self.port) == 0: warnings.warn("Resonator has no coupling ports!", UserWarning)

    def get_Qi(self):
        '''
        based on the number of ports and the corresponding measurements
        it calculates the internal losses
        '''
        pass

    def get_single_photon_limit(self, port):
        '''
        returns the amout of power necessary to maintain one photon 
        on average in the cavity
        '''
        pass

    def get_photons_in_resonator(self, power, port):
        '''
        returns the average number of photons
        for a given power
        '''
        pass

    def add_transm_meas(self, port1, port2):
        '''
        input: port1
        output: port2
        adds a transmission measurement 
        connecting two direct ports S21
        '''
        key = port1 + " -> " + port2
        self.port.update({key: transm()})
        pass
