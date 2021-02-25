#!/bin/python3

"""
Created on Wed Jan 08 20:00:39 2014

@author: nbaiboun
"""

import scipy.signal as ss
import scipy.optimize as so
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial as facto
import sys

class LinSys(object):

    def __init__(self, *args, **kwargs):
        self.linsys = ss.lti(*args)
        if "dt" in kwargs:
            self._dtime = kwargs.get("dt")
        else:
            self._dtime = 0

    #def __repr__(self):

    @property
    def dtime(self):
        return self._dtime

    @dtime.setter
    def dtime(self, value):
        self._dtime = value

    def __add__(self, system):
        return parallel(self, system)

    def __mul__(self, system):
        return series(self, system)

def series(h1, h2):
    """
    This function computes and returns the equivalent
    transfer function of 2 transfer functions in series.
    Inputs:
        h1: first transfer function
        h2: second transfer function
        h1, h2 can be both a LinSys instance or a numeric
    Output:
        hs: LinSys instance representing the equivalent transfer function
    """

    if isinstance(h1, (int, float)):
        h1 = LinSys(h1, [1])
    if isinstance(h2, (int, float)):
        h2 = LinSys(h2, [1])

    if isinstance(h1, LinSys) and isinstance(h2, LinSys):
        num = np.polymul(h1.linsys.num, h2.linsys.num)
        den = np.polymul(h1.linsys.den, h2.linsys.den)
        dtime = h1.dtime + h2.dtime
        hs = LinSys(num, den, dt=dtime)
    else:
        print("The entries have to be LinSys instances!")
        hs = None

    return hs

def parallel(h1, h2, pade_order=2):
    """
    This function computes and returns the equivalent
    transfer function of 2 transfer functions in parallel.
    Inputs:
        h1: first transfer function
        h2: second transfer function
        pade_order: pade_order approximation of time delays
        h1, h2 can be both a LinSys instance or a numeric
    Output:
        hp: LinSys instance representing the equivalent transfer function
    """

    if isinstance(h1, (int, float)):
        h1 = LinSys(h1, [1])
    if isinstance(h2, (int, float)):
        h2 = LinSys(h2, [1])

    if isinstance(h1, LinSys) and isinstance(h2, LinSys):
        if np.isclose(h1.dtime, h2.dtime):
            num1 = np.polymul(h1.linsys.num, h2.linsys.den)
            num2 = np.polymul(h1.linsys.den, h2.linsys.num)
            num = np.polyadd(num1, num2)
            den = np.polymul(h1.linsys.den, h2.linsys.den)
            hp = LinSys(num, den, dt=h1.dtime)
        else:
            print("The result isn't a linear system!")

            h1_lin = pade(h1, N=pade_order)
            h2_lin = pade(h2, N=pade_order)

            num1 = np.polymul(h1_lin.linsys.num, h2_lin.linsys.den)
            num2 = np.polymul(h1_lin.linsys.den, h2_lin.linsys.num)
            num = np.polyadd(num1, num2)
            den = np.polymul(h1_lin.linsys.den, h2_lin.linsys.den)
            hp = LinSys(num, den)
    else:
        print('Dead times have been approximated using second order pade approximation!')
        hp = None

    return hp

def feedback(chain_d, chain_f, sign=-1, pade_order=2):
    """
    This function computes and returns the equivalent
    transfer function of 2 transfer functions in feedback system.
    Inputs:
        chain_d: transfer function of the direct chain
        chain_f: transfer function of the feedback chain
        sign: sign of the feedback chain
        pade_order: pade_order approximation of time delays
        chain_d, chain_f can be both a LinSys instance or a numeric
    Output:
        hf: LinSys instance representing the equivalent transfer function
    """

    if isinstance(chain_d, (int, float)):
        chain_d = LinSys(chain_d, [1])
    if isinstance(chain_f, (int, float)):
        chain_f = LinSys(chain_f, [1])

    if isinstance(chain_d, LinSys) and isinstance(chain_f, LinSys):
        if chain_d.dtime != 0 or chain_f.dtime != 0:
            print('Dead times have been approximated using second order pade approximation!')

        chain_d_lin = pade(chain_d, N=pade_order)
        chain_f_lin = pade(chain_f, N=pade_order)

        num = np.polymul(chain_d_lin.linsys.num, chain_f_lin.linsys.den)
        den1 = sign*np.polymul(chain_d_lin.linsys.num, chain_f_lin.linsys.num)
        den2 = np.polymul(chain_d_lin.linsys.den, chain_f_lin.linsys.den)
        den = np.polysub(den2, den1)

        hf = LinSys(num, den)
    else:
        print("The entries have to be LinSys instances!")
        hf = None

    return hf

def impulse(system, X0=None, T=None, N=None, **kwargs):
    """
    This function computes the impulse response of a linear system
    Inputs:
        system: LinSys instance
        X0: vector containing initial states of the system
        T: vector containing time values for the output to be computed
        N: number of points to be computed
    Outputs:
        t: time vector
        yout: output vector
    """

    if isinstance(system, LinSys):
        t, y = ss.impulse2(system.linsys, X0, T, N, **kwargs)
        if system.dtime != 0:
            i = np.nonzero(t > system.dtime)[0][0]
            yout = np.concatenate((np.zeros(i), y[:-i]))
        else:
            yout = y
    else:
        print("System must be a LinSys instance!")
        t, yout = None, None

    return t, yout

def step(system, X0=None, T=None, N=None, **kwargs):
    if isinstance(system, LinSys):
        t, y = ss.step2(system.linsys, X0, T, N, **kwargs)
        if system.dtime != 0:
            i = np.nonzero(t > system.dtime)[0][0]
            yout = np.concatenate((np.zeros(i), y[:-i]))
        else:
            yout = y
    else:
        print("System must be a LinSys instance!")
        t, yout = None, None

    return t, yout

def getZPK(system):
    """
    """

    if isinstance(system, LinSys):
        zpk = system.linsys.to_zpk()
        z = zpk.zeros
        p = zpk.poles
        k = zpk.gain
    else:
        print("System must be a LinSys instance!")
        z, p, k = None, None, None

    return z, p, k

def freqrange(system):
    """Computes the frequency range of the bode diagrams"""

    if isinstance(system, LinSys):
        z,p,k = getZPK(system)
        zp = np.concatenate((z,p))

        #Calcul des bornes inférieures et supérieures du spectre de fréquence utile
        f = zp[np.nonzero(zp)]
        if len(f) == 0:
            Wmin = abs(k)*0.01
            Wmax = abs(k)*10**(len(zp)+1)
        else:
            Wmin = min(abs(f))*0.01
            Wmax = max(abs(f))*100

        if system.dtime != 0:
            Wmax = max([Wmax, 18.0/system.dtime])

        w = np.logspace(np.floor(np.log10(Wmin)), np.ceil(np.log10(Wmax)), 10000)
    else:
        print("System must be a LinSys instance!")
        w = None

    return w

def bode(system, w=None, n=100, plot=True):

    if isinstance(system, LinSys):
        if w is None:
            w = freqrange(system)
        w, mag, phase = ss.bode(system.linsys, w=w, n=n)
        phase -= w * system.dtime * 180.0 / np.pi

        #Affichage
        if plot:
            fig, axs = plt.subplots(2, 1, constrained_layout=True)
            axs[0].semilogx(w, mag)
            axs[1].semilogx(w, phase)
            axs[0].set_title('Diagramme des gains')
            axs[1].set_title('Diagramme des phases')
            axs[0].grid(True, which='both')
            axs[1].grid(True, which='both')
            axs[0].set_xlim(w[0], w[-1])
            axs[0].set_ylim(min(mag)-1, max(mag)+1)
            axs[1].set_xlim(w[0], w[-1])
            axs[1].set_ylim(min(phase)-1, max(phase)+1)
            axs[0].set_ylabel('Gain [dB]')
            axs[1].set_ylabel('Phase [degrés]')
            axs[1].set_xlabel('Pulsation [rad/s]')
    else:
        print("System must be a LinSys instance!")
        w, mag, phase = None, None, None

    return w, mag, phase

def bodeasympt(system, approx=0, plot=True):

    if isinstance(system, LinSys):
        z, p, k = getZPK(system)
        zp = np.concatenate((z,p))
        zpc = np.concatenate(((np.ones(len(z))),(-1*np.ones(len(p)))))

        w = freqrange(system)
        mag = np.zeros(len(w))
        phase = np.zeros(len(w))

        #Calcul de l'effet du gain sur le diagramme de Bode
        K = abs(k*np.prod(z[np.nonzero(z)])/np.prod(p[np.nonzero(p)]))
        mag += 20*np.log10(abs(K))
        phase += np.arccos(np.sign(K))*180.0/np.pi

# Méthode de superposition appliquée sur le diagramme de Bode
        if len(zp) > 0:
            for i in range(len(zp)):
                j = np.nonzero(w > abs(zp[i]))
                if (zp[i] == 0):
                    mag[j] += zpc[i]*20*np.log10(w[j])
                else:
                    mag[j] += zpc[i]*20*np.log10(w[j]/abs(zp[i]))

                if approx == 0:
                    if zp[i] > 0:
                        phase = phase + zpc[i]*180
                        phase[j] = phase[j] - zpc[i]*90
                    else:
                        phase[j] = phase[j] + zpc[i]*90
                elif approx == 1:
                    l = np.nonzero(w >= abs(zp[i]*0.1))[0]
                    m = np.nonzero(w >= abs(zp[i]*10))[0]
                    if zp[i] > 0:
                        phase += zpc[i]*180
                        phase[l[0]:m[0]] -= zpc[i]*45*np.log10(w[l[0]:m[0]]/abs(zp[i]*0.1))
                        phase[m] -= zpc[i]*90
                    else:
                        if zp[i] == 0:
                            phase[l[0]:m[0]] += zpc[i]*45*np.log10(w[l[0]:m[0]])
                        else:
                            phase[l[0]:m[0]] += zpc[i]*45*np.log10(w[l[0]:m[0]]/abs(zp[i]*0.1))

                        phase[m] += zpc[i]*90

        #Affichage
        if plot:
            fig, axs = plt.subplots(2, 1, constrained_layout=True)
            axs[0].semilogx(w, mag)
            axs[1].semilogx(w, phase)
            axs[0].set_title('Diagramme des gains')
            axs[1].set_title('Diagramme des phases')
            axs[0].grid(True, which='both')
            axs[1].grid(True, which='both')
            axs[0].set_xlim(w[0], w[-1])
            axs[0].set_ylim(min(mag)-1, max(mag)+1)
            axs[1].set_xlim(w[0], w[-1])
            axs[1].set_ylim(min(phase)-1, max(phase)+1)
            axs[0].set_ylabel('Gain [dB]')
            axs[1].set_ylabel('Phase [degrés]')
            axs[1].set_xlabel('Pulsation [rad/s]')
    else:
        w, mag, phase = None, None, None

    return w, mag, phase

def margin(system):
    if isinstance(system, LinSys):
        w = freqrange(system)
        w0 = 10**((np.log10(w[-1]) + np.log10(w[0])) / 2.0)
        print(np.log10(w[-1]), np.log10(w[0]), w0)
        x, info, ier, msg = so.fsolve(lambda wi: bode(system, w=wi)[2] + 180,
                      w0, full_output=True)

        if ier == 1:
            wg = abs(x)
            Gm = -bode(system, w=wg)[1]
            if len(Gm) == 1:
                wg = wg[0]
                Gm = Gm[0]
        else:
            Gm = np.inf
            wg = None
        x, info, ier, msg = so.fsolve(lambda wi: bode(system, w=wi)[1],
                      w0, full_output=True)

        if ier == 1:
            wp = abs(x)
            Pm = bode(system, w=wp)[2] + 180.0
            if len(Pm) == 1:
                wp = wp[0]
                Pm = Pm[0]
        else:
            Pm = 180.0
            wp = None
    else:
        print("System must be a LinSys instance!")
        Gm, Pm, wg, wp = None, None, None, None

    return Gm, Pm, wg, wp

def lsim(system, U=None, T=None, X0=None, **kwargs):
    if isinstance(system, LinSys):
        t, y, x = ss.lsim2(system.linsys, U=U, T=T, X0=X0, **kwargs)
        if system.dtime != 0:
            i = np.nonzero(t > system.dtime)[0][0]
            yout = np.concatenate((np.zeros(i), y[:-i]))
        else:
            yout = y
    else:
        print("System must be a LinSys instance!")
        t, yout = None, None

    return t, yout

def pade(T, N=2):
##    if isinstance(T, LinSys):
##        Tm = T.dtime
##        num = [Tm**2/12., -Tm/2., 1]
##        den = [Tm**2/12., Tm/2., 1]
##
##        h = series(LinSys(T.linsys.num, T.linsys.den), LinSys(num, den))
##    else:
##        num = [T**2/12., -T/2., 1]
##        den = [T**2/12., T/2., 1]
##
##        h = LinSys(num, den)

    num = []
    den = []

    sign = 1
    fact_n = facto(N)
    fact_2n = facto(2*N)

    if isinstance(T, LinSys):
        Tm = T.dtime

        for i in range(0, N+1):
            p = (facto(2*N-i) * fact_n) / (fact_2n * facto(i) * facto(N-i))
            pi = p * Tm**i
            num.insert(0, sign*pi)
            den.insert(0, pi)
            sign = -sign

        h = series(LinSys(T.linsys.num, T.linsys.den), LinSys(num, den))
    elif isinstance(T, (int, float)):
        for i in range(0, N+1):
            p = (facto(2*N-i) * fact_n) / (fact_2n * facto(i) * facto(N-i))
            pi = p * T**i
            num.insert(0, sign*pi)
            den.insert(0, pi)
            sign = -sign

        h = LinSys(num, den)
    else:
        print("T has to be a Number or a LinSys instance!")

    return h

if __name__ == "__main__":
##    h1 = ss.lti(1, [1, 1])
##    h2 = ss.lti(2, [10, 1])
##
##    h3 = serial(h1, h2)
##
##    print(h3.num, h3.den)
##
##    ht = ss.lti(10, 1)
##    print(ht.num, ht.den)
##
##    H = [ht, ht, h1, ht, ht, h2]
##    hs = reduce(serial, H)
##
##    print(hs.num, hs.den)
##
##    h4 = parallel(h1, h2)
##
##    print(h4.num, h4.den)
##
##    h5 = feedback(h1, ss.lti(1,1), 1)
##
##    print(h5.num, h5.den)
##
##    signs = [1, -1, 1, -1, 1, 1]
##    hp = reduce(parallel, map(lambda h, sign: ss.lti(sign*h.num, h.den),
##                H, signs))
##
##    print(hp.num, hp.den)

    #Création des variables
    R1 = 10;
    R2 = 10000;
    L1 = 0.001;
    C1 = 100*10**-6

    #Calcul des paramètres de la fonction de transfert
    K = 1/(R2*C1)
    T1 = (R1+R2)/L1
    T2 = (R1*R2*C1+L1)/(L1*R2*C1)
    T3 = (R1+R2)/(L1*R2*C1)

    #Création de la fonction de transfert
    Num = [1,T1]
    Den = [1,T2,T3]
    H = LinSys(K*np.array(Num), Den, dt=0)

    my_lti = LinSys(2, [1, 1], dt=2)
    my_pade = pade(1)
    H1 = my_lti * my_pade * 2
    chain_f = LinSys(1, 1)
    Hf = feedback(my_lti, chain_f, pade_order=10)
    #H1 = LinSys(np.polymul([3, 1],[1.5, 1]), np.polymul(np.polymul([10, 7, 1], [4, -1]), [1, 1, 1, 0]))
    t = np.linspace(0, 10, 1001)
    t, yout = step(Hf, T=t)
    mp.figure(1)
    mp.plot(t, yout)
    mp.grid(b=True, axis='both')

    mp.figure(3)
    print(margin(my_lti))
##    #T = np.linspace(0, 10, 101)
    w, m, p = bode(my_lti)
##    #Gm, Pm, wg, wp = margin(H)
##    #print(Gm, Pm, wg, wp)
##    #print(my_lti)
##    #my_sys = LinSys(10, [10, 1], dt=2)
##    #my_lti.dtime = 0
##    #my_sys.dtime = 0
##    #hs = feedback(my_lti, my_sys)
##    #print(hs.linsys.num, hs.linsys.den, hs.dtime)
##    #t, y = ss.step2(hs.linsys, T=np.linspace(0, 10, 1001))
    mp.subplot(211)
    mp.semilogx(w, m)
    mp.subplot(212)
    mp.semilogx(w, p)
##    mp.figure(2)
##    bodeasympt(H1, 0)
    mp.show()
