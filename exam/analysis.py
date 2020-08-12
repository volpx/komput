#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

DATA_FOLDER = "data/"
# BFILE = DATA_FOLDER+"B.dat"
# BFILE = DATA_FOLDER+"B_MC.dat"
# BFILE = DATA_FOLDER+"bak/B_lotsofpoints_good_thermal5000.dat"
# BFILE = DATA_FOLDER+"bak/B_nose_final.dat"
BFILE = DATA_FOLDER+"bak/B_final2.dat"


def main():
    # %% Get B
    Ts, rhos, B, stdB = getB(BFILE)
    start, stop = 0, rhos.size
    rhos = rhos[start:stop]
    B = B[:, start:stop]
    stdB = stdB[:, start:stop]
    print("Points:", rhos.size)

    stdB[0, :] *= 2
    stdB[1:, :] *= 1.4

    # %% Work on B
    # # Smooth out the curve
    n = 10
    B1 = np.empty((Ts.size, rhos.size-n+1))
    stdB1 = np.empty_like(B1)
    for i in range(Ts.size):
        rhos1, B1[i, :], stdB1[i, :] = flexible_mean_smoother(
            rhos, B[i, :], n, True)

    # fit
    npar = 4
    parameters = np.empty((Ts.size, npar))
    Dparameters = np.empty((Ts.size, npar, npar))
    Bfit = np.empty_like(B)
    stdBfit = np.empty_like(Bfit)
    for i in range(Ts.size):
        if i == 0:
            parameters[i, :],  Dparameters[i, :, :] = curve_fit(
                ffit, rhos[:23], B[i, :23], np.ones((npar,)), stdB[i, :23],
                method="lm", maxfev=5000)
            Bfit[i, :] = ffit(rhos, *parameters[i, :])
            stdBfit[i, :] = ffit_ydev(
                rhos, parameters[i, :],  Dparameters[i, :, :])
            chi = chi2red(
                B[i, :23], stdB[i, :23], Bfit[i, :23], npar)
        else:
            parameters[i, :],  Dparameters[i, :, :] = curve_fit(
                ffit, rhos[:], B[i, :], np.ones((npar,)), stdB[i, :],
                method="lm", maxfev=5000)
            Bfit[i, :] = ffit(rhos, *parameters[i, :])
            stdBfit[i, :] = ffit_ydev(
                rhos, parameters[i, :],  Dparameters[i, :, :])
            chi = chi2red(
                B[i, :], stdB[i, :], Bfit[i, :], npar)
        print("Chi2red:", chi, "T:", Ts[i])

    # fit by pieces
    # npar = 3
    # pieces_length = 10
    # rho_length = rhos.size
    # pieces = rho_length//pieces_length
    # Bfit2 = np.empty_like(B)
    # stdBfit2 = np.empty_like(Bfit2)
    # parameters2 = np.empty((Ts.size, pieces, npar))
    # Dparameters2 = np.empty((Ts.size, pieces, npar, npar))
    # for i in range(Ts.size):
    #     for j in range(pieces):
    #         parameters2[i, j, :],  Dparameters2[i, j, :, :] = curve_fit(
    #             ffit, rhos[j*pieces_length:(j+1)*pieces_length],
    #             B[i, j*pieces_length:(j+1)*pieces_length],
    #             np.ones((npar,)), stdB[i, j*pieces_length:(j+1)*pieces_length],
    #             method="lm", maxfev=5000)
    #         Bfit2[i, j*pieces_length:(j+1) *
    #               pieces_length] = ffit(
    #                   rhos[j*pieces_length:(j+1)*pieces_length],
    #                   *parameters2[i, j, :])
    #         stdBfit2[i, j*pieces_length:(j+1)*pieces_length] = ffit_ydev(
    #             rhos[j*pieces_length:(j+1)*pieces_length],
    #             parameters2[i, j, :],  Dparameters2[i, j, :, :])
    #     print("Chi2red:", chi2red(
    #         B[i, :], stdB[i, :], Bfit2[i, :], 3*pieces), "T:", Ts[i])

    figB = plt.figure()
    figB.suptitle("B")
    ax = figB.subplots()
    ax.set_xlabel("rho")
    ax.set_ylabel("B")
    ax.grid()
    for i, T in enumerate(Ts):
        ax.errorbar(rhos, B[i, :],
                    yerr=stdB[i, :], fmt='.-',
                    label="T: {}".format(T))
        ax.errorbar(rhos,
                    Bfit[i, :],
                    yerr=stdBfit[i, :],
                    label="T:{} fit".format(T))
        # ax.errorbar(rhos,
        #             Bfit2[i, :],
        #             yerr=stdBfit2[i, :],
        #             label="T:{} mfit".format(T))
        # ax.errorbar(rhos1, B1[i, :],
        #             yerr=stdB1[i, :], fmt='.',
        #             label="T: {} smoothed".format(T))
    ax.legend()
    # figB.savefig("data/B.pdf")

    # %% Compute A
    As = np.empty((B.shape[0], B.shape[1]-1))
    stdAs = np.empty_like(As)
    for i in range(Ts.size):
        rhods, As[i, :], stdAs[i, :] = dumb_derivative(
            rhos, B[i, :], stdB[i, :])

    As_fit = np.empty_like(Bfit)
    stdAs_fit = np.empty_like(As_fit)
    for i in range(Ts.size):
        As_fit[i, :] = dffit(rhos, *parameters[i, :])
        stdAs_fit[i, :] = dffit_ydev(
            rhos, parameters[i, :],  Dparameters[i, :])

    # As_fit2 = np.empty_like(Bfit2)
    # stdAs_fit2 = np.empty_like(As_fit2)
    # for i in range(Ts.size):
    #     for j in range(pieces):
    #         As_fit2[i, j*pieces_length:(j+1)*pieces_length] = dffit(
    #             rhos[j*pieces_length:(j+1)*pieces_length], *parameters2[i, j, :])
    #         stdAs_fit2[i, j*pieces_length:(j+1)*pieces_length] = dffit_ydev(
    #             rhos[j*pieces_length:(j+1)*pieces_length],
    #             parameters2[i, j, :],  Dparameters2[i, j, :])

    # As_mp = np.empty_like(Bfit)
    # stdAs_mp = np.empty_like(As_mp)
    # for i in range(Ts.size):
    #     As_mp[i, :], stdAs_mp[i, :] = derivative_5points(
    #         rhos, B[i, :], stdB[i, :])

    # # On smooothed
    As1 = np.empty((B1.shape[0], B1.shape[1]-1))
    stdAs1 = np.empty_like(As1)
    As11 = np.empty((As1.shape[0], As1.shape[1]-n+1))
    stdAs11 = np.empty_like(As11)
    As_smmp = np.empty_like(B1)
    stdAs_smmp = np.empty_like(As_smmp)
    for i in range(Ts.size):
        rhods1, As1[i, :], stdAs1[i, :] = dumb_derivative(
            rhos1, B1[i, :], stdB1[i, :])
        rhods11, As11[i, :], stdAs11[i, :] = flexible_mean_smoother(
            rhods1, As1[i, :], n, True)
        As_smmp[i, :], stdAs_smmp[i, :] = derivative_5points(
            rhos1, B1[i, :], stdB1[i, :])

    theorA_df = pd.read_csv(DATA_FOLDER+"TheorA.dat",
                            sep='\s+', skiprows=1, header=None).values
    Ts_th = theorA_df[0, 1:]
    rhos_th = theorA_df[1:, 0]
    As_th = theorA_df[1:, 1:]

    # %% Work on A
    figdA = plt.figure()
    figdA.suptitle("A")
    ax = figdA.subplots()
    ax.set_xlabel("rho")
    ax.set_ylabel("A")
    ax.grid()
    ax.set_xlim([0.64, 0.91])
    ax.set_ylim([20, 90])
    # Theorical from paper
    for i, T in enumerate(Ts_th):
        ax.errorbar(rhos_th, As_th[i, :],
                    fmt='.-',
                    label="th T: {}".format(T))
    # # From original
    # for i, T in enumerate(Ts):
    #     ax.errorbar(rhods, As[i, :],
    #                 yerr=stdAs[i, :], fmt='.-',
    #                 label="T: {}".format(T))
    # From fitted
    for i, T in enumerate(Ts):
        ax.errorbar(rhos, As_fit[i, :],
                    yerr=stdAs_fit[i, :], fmt='.-',
                    label="fit T: {}".format(T))
    # From fitted 2
    # for i, T in enumerate(Ts):
    #     ax.errorbar(rhos, As_fit2[i, :],
    #                 yerr=stdAs_fit2[i, :], fmt='.-',
    #                 label="fitm T: {}".format(T))
    # # From smoothed
    # for i, T in enumerate(Ts):
    #     ax.errorbar(rhods1, As1[i, :],
    #                 yerr=stdAs1[i, :],
    #                 fmt='.-',
    #                 label="sm T: {}".format(T))
    #     # Double
    #     ax.errorbar(rhods11, As11[i, :],
    #                 yerr=stdAs11[i, :],
    #                 fmt='.-',
    #                 label="smsm T: {}".format(T))
    # From multipoint der
    # for i, T in enumerate(Ts):
    #     ax.errorbar(rhos, As_mp[i, :],
    #                 yerr=stdAs_mp[i, :], fmt='.-',
    #                 label="mp T: {}".format(T))
    # From multipoint der on smoothed
    # for i, T in enumerate(Ts):
    #     ax.errorbar(rhos1, As_smmp[i, :],
    #                 yerr=stdAs_smmp[i, :], fmt='.-',
    #                 label="smmp T: {}".format(T))
    ax.legend(loc="upper left")
    figdA.savefig("data/A.pdf", dpi=1000)

    # plt.show()
    # return

    # %% P corrections

    def p1(rho, T, A):
        return 2*rho/T**2*A

    p1s = np.empty_like(As_fit)
    stdp1s = np.empty_like(p1s)
    for i, T in enumerate(Ts):
        p1s[i, :] = p1(rhos, T, As_fit[i, :])
        stdp1s[i, :] = p1(rhos, T, stdAs_fit[i, :])

    figdP1 = plt.figure()
    figdP1.suptitle("first quantum corrections to the compressibility factor")
    ax = figdP1.subplots()
    ax.set_xlabel("rho")
    ax.set_ylabel("p1/nT")
    ax.grid()
    # From fitted
    for i, T in enumerate(Ts):
        ax.errorbar(rhos, p1s[i, :]/rhos/T,
                    yerr=stdp1s[i, :]/rhos/T, fmt='.-',
                    label="fit T: {}".format(T))
    # Theorical from paper
    # for i, T in enumerate(Ts_th):
    #     ax.errorbar(rhos_th, As_th[i, :],
    #                 fmt='.-',
    #                 label="th T: {}".format(T))
    ax.legend()

    def p(rho, T, A):
        return

    plt.show()


def chi2red(data, dy, model=0, ddof=0):
    return np.sum(((model - data)**2)/(dy**2)) / (data.size - ddof)


def flexible_mean_smoother(datax, datay, n, datayerr=None):
    x = np.empty((datax.size-n+1,))
    y = np.empty_like(x)
    if datayerr is None:
        for i in range(y.size):
            x[i] = np.mean(datax[i:i+n])
            y[i] = np.mean(datay[i:i+n])
        return x, y
    else:
        yerr = np.empty_like(y)
        for i in range(y.size):
            x[i] = np.mean(datax[i:i+n])
            y[i] = np.mean(datay[i:i+n])
            yerr[i] = np.std(datay[i:i+n], ddof=1)/np.sqrt(n)
        return x, y, yerr


def dumb_derivative(x, y, yerr=None):
    dx_s = x[1:]-x[:-1]
    dy_s = y[1:]-y[:-1]
    dy_dx_y_s = dy_s/dx_s
    dy_dx_x_s = (x[1:]+x[:-1])/2

    if yerr is None:
        return dy_dx_x_s, dy_dx_y_s
    else:
        dy_dx_y_s_err = np.sqrt(yerr[1:]**2+yerr[:-1]**2)/dx_s
        return dy_dx_x_s, dy_dx_y_s, dy_dx_y_s_err


def derivative_5points(x, y, yerr=None):
    # x must be uniformly spaced and of length >= 4
    h = x[1]-x[0]
    y_x = np.empty_like(y)
    y_x[0] = (-3*y[0+0]+4*y[0+1]-1*y[0+2])/(2*1.0*h**1)
    y_x[1] = (-2*y[1-1]-3*y[1+0]+6*y[1+1]-1*y[1+2])/(6*1.0*h**1)
    y_x[-1] = (1*y[-1-2]-4*y[-1-1]+3*y[-1+0])/(2*1.0*h**1)
    y_x[-2] = (1*y[-2-2]-6*y[-2-1]+3*y[-2+0]+2*y[-2+1])/(6*1.0*h**1)
    for i in range(2, y_x.size-2):
        y_x[i] = (1*y[i-2]-8*y[i-1]+0*y[i+0]+8 *
                  y[i+1]-1*y[i+2])/(12*1.0*h**1)

    if yerr is None:
        return y_x
    else:
        stdy_x = np.empty_like(y_x)

        stdy_x[0] = np.sqrt((-3*yerr[0+0])**2+(4*yerr[0+1])
                            ** 2+(1*yerr[0+2])**2)/(2*1.0*h**1)
        stdy_x[1] = np.sqrt((-2*yerr[1-1])**2+(3*yerr[1+0])**2+(6 *
                                                                yerr[1+1])**2+(1*yerr[1+2])**2)/(6*1.0*h**1)
        stdy_x[-1] = np.sqrt((1*yerr[-1-2])**2+(4*yerr[-1-1])
                             ** 2+(3*yerr[-1+0])**2)/(2*1.0*h**1)
        stdy_x[-2] = np.sqrt((1*yerr[-2-2])**2+(6*yerr[-2-1])**2+(3 *
                                                                  yerr[-2+0])**2+(2*yerr[-2+1])**2)/(6*1.0*h**1)
        for i in range(2, stdy_x.size-2):
            stdy_x[i] = np.sqrt((1*yerr[i-2])**2+(8*yerr[i-1])**2+(0*yerr[i+0])**2+(8 *
                                                                                    yerr[i+1])**2+(1*yerr[i+2])**2)/(12*1.0*h**1)
        return y_x, stdy_x


def get_different_values(a):
    res = []
    for el in a:
        if not el in res:
            res.append(el)
    return np.array(res)


def getB(filen):
    dfB = pd.read_csv(filen, sep=' ',
                      names=["T", "rho", "B", "stdB"], header=0)

    Ts = get_different_values(dfB["T"].values)
    rhos = get_different_values(dfB["rho"].values)

    B = dfB["B"].values.reshape((Ts.size, rhos.size))
    stdB = dfB["stdB"].values.reshape((Ts.size, rhos.size))

    return Ts, rhos, B, stdB

# # FIT FUNCTIONS poly
# def ffit(x, a, b, c, d, e):
#     return a+b*x+c*x**2+d*x**3+e*x**4


# def Jacobian_parameters_ffit(x, a, b, c, d, e):
#     return np.array([x**0, x, x**2, x**3, x**4]).T


# def dffit(x, a, b, c, d, e):
#     return b+2*c*x+3*d*x**2+4*e*x**3


# def Jacobian_parameters_dffit(x, a, b, c, d, e):
#     return np.array([np.zeros_like(x), np.ones_like(x),
#                      2*x, 3*x**2, 4*x**3]).T


# FIT FUNCTIONS exp
# def ffit(x, a, b, c):
#     return a*np.exp(b*x)+c


# def Jacobian_parameters_ffit(x, a, b, c):
#     return np.array([np.exp(b*x), a*np.exp(b*x)*x, np.ones_like(x)]).T


# def dffit(x, a, b, c):
#     return a*b*np.exp(b*x)


# def Jacobian_parameters_dffit(x, a, b, c):
#     return np.array([b*np.exp(b*x),
#                      a*np.exp(b*x)+a*b*np.exp(b*x)*x,
#                      np.zeros_like(x)]).T

# FIT FUNCTIONS 2nd degree
# def ffit(x, a, b, c):
#     return a+b*x+c*x**2


# def Jacobian_parameters_ffit(x, a, b, c):
#     return np.array([np.ones_like(x), x, x**2]).T


# def dffit(x, a, b, c):
#     return b+2*c*x


# def Jacobian_parameters_dffit(x, a, b, c):
#     return np.array([np.zeros_like(x),
#                      np.ones_like(x),
#                      2*x]).T


# FIT FUNCTIONS 3rd degree
def ffit(x, a, b, c, d):
    return a+b*x+c*x**2+d*x**3


def Jacobian_parameters_ffit(x, a, b, c, d):
    return np.array([np.ones_like(x), x, x**2, x**3]).T


def dffit(x, a, b, c, d):
    return b+2*c*x+3*d*x**2


def Jacobian_parameters_dffit(x, a, b, c, d):
    return np.array([np.zeros_like(x),
                     np.ones_like(x),
                     2*x, 3*x**2]).T


def dffit_ydev(x, popt, pcov):
    jac = Jacobian_parameters_dffit(x, *popt)
    return np.diag(jac@pcov@jac.T)


def ffit_ydev(x, popt, pcov):
    jac = Jacobian_parameters_ffit(x, *popt)
    return np.diag(jac@pcov@jac.T)


if __name__ == "__main__":
    main()
