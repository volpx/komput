#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA_FOLDER = "output_data/"


def dumb_derivative(x, y, yerr=None):
    dx_s = x[1:]-x[:-1]
    dy_s = y[1:]-y[:-1]
    dy_dx_y_s = dy_s/dx_s
    dy_dx_x_s = (x[1:]+x[:-1])/2

    if(yerr is None):
        return dy_dx_x_s, dy_dx_y_s
    else:
        dy_dx_y_s_err = np.sqrt(yerr[1:]**2+yerr[:-1]**2)/dx_s
        return dy_dx_x_s, dy_dx_y_s, dy_dx_y_s_err


def main():
    dfB = pd.read_csv(DATA_FOLDER+"B.dat", sep=' ',
                      names=["rho", "B", "stdB"], header=0)

    figB = plt.figure()
    figB.suptitle("B")
    ax = figB.subplots()
    ax.set_xlabel("rho")
    ax.set_ylabel("B")
    ax.grid()
    ax.errorbar(dfB["rho"], dfB["B"], yerr=dfB["stdB"], fmt='.')

    rho1_s, A_s, dA_s = dumb_derivative(
        dfB["rho"].values, dfB["B"].values, dfB["stdB"].values)

    theorA_df = pd.read_csv(DATA_FOLDER+"TheorA.dat",
                            sep='\s+', skiprows=1, header=None).values
    Ts = theorA_df[0, 1:]
    rhos = theorA_df[1:, 0]
    As = theorA_df[1:, 1:]
    # print(As)

    figdB = plt.figure()
    figdB.suptitle("A")
    ax = figdB.subplots()
    ax.set_xlabel("rho")
    ax.set_ylabel("dB/drho")
    ax.grid()
    ax.errorbar(rho1_s, A_s, yerr=dA_s, fmt='.')
    ax.errorbar(rhos, As[:, 2], label="T=1")

    plt.show()


if __name__ == "__main__":
    main()
