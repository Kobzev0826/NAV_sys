import numpy as np
import math
import os
import sys
import pandas as pd
import datetime
from scipy.optimize import fsolve


def E_k_funk(x, e, M_k):
    return x - e * np.sin(x) - M_k


def GEO_calc(e, t_oe, i, Omega_doc, A, Omega_0, w, af0, af1, M_0, year, month,day, N_var):
    date_time = datetime.datetime(year, month, day, N_var, 00, 00)
    #print(date_time.weekday())
    t = (date_time.weekday() * 24 - 3 + N_var) * 60 * 60 + 5*60


    # ---constants--------------------------------------------------------------------
    mu = 3.986004418 * 10 ** 14  # [m/s] Geocentric gravitational constant of BDCS
    Omega = 7.2921150 * 10 ** (-5)  # [rad/s] Earth’s rotation rate of BDCS
    pi = 3.1415926535898
    i_0 = pi * 0.3  # IGEO/GEO
    # ----------------------------------------------------------------------------------
    t_k = t - t_oe
    if t_k > 302400:
        t_k -= 604800
    elif t_k < -302400:
        t_k += 604800

    A = A ** 2

    n_0 = (mu / (A ** 3)) ** 0.5
    M_k = M_0 + n_0 * t_k

    E_k = fsolve(E_k_funk, x0=1, args=(e, M_k))

    v_k = np.arcsin((np.sin(E_k) * ((1 - e ** 2) ** 0.5)) / (1 - e * np.cos(E_k)))
    fi_k = v_k + w

    r_k = A * (1 - e * np.cos(E_k))  # Corrected radius

    # -----Position in orbital plane-------------------------------------
    x_k = r_k * np.cos(fi_k)
    y_k = r_k * np.sin(fi_k)
    # -----------------------------------------------------------------------
    Omega_k = Omega_0 + (Omega_doc - Omega) * t_k - Omega * t_oe  # Corrected longitude of ascending node
    i_k = i_0 + i  # Inclination at reference time

    X_k = x_k * np.cos(Omega_k) - y_k * np.cos(i_k) * np.sin(Omega_k)
    Y_k = x_k * np.sin(Omega_k) + y_k * np.cos(i_k) * np.cos(Omega_k)
    Z_k = y_k * np.sin(i_k)
    return X_k, Y_k, Z_k


def Medvedev(x_k, y_k, z_k, x, y, z):
    delta_r_x = x_k - x  # (m)

    delta_r_y = y_k - y  # (m)

    delta_r_z = z_k - z  # (m)

    a = 6378137  # semi-major axis of an ellipsoid

    LL = math.atan2(y, x)  # longotude
    ff = 298.257223  # compression

    k0 = (ff - 1) / ff
    k1 = a * (2 * ff - 1) / ff / (ff - 1)
    k2 = k0 * k1

    R = math.sqrt(x ** 2 + y ** 2)
    U = math.atan((k1 / math.sqrt(z ** 2 + (k0 * R) ** 2) + 1) * k0 * z / R)
    BB = math.atan2(z + k1 * (math.sin(U)) ** 3, R - k2 * (math.cos(U)) ** 3)  # latitude

    A1 = np.array([[delta_r_x], [delta_r_y], [delta_r_z]])

    # Move on topografic rectangular coorsinates

    A2 = np.array([[-math.sin(BB) * math.cos(LL), -math.sin(BB) * math.sin(LL), math.cos(BB)],

                   [-math.sin(LL), math.cos(LL), 0.],

                   [math.cos(BB) * math.cos(LL), math.cos(BB) * math.sin(LL), math.sin(BB)]])

    Top_rec_ccor = A2.dot(A1)

    # zenith distance:

    Z_r = math.atan2(np.sqrt(Top_rec_ccor[0] ** 2 + Top_rec_ccor[1] ** 2), Top_rec_ccor[2])

    elevation = math.pi / 2 - Z_r  # evidence (rad)
    elevation_deg = elevation * 180 / math.pi  # evidence (degree)

    azimut = math.atan2(Top_rec_ccor[1], Top_rec_ccor[0])  # azimuth (rad)
    azimut_deg = azimut * 180 / math.pi  # azimuth (grad)

    if azimut_deg < 0:
        azimut_deg = azimut_deg + 360
    elif azimut_deg > 360:
        azimut_deg = azimut_deg - 360

    return azimut_deg,elevation_deg


def pars_file(name, N_var):
    book = pd.read_excel(name, header=1)

    data = book.loc[book['PRN'] == 'C' + f"{N_var:02d}"]
    return data


if __name__ == "__main__":
    #name = input('input name xl file with exception: ')#30_09_2021.xls'
    name = '30_09_2021.xls'
    #N_var = int(input('N_var(Time(H)): '))
    #N_satellite = int(input('satellite: '))  # N_var

    dir = os.path.abspath(os.curdir)

    #print(data)

    X_obser = 2846228.896
    Y_obser = 2198658.103
    Z_obser = 5249983.343

    lat_observer = 55.76581124
    long_observer = 37.68540894

    #choose = 'y'
    while True:
        k=0

        N_var = int(input('N_var(Time(H)): '))
        N_satellite = int(input('satellite: '))  # N_var

        data = pars_file(str(dir) + '\\' + str(name), N_satellite)

        X, Y, Z = GEO_calc(e=data.iloc[0]['e'], t_oe=data.iloc[0]['t'], i=data.iloc[0]['δi'], Omega_doc=data.iloc[0]['Ω'],
                           A=data.iloc[0]['A'], Omega_0=data.iloc[0]['Ω0'], w=data.iloc[0]['ω'], M_0=data.iloc[0]['m'],
                           af0=data.iloc[0]['af0'], af1=data.iloc[0]['af1'], year=2021, month=9, day=30, N_var=N_var)

        print("coordinats=",X, Y, Z)

        azimut, elevation = Medvedev(X[0],Y[0], Z[0], X_obser, Y_obser, Z_obser)
        print("AZ: ", azimut,'\n'+' EL: ', elevation)

        while k==0:
            choose = input('Continue?[y/n]' )
            if choose == 'n':
                sys.exit()
            elif choose != 'y' and choose != 'n':
                print("the answer could not be recognized, please try again")
                k=0
            else:
                k=1