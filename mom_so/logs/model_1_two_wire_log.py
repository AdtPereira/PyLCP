import datetime
import numpy as np
from scipy.constants import mu_0 as MU0, epsilon_0 as E0

class MainParameters():
    def __init__(self):
        self.file_path = 'C:\\Users\\adilt\\OneDrive\\01 ACADEMIA\\06 MODELS\\7.MoM-SO\\logs\\'

    def my_complex_form(self, decimals):
        def eng_format_complex(x):
            def eng_format(x):
                magnitude = 0
                while abs(x) >= 1000:
                    magnitude += 1
                    x /= 1000.0
                return '{:.{}f}e{}'.format(x, decimals, magnitude*3)

            real = eng_format(x.real)
            imag = eng_format(x.imag)

            if x.imag < 0:
                complex = f'{real}{imag}j'
            elif x.imag > 0:
                complex = f'{real}+{imag}j'
            elif x.imag == 0:
                complex = f'{real}'

            return complex

        return eng_format_complex


    def principal(self, system_name, geo, tl, mom_so, green):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        filename = f"{self.file_path}{system_name}_{timestamp}.log"

        # Writing the log file
        with open(filename, 'w') as Log_File:
            
            # Writing the header
            Log_File.write(                               
                f"==================================================================================\n"
                f"Universidade Federal de Minas Gerais (UFMG)\n"
                f"Programa de Pos-Graduacao em Engenharia Eletrica (PPGEE)\n\n"                
                
                f"MoM-SO-v1\n"
                f"Two-Wire Transmission Line Analysis using the \n"
                f"Method of Moments - Superficial Admittance Operator (MoM-SO)\n\n"

                f"Author: Adilton Pereira (c)\n"
                f"Belo Horizonte, {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                f"Brazil"

                f"\n\n"
                f"References:\n"
                f"[1]   PATEL, Utkarsh R. A Surface Admittance Approach For Fast Calculation of the \n" 
                f"      Series Impedance of Cables Including Skin, Proximity, and Ground Return     \n"
                f"      Effects. 2014. University of Toronto, Graduate Department of The Edward S.  \n" 
                f"      Rogers Sr. Department of Electrical & Computer Engineering.                 \n"
                f"[2]   U. R. Patel, B. Gustavsen and P. Triverio, 'An Equivalent Surface Current   \n"
                f"      Approach Computation of the Series Impedance of Power Cables with Inclusion \n"
                f"      of Skin and Proximity Effects,' IEEE Transactions on Power Delivery, vol.28,\n"
                f"      no. 4, pp. 2474-2482, Oct. 2013, doi: 10.1109/TPWRD.2013.2267098            \n\n" 
                f"==================================================================================\n"
            )

            # Writing the main parameters
            Log_File.write(
                f"\nUNIVERSAL PHYSICAL CONSTANTS.....................................\n"
                f" Vacuum magnetic permeability. mu_0 = {MU0:.4e} H/m.\n"
                f" Vacuum electrical permittivity. epsilon_0 = {E0:.4e} F/m.\n" 


                f"\nGEOMETRY OF THE SYSTEM...........................................\n"
                f" Number of points in conductors surface: {geo.surface_points}.\n"
                f" Distance matrix between the center points of the conductors [m]:\n"
                f" {geo.distance_matrices()[0]}\n"
                f" x-coord. matrix between the center points of the conductors [m]:\n"
                f" {geo.distance_matrices()[1]}\n"
                f" y-coord. matrix between the center points of the conductors [m]:\n"
                f" {geo.distance_matrices()[2]}\n"
                f" theta angle between the center points of the conductors [rad]:\n"
                f" {geo.distance_matrices()[3]}\n"
            )

            # Writing the transmission line parameters
            Log_File.write(
                f"\nFREQUENCY-DEPENDENT PARAMETERS...................................\n"
                f"Conductivity of the conductors: "  
                f"{np.array2string(
                    tl.sigma, precision=2, formatter={'float_kind':'{:.2e}'.format})} [S/m]\n"
                f"Permittivity of the medium: {tl.epsilon}.\n"
                f"Permeability of the medium: {tl.mu}.\n"
                f"Conductors wavenumber: "
                f"{np.array2string(
                    tl.k, precision=2, formatter={'float_kind':'{:.2e}'.format})} [1/m]\n"
                f"Free-space wavenumber: "
                f"{np.array2string(
                    tl.kout, precision=2, formatter={'float_kind':'{:.2e}'.format})} [1/m]\n"
            )

            # Writing the MoM-SO parameters
            Log_File.write(
                f"\nMoM-SO INITIALIZATION............................................\n"
                f" Number of elements in the MoM vectors: N = {mom_so.big_n}.\n"
                f" Matrix U:\n"
                f" {mom_so.matrix_u()}\n\n"
                f" Surface admittance operator: Ys:\n" 
                f" {np.array2string(
                    mom_so.matrix_ys(), 
                    formatter={'complex_kind': self.my_complex_form(2)})} \n\n"
            )

            # Writing the Green's elements
            Log_File.write(
                f"\nELEMENTS OF GREENS SUB-MATRICES...............................\n\n"

                f" 1st Case: Self-Terms (p = q = 0)\n"
                f"  Matrix G_(n' =n, 0)^(p,p): {green.ieee_paper(0, 0, 0, 0)}\n"
                f"  Matrix G_(n' =n, n)^(p,p): {green.ieee_paper(1, 1, 0, 0)}\n"
                f"  Matrix G_(n'!=n, n)^(p,p): {green.ieee_paper(-1, 1, 0, 0)}\n\n"

                f" 2nd Case: Mutual-Terms (p != q)\n"
                f"  1st Sub-Case: Central Column [n=0] of G_(n'n)^(p,q)\n"
                f"  Matrix G_( n' =n, 0)^(p,q): {green.ieee_paper( 0, 0, 0, 1)}\n"
                f"  Matrix G_(-n'!=n, 0)^(p,q): {green.ieee_paper(-1, 0, 0, 1)}\n"
                f"  Matrix G_(+n'!=n, 0)^(p,q): {green.ieee_paper( 1, 0, 0, 1)}\n\n"
            )

            # Writing the Green sub-matrices
            Log_File.write(
                f"GREEN SUB-MATRIX, G_(n',n)^(0, 0).......................................\n"
                f" {green.sub_matrices()[0]} \n\n"    

                f"GREEN SUB-MATRIX, G_(n',n)^(0, 1).......................................\n"
                f" {green.sub_matrices()[1]} \n\n" 

                f"GREEN SUB-MATRIX, G_(n',n)^(1, 0).......................................\n"
                f" {green.sub_matrices()[2]} \n\n" 

                f"GREEN SUB-MATRIX, G_(n',n)^(1, 1).......................................\n"
                f" {green.sub_matrices()[3]} \n\n"  
            )

            # Writing the full Green matrix
            Log_File.write(
                f"\nFULL GREEN MATRIX, G .................................................\n"
                f" {mom_so.matrix_g()} \n\n"    

            )

            # Writing the impedance matrix
            Log_File.write(
                f"\nMATRIX IMPEDANCE, Z............................................\n"
                f" {mom_so.matrix_z()}\n\n"           
            )
