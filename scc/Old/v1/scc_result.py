from scipy.constants import mu_0 as MUO, epsilon_0 as EO

class Result:
    @staticmethod
    def show_results(scc, syst, Gri, zG, zI, A, ZL, Zs):
        print("PHYSICAL CONSTANTS")
        print(f"  Vacuum magnetic permeability. mu_0 = {MUO:.4e} H/m")
        print(f"  Vacuum electrical permittivity. epsilon_0 = {EO:.4e} F/m")

        print("\nSINGLE CORE CABLES PARAMETERS")
        for param, value in vars(scc).items():
            print(f"  {param} = {value} m")

        print("\nSOIL PARAMETERS")
        print(f"  Ground resistivity: rho1 = {syst.rho1} ohm.m")
        print(f"  Ground relative permittivity: eps_r1 = {syst.eps_r1}")

        print("\nSCC INTERNAL IMPEDANCE")
        for param, value in vars(zI).items():
            if value != 0:
                print(f"  {param} = {value:.4e} ohm/m")

        print("\nMATRIX DISTANCES")
        print("  \nd =\n", Gri.d, "m")
        print("  \nD =\n", Gri.D, "m")

        print("\nGROUND RETURN IMPEDANCE")
        print(f"  Air Propagation constant. y0 = {Gri.y0:.4e}")
        print(f"  Ground Propagation constant. y1 = {Gri.y1:.4e}")

        print("\n  De Conti et al. (2023) Closed-Form Expressions")
        print(f"  Ground return impedance Matrix.\n zG = \n {zG} ohm/m")

        print("\nSERIES IMPEDANCE MATRIX (NODA,2008)")
        print(f"  Transformation Matrix.\n A = \n {A}\n")
        print(f"  Loop impedance.\n ZL = \n {ZL} ohm/m\n")
        print(f"  Series Impedance.\n Zs = \n {Zs} ohm/m\n")
