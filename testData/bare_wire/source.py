import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, Any

from utils.case_utils import *
from mtl_main.graphics import IsolatedMTLRepresentation
from mtl_main.source import MulticonductorTransmissionLine
from mtl_paul.py_fortran import FortranRunner
from analytical_forms.isolated_wires import WiresHomogeneousMedia
from mom_so.quasi_static_green import QuasiStatic
from mom_so.lossless_medium import HomogeneousLosslessMedium, LosslessPostProcessing
from mom.bare_wire_systems import BareWireMoMSolver
from plotter.mom_models import MoMVisualizer

class BifilarBareWirePULParameters:
    """
    The ConvergenceAnalyzer class is a tool designed to perform and visualize a convergence analysis 
    for the electrical parameters of multiconductor transmission lines (MTLs). 
    
    It systematically runs simulations with an increasing number of Fourier coefficients to observe 
    how the calculated inductance and capacitance values stabilize. The class facilitates a comparison
    between a Fortran-based simulation (RIBBON.FOR) and a Python-based Method of Moments (MoM) implementation.
    """
    def __init__(self, project_root: Path, case_name: str, mtl: Dict[str, Any], SUM_MAX: int = 10, comsol_data: Dict[str, pd.DataFrame] = {}):
        """
        Initializes the convergence analyzer.

        Args:
            mtl_config (Dict[str, Any]): Dictionary with the MTL model configuration.
            nf_max (int): Maximum number of coefficients / harmonic order to test.
        """

        assert len([key for key in mtl.keys() if isinstance(key, int)]) == 2, "The bifilar line must contain exactly two conductors."

        self.project_root = project_root
        self.case_name = case_name
        self.comsol_data = comsol_data
        self.mtl_copy = copy.deepcopy(mtl)
        self.freq_range = {'ana': np.logspace(0, 6, num=200), 'mom': np.logspace(0, 6, num=30)}
        self.srw_ratios = {'ana': np.linspace(2.1, 8, num=300), 'mom': np.linspace(2.1, 8, num=20)}
        self.N = len([key for key in mtl.keys() if isinstance(key, int)])
        
        self.sum_max = SUM_MAX
        self.results_df = None

        # Extract parameters and prepare the Fortran runner
        self._bifilar_analytical_solution()

        # Data parameters
        self.srw_data = {}
        self.srw_mum_data = {}
        self.analytical_data = {}
        self.mom_collocation_data = {}
        self.mom_galerkin_data = {}
        self.mom_so_data = {}
        self.ribbon_data = {}

        # Additional parameters
        self.c_factor = 1e12  # F/m to nF/km
        self.l_factor = 1e6   # H/m to mH/km
        self.r_factor = 1e3   # Ohm/m to Ohm/km
        self.figsize = (12, 5)
        self.pt1 = 63
        self.pt2 = 10

        # Plotting parameters
        self.plot_params = {
            'linestyles': [':', '-.', '--', '-', ':', '-.', '--'],
            'markers': ['o', 's', '^', 'd', 'v', '<', '>'],
            'colors': ['black', 'gray', 'lightgray', 'darkgray', 'dimgray', 'silver', 'gainsboro']
        }

        # Assumes the script is run from the project's root directory.
        self.results_dir = os.path.join('testData', self.case_name, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)
        
    def _prepare_fortran_runner(self, mtl):
        """Prepares the parameters and the runner for the Fortran simulation."""

        # Get the 'sheath' dictionary safely.
        #    If 'sheath' does not exist or is None, use an empty dictionary {} as a fallback.
        refIdx = mtl['idx_ref_conductor']
        sheath_dict = mtl[refIdx].get('sheath') or {}

        self.fortran_base_params = {
            'N':    len([key for key in mtl.keys() if isinstance(key, int)]),
            'NF':   mtl[refIdx]['fourier_order'] + 1,
            'IREF': refIdx,
            'RW':   mtl[refIdx]['radius'][1],
            'TD':   sheath_dict.get('thickness', 0.0),
            'ER':   sheath_dict.get('relative_permittivity', 1.0),
            'S':    np.linalg.norm(np.array(mtl[0]['center_point']) - np.array(mtl[1]['center_point'])),
        }

        fortran_exe_path = self.project_root / 'mtl_paul' / 'RIBBON' / 'RIBBON.EXE'
        self.runner = FortranRunner(exe_path=str(fortran_exe_path), silent=True)

    def _extract_matrix_element(self, column_name: str, row: int, col: int) -> pd.Series:
        """
        Extracts and processes a specific element from a matrix column in the results DataFrame.
        The conversion factor (for inductance or capacitance) is determined automatically
        based on the column name.
        """
        sign = 1.0

        # Decide which conversion factor to use based on the column name
        if column_name.startswith('L'):
            factor = self.l_factor
        elif column_name.startswith('C'):
            factor = self.c_factor
            if row != col:
                sign = -1.0
        else:
            # Raise an error if the column is not Inductance ('L') or Capacitance ('C')
            raise ValueError(f"Could not determine the conversion factor for the column: '{column_name}'")
        
        extractor = lambda matrix: (
            sign * matrix[row, col] * factor
            if isinstance(matrix, np.ndarray)
            else np.nan
        )
        return self.results_df[column_name].apply(extractor)

    def _bifilar_analytical_solution(self):
        """Computes the analytical solution for bare wires as a reference."""
        R = self.mtl_copy[0]['radius'][1]
        D = np.linalg.norm(np.array(self.mtl_copy[0]['center_point']) - np.array(self.mtl_copy[1]['center_point']))
        self.DR_ratio = D/R
        from scipy.constants import epsilon_0
        self.analytical_bifilar_capacitance = (np.pi * epsilon_0) / np.arccosh(0.5*self.DR_ratio)

    def _configure_plot_appearance(self, ax, ylabel, data_to_plot, yscale='log'):
        """
        Helper function to configure a single subplot.

        Args:
            ax (matplotlib.axes.Axes): The subplot axis to configure.
            title (str): Subplot title.
            ylabel (str): Y-axis label.
            data_to_plot (dict): Main data for plotting.
            ref_data (tuple, optional): Reference data for plotting.
        """
        # Iterate over the data for plotting
        for key, data in data_to_plot.items():
            freq, value = data['data']
            label = data['label']
            if freq is not None and value is not None:
                if key == 'mom-so':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='o', markersize=3, zorder=2)
                elif key == 'mom':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='s', markersize=10, fillstyle='none', markeredgecolor='k', zorder=2)
                elif key == 'ribbon':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='o', markersize=8, fillstyle='none', markeredgecolor='k', zorder=2)
                elif key == 'analytical':
                    ax.plot(freq, value, label=label, color='k', linestyle='--', linewidth=1.0)
                elif key == 'wires':
                    ax.plot(freq, value, label=label, color='k', linestyle='-.', linewidth=1.0)
                elif key == 'bifilar':
                    ax.plot(freq, value, label=label, color='r', linestyle=':', linewidth=1.0)
                elif key == 'exactly':
                    ax.plot(freq, value, label=label, color='k', linestyle='-', linewidth=1.0)
                elif key == 'approx':
                    ax.plot(freq, value, label=label, color='k', linestyle=':', linewidth=1.0)

        ax.set_xscale('log')
        ax.set_yscale(yscale)
        ax.set_xlim(1E0, 1E6)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(False)

    def show_header(self):
        """Displays the script header."""
        print("\n")
        print("="*self.pt2 + " BIFILAR BARE-WIRE RIBBON CABLE SIMULATION " + "="*self.pt2)
        print(f"Project: {self.project_root}")
        print(f"Model: {self.mtl_copy['name']}")
        print(f"D/R = {self.DR_ratio:.3f}. Fourier Order (k) = {self.mtl_copy[0]['fourier_order']}.")
        print(f"RIBBON Fourier Coef./cond. (NF) = {self.mtl_copy[0]['fourier_order']+1}.")
        print(f"PYTHON Fourier Coef./cond. (NF) = {2*self.mtl_copy[0]['fourier_order']+1}.")
        print(f"Exact Bifilar Bare Wire Capacitance: {self.analytical_bifilar_capacitance * 1E12:.4f} pF/m")
        print("="*self.pt1)

    def run_single_fortran(self):
        """Runs a single simulation for a specific value of k."""
        self._prepare_fortran_runner(self.mtl_copy)
        self.runner.run_fortran(self.fortran_base_params)

        print("\n")
        print("="*self.pt2 + "                 RIBBON.FOR                " + "="*self.pt2)
        if self.runner.CAP_matrix is not None:
            if self.runner.A_matrix.shape[0] < 6:
                matrix_viewer(self.runner.A_matrix,     "Block Matrix A")
                matrix_viewer(self.runner.B_matrix,     "Block Matrix B")
                matrix_viewer(self.runner.C_matrix,     "Block Matrix C")
                matrix_viewer(self.runner.D_matrix,     "Block Matrix D")
            else:
                print(f"\nBlock matrices Shape: {self.runner.A_matrix.shape}.")

            matrix_viewer(self.runner.IND_matrix,   "External Inductance Matrix, Le (H/m)")
            matrix_viewer(self.runner.CAP_matrix,   "Capacitance Matrix, C (F/m)")
            matrix_viewer(self.runner.CAP0_matrix,  "Capacitance Matrix in Vacuum, C0 (F/m)")
            matrix_viewer(self.runner.CGEN0_matrix, "Free-Space Generalized Capacitance Matrix, CGEN0 (F/m)")
            matrix_viewer(self.runner.CGEN_matrix,  "Generalized Capacitance Matrix, CGEN (F/m)")
        else:
            print("\nNo results found.")

    def run_fortran(self):
        """Runs a single simulation for a specific value of k."""
        self._prepare_fortran_runner(self.mtl_copy)
        self.runner.run_fortran(self.fortran_base_params)

        self.ribbon_data = {
            freq: {'c': self.c_factor * self.runner.CAP0_matrix.item(), 
                   'le': self.l_factor * self.runner.IND_matrix.item()} for freq in self.freq_range['mom']}

    def run_mom_methods(self, case_name, autoPlots=False):
        """
        Runs the classical Method of Moments (MoM) simulation for the bifilar transmission line.

        Returns:
            BifilarMoM: Instance of the configured BifilarMoM object.
        """
        print("\n============== pyMoM MulticonductorBareWireSystems =============")

        mtl_model = MulticonductorTransmissionLine(self.mtl_copy)
        # IsolatedMTLRepresentation(mtl_model, self.case_name, units='millimeter').system_schematic()
        solver = BareWireMoMSolver(mtl_model)
        solver.run_collocation_method()
        solver.run_galerkin_method()
        
        visualizer = MoMVisualizer(solver, case_name)
        visualizer.print_terminal_results()
        if autoPlots:            
            visualizer.plot_collocation_points()
            visualizer.plot_harmonic_coefficients()
            visualizer.plot_surface_charge_density(comsol_data=self.comsol_data)
            visualizer.plot_convergence_rates(self.mtl_copy, nf_max=12)

        self.mom_collocation_data = {
            freq: {'c': self.c_factor * solver.mom_data['collocation']['maxwellian_capacitance'].item()} for freq in self.freq_range['mom']}
        
        self.mom_galerkin_data = {
            freq: {'c': self.c_factor * solver.mom_data['galerkin']['maxwellian_capacitance'].item()} for freq in self.freq_range['mom']}

    def run_analytical(self):
        """
        Runs the analytical simulation of the transmission line impedance.
        """
        print("\n==============         Analytical Processing       =============")

        f_ana = self.freq_range['ana']
        mtl_obj = MulticonductorTransmissionLine(self.mtl_copy)
        bifilar_wires = WiresHomogeneousMedia(mtl_obj, f_ana)

        le_wires = bifilar_wires.n_wires_external_inductance()['L_ext']
        c_wires  = bifilar_wires.n_wires_capacitance(le_wires)['C_pul']
        pul_bifilar = bifilar_wires.bifilar_static_params()

        impedance = bifilar_wires.bifilar_series_impedance()
        Zs  = impedance['series_impedance_matrix']   # (N_freq, 1, 1)
        Rhf = impedance['high_frequency_limit']      # (N_freq, 1, 1): R_hf + jw*L_ext

        for i, freq in enumerate(f_ana):
            self.analytical_data[freq] = {
                'rhf':       self.r_factor * np.real(Rhf[i, 0, 0]),
                'rs':        self.r_factor * np.real(Zs[i, 0, 0]),
                'ls':        self.l_factor * np.imag(Zs[i, 0, 0]) / (2 * np.pi * freq),
                'le_wires':  self.l_factor * le_wires.item(),
                'le_exact':  self.l_factor * pul_bifilar['inductance']['exact'],
                'le_bifilar':self.l_factor * pul_bifilar['inductance']['approximate'],
                'c_exact':   self.c_factor * pul_bifilar['capacitance']['exact'],
                'c_approx':  self.c_factor * pul_bifilar['capacitance']['approximate'],
                'c_wires':   self.c_factor * c_wires.item(),
            }

    def run_mom_so(self):
        """
        Runs the impedance simulation using the Method of Moments (MoM-SO).

        Args:
            mtl_config (dict): Transmission line configuration dictionary.
            frequencies (np.ndarray): Array of frequencies for the analysis.
            green_mode (GreenFunctionMode): The computation mode for the Green's function.

        Returns:
            list: A list containing the total series impedances computed via MoM.
        """
        print("\n=============   MoM-SO HomogeneousLosslessMedium   =============")

        mtl_model = MulticonductorTransmissionLine(self.mtl_copy)
        numerical_freqs = self.freq_range['mom']

        # The Green's matrix is frequency-independent
        green_matrix = QuasiStatic(mtl_model).green_matrix()
        
        # Instantiate the model ONCE with the full array of numerical frequencies
        mom_so = HomogeneousLosslessMedium(mtl_model, numerical_freqs)
        post_processor = LosslessPostProcessing(mtl_model)

        # Calculate partial impedance for all frequencies
        z_partial = mom_so.z_partial(green_matrix)
        
        # Calculate total series impedance for all frequencies
        zs_stack = post_processor.z_total(z_partial)
        
        # Calculate series resistance and inductance for all frequencies
        rs_stack = post_processor.rs_matrix(zs_stack)
        ls_stack = post_processor.ls_matrix(zs_stack, numerical_freqs)
        general_cap = mom_so.generalized_capacitance_matrix(green_matrix)
        maxwell_cap = mom_so.maxwellian_capacitance_matrix(general_cap)
            
        for i, freq in enumerate(numerical_freqs):
            self.mom_so_data[freq] = {
                'zs': zs_stack[i].item(),
                'rs': rs_stack[i].item() * self.r_factor,
                'ls': ls_stack[i].item() * self.l_factor,
                'c': np.real(maxwell_cap.item()) * self.c_factor,
            }

        if green_matrix.shape[0] < 6:
            from scipy.constants import epsilon_0
            print(f"\nGreen's Matrix (Dim: {green_matrix.shape}):\n{green_matrix}")
            print(f"\n2*pi*e0*G:\n{- 2 * np.pi * epsilon_0 * np.real(green_matrix)}")

        print(f"\nGreen's Matrix Shape: {green_matrix.shape}.")
        matrix_viewer(np.real(general_cap), "MoM-SO Generalized Capacitance Matrix (F/m)")
        matrix_viewer(np.real(maxwell_cap), "MoM-SO Bifilar Capacitance (F/m)")

    def run_srw_rates(self):
        """
        Runs the analytical simulation of the transmission line impedance.

        Args:
            mtl_config (dict): Transmission line configuration dictionary.
            frequencies (np.ndarray): Array of frequencies for the analysis.

        Returns:
            tuple: A tuple containing three lists: series impedances,
                high-frequency resistances and external inductances.
        """
        print("\n==============         SRW RATES EVALUATION        =============")

        mtl_local = copy.deepcopy(self.mtl_copy)
        _dummy_f = np.array([1.0])  # static params don't depend on frequency

        for ratio in self.srw_ratios['ana']:
            separation = ratio * mtl_local[0]['radius'][1]
            mtl_local[1]['center_point'] = (separation, 0.0)

            mtl_obj = MulticonductorTransmissionLine(mtl_local)
            wires = WiresHomogeneousMedia(mtl_obj, _dummy_f)
            le_wires = wires.n_wires_external_inductance()['L_ext']
            c_wires  = wires.n_wires_capacitance(le_wires)['C_pul']
            pul_bifilar = wires.bifilar_static_params()

            self.srw_data[ratio] = {
                'le_wires':  self.l_factor * le_wires.item(),
                'le_exact':  self.l_factor * pul_bifilar['inductance']['exact'],
                'le_bifilar':self.l_factor * pul_bifilar['inductance']['approximate'],
                'c_exact':   self.c_factor * pul_bifilar['capacitance']['exact'],
                'c_approx':  self.c_factor * pul_bifilar['capacitance']['approximate'],
                'c_wires':   self.c_factor * c_wires.item(),
            }
        
        for ratio in self.srw_ratios['mom']:
            separation = ratio * mtl_local[0]['radius'][1]
            mtl_local[1]['center_point'] = (separation, 0.0)

            mtl_model = MulticonductorTransmissionLine(mtl_local)
            solver = BareWireMoMSolver(mtl_model)
            solver.run_collocation_method()

            self._prepare_fortran_runner(mtl_local)
            self.runner.run_fortran(self.fortran_base_params)

            green_matrix = QuasiStatic(mtl_model).green_matrix()
            mom_so = HomogeneousLosslessMedium(mtl_model, frequencies=np.array([0]))
            gen_cap = mom_so.generalized_capacitance_matrix(green_matrix)
            capacitance = mom_so.maxwellian_capacitance_matrix(gen_cap)

            self.srw_mum_data[ratio] = {
                'c_mom-so': self.c_factor * np.real(capacitance.item()),
                'c_ribbon': self.c_factor * self.runner.CAP0_matrix.item(),
                'c_mom': self.c_factor * solver.mom_data['collocation']['maxwellian_capacitance'].item(),
            }

    def run_convergence(self):
        """Runs the convergence loop for both simulations and stores the results."""
        print(f"\nRunning Convergence Rate until k = {self.sum_max}!")
        
        results = []
        for k in range(0, self.sum_max):
            temp_mtl = copy.deepcopy(self.mtl_copy)
            for key in temp_mtl.keys():
                if isinstance(key, int):
                    temp_mtl[key]['fourier_order'] = k  
                    if temp_mtl['type'] == 'coated_wires':
                        temp_mtl[key]['sheath']['fourier_order'] = k

            # === Fortran Instance ===
            self._prepare_fortran_runner(temp_mtl)
            self.runner.run_fortran(self.fortran_base_params)

            # === MoM Instance ===
            temp_mtl['type'] = 'bare_wires'
            for key in temp_mtl.keys():
                if isinstance(key, int):
                    temp_mtl[key]['fourier_order'] = k
                    temp_mtl[key]['sheath'] = None

            mom_bare_wires_model = MulticonductorTransmissionLine(temp_mtl)
            solver = BareWireMoMSolver(mom_bare_wires_model)
            solver.run_collocation_method()

            # === MoM-SO Instance ===
            green_matrix = QuasiStatic(mom_bare_wires_model).green_matrix()
            mom_so = HomogeneousLosslessMedium(mom_bare_wires_model, frequencies=np.array([0]))
            general_cap = mom_so.generalized_capacitance_matrix(green_matrix)
            maxwell_cap = mom_so.maxwellian_capacitance_matrix(general_cap)
            
            # Collect results
            results.append({
                'k': k,
                'L (RIBBON.FOR)':       self.runner.IND_matrix,
                'C (RIBBON.FOR)':       self.runner.CAP_matrix,
                'C0 (RIBBON.FOR)':      self.runner.CAP0_matrix,
                'CGEN (RIBBON.FOR)':    self.runner.CGEN0_matrix,
                'C0 (BARE-WIRE.PY)':    solver.mom_data['collocation']['maxwellian_capacitance'],
                'CGEN (BARE-WIRE.PY)':  solver.mom_data['collocation']['generalized_capacitance'],
                'C0 (MOM-SO.PY)':       np.real(maxwell_cap),
                'CGEN (MOM-SO.PY)':     np.real(general_cap),
            })
            print(f"  Complete for k = {k}.")

        self.results_df = pd.DataFrame(results).set_index('k')

    def plot_resistance_results(self):
        """
        Generates and displays the simulation result plots in a flexible way,
        organized in subplots.

        Args:
            mtl (dict): Transmission line configuration dictionary.
            freqs (dict): Dictionary containing the frequency arrays for each simulation.
            analytical (dict): Dictionary with the analytical simulation results.
            mom_so (dict): Dictionary with the MoM-SO simulation results.
        """
        resistance_data = {
            'mom-so':       {'data': (self.freq_range.get('mom'), [data['rs']  for data in self.mom_so_data.values()]),     'label': 'MoM-SO'},
            'analytical':   {'data': (self.freq_range.get('ana'), [data['rs']  for data in self.analytical_data.values()]), 'label': '$R_i$'},
            'approx':       {'data': (self.freq_range.get('ana'), [data['rhf'] for data in self.analytical_data.values()]), 'label': '$R_{HF}$'},
        }

        fig, ax = plt.subplots(figsize=self.figsize)
        self._configure_plot_appearance(ax, r'Series Resistance p.u.l. ($\Omega$/km)', resistance_data)
        save_figure(fig, self.results_dir, base_filename='pul_series_resistance')
        plt.tight_layout()

    def plot_inductance_results(self):
        """
        Generates and displays the simulation result plots in a flexible way,
        organized in subplots.

        Args:
            mtl (dict): Transmission line configuration dictionary.
            freqs (dict): Dictionary containing the frequency arrays for each simulation.
            analytical (dict): Dictionary with the analytical simulation results.
            mom_so (dict): Dictionary with the MoM-SO simulation results.
        """
        inductance_data = {
            'ribbon':       {'data': (self.freq_range.get('mom'), [data['le']           for data in self.ribbon_data.values()]),     'label': r'$\ell_{e,RIBBON.FOR}$'},
            'mom-so':       {'data': (self.freq_range.get('mom'), [data['ls']           for data in self.mom_so_data.values()]),     'label': 'MoM-SO'},
            'wires':        {'data': (self.freq_range.get('ana'), [data['le_wires']     for data in self.analytical_data.values()]), 'label': r'$\ell_{e,n+1 \; wires}$'},
            'exactly':      {'data': (self.freq_range.get('ana'), [data['le_exact']     for data in self.analytical_data.values()]), 'label': r'$\ell_{e,Bifilar (Exactly)}$'},
            'bifilar':      {'data': (self.freq_range.get('ana'), [data['le_bifilar']   for data in self.analytical_data.values()]), 'label': r'$\ell_{e,Bifilar (Approx.)}$'},
            'analytical':   {'data': (self.freq_range.get('ana'), [data['ls']           for data in self.analytical_data.values()]), 'label': '$\ell_s$'},
        }

        fig, ax = plt.subplots(figsize=self.figsize)
        self._configure_plot_appearance(ax, 'Series Inductance p.u.l. (mH/km)', inductance_data, yscale='linear')
        save_figure(fig, self.results_dir, base_filename='pul_series_inductance')
        plt.tight_layout()

    def plot_capacitance_results(self):
        """
        Generates and displays the simulation result plots in a flexible way,
        organized in subplots.

        Args:
            mtl (dict): Transmission line configuration dictionary.
            freqs (dict): Dictionary containing the frequency arrays for each simulation.
            analytical (dict): Dictionary with the analytical simulation results.
            mom_so (dict): Dictionary with the MoM-SO simulation results.
        """
        
        capacitante_data = {
            'mom':      {'data': (self.freq_range.get('mom'), [data['c']        for data in self.mom_collocation_data.values()]),        'label': 'MoM'},
            'ribbon':   {'data': (self.freq_range.get('mom'), [data['c']        for data in self.ribbon_data.values()]),     'label': 'RIBBON.FOR'},
            'mom-so':   {'data': (self.freq_range.get('mom'), [data['c']        for data in self.mom_so_data.values()]),     'label': 'MoM-SO'},
            'wires':    {'data': (self.freq_range.get('ana'), [data['c_wires']  for data in self.analytical_data.values()]), 'label': 'n+1 wires'},
            'exactly':  {'data': (self.freq_range.get('ana'), [data['c_exact']  for data in self.analytical_data.values()]), 'label': 'Bifilar (Exactly)'},
            'bifilar':  {'data': (self.freq_range.get('ana'), [data['c_approx'] for data in self.analytical_data.values()]), 'label': 'Bifilar (Approx.)'},
        }

        fig, ax = plt.subplots(figsize=self.figsize)
        self._configure_plot_appearance(ax, 'Capacitance p.u.l. (nF/km)', capacitante_data, yscale='linear')
        save_figure(fig, self.results_dir, base_filename='pul_shunt_capacitance')
        plt.tight_layout()

    def plot_srw_rates(self):
        """
        Generates and displays the simulation result plots in a flexible way,
        organized in subplots.

        Args:
            mtl (dict): Transmission line configuration dictionary.
            freqs (dict): Dictionary containing the frequency arrays for each simulation.
            analytical (dict): Dictionary with the analytical simulation results.
            mom_so (dict): Dictionary with the MoM-SO simulation results.
        """

        capacitante_data = {
            'mom':      {'data': (self.srw_ratios.get('mom'), [data['c_mom']    for data in self.srw_mum_data.values()]),   'label': 'MoM'},
            'ribbon':   {'data': (self.srw_ratios.get('mom'), [data['c_ribbon'] for data in self.srw_mum_data.values()]),   'label': 'RIBBON.FOR'},
            'mom-so':   {'data': (self.srw_ratios.get('mom'), [data['c_mom-so'] for data in self.srw_mum_data.values()]),   'label': 'MoM-SO'},
            'exactly':  {'data': (self.srw_ratios.get('ana'), [data['c_exact']  for data in self.srw_data.values()]),       'label': 'Exactly'},
            'approx':   {'data': (self.srw_ratios.get('ana'), [data['c_approx'] for data in self.srw_data.values()]),       'label': 'Approx.'},
        }

        fig, ax = plt.subplots(figsize=self.figsize)
        for key, data in capacitante_data.items():
            freq, value = data['data']
            label = data['label']
            if freq is not None and value is not None:
                if key == 'ribbon':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='o', markersize=8, fillstyle='none', markeredgecolor='k', zorder=2)
                elif key == 'mom':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='s', markersize=10, fillstyle='none', markeredgecolor='k', zorder=2)
                elif key == 'mom-so':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='o', markersize=3, zorder=2)
                elif key == 'approx':
                    ax.plot(freq, value, label=label, color='k', linestyle='--', linewidth=1.0, zorder=1)
                elif key == 'exactly':
                    ax.plot(freq, value, label=label, color='k', linestyle=':', linewidth=1.0, zorder=1)

        ax.set_xlabel('Ratio of separation to wire radius, s/r$_w$')
        ax.set_ylabel('Per-unit-length capacitance (pF/m)')
        ax.set_xlim(2, 8)
        ax.set_ylim(10, 90)
        ax.legend()
        ax.grid(True, linestyle='--', linewidth=0.5)
        save_figure(fig, self.results_dir, base_filename='ratio_srw_capacitance')
        plt.tight_layout()

    def plot_generalized_capacitance_convergence(self):
        """
        Generates the convergence plot of the generalized capacitance, with separate
        subplots for the CGEN_00 and CGEN_01 terms.
        """
        if self.results_df is None:
            print("Run the simulations first with 'run_convergence()'.")
            return

        assert self.results_df.index.max() == self.sum_max - 1, \
            f"The simulation did not run up to the expected maximum value of k={self.sum_max - 1}"

        ribbon, mom, mom_so = {}, {}, {}
        try:
            for i in range(self.N):
                for j in range(self.N):
                    ribbon[f'c_{i}{j}'] =   self._extract_matrix_element('CGEN (RIBBON.FOR)', row=i, col=j)
                    mom_so[f'c_{i}{j}'] =   self._extract_matrix_element('CGEN (MOM-SO.PY)', row=i, col=j)
                    mom[f'c_{i}{j}'] =      self._extract_matrix_element('CGEN (BARE-WIRE.PY)', row=i, col=j)
        except (TypeError, IndexError, ValueError, KeyError) as e:
            print(f"Error while extracting elements of the capacitance matrix: {e}")
            print("Check that the simulations have been run and that the 'CGEN' matrices were populated with the correct dimensions.")
            return

        plt.style.use('default')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=True)
        # fig.suptitle('')

        fortran_nf_axis = self.results_df.index + 1
        python_nf_axis = 2 * self.results_df.index + 1
        max_nf_fortran = fortran_nf_axis.max()
        mask_py = python_nf_axis <= max_nf_fortran
        
        # --- Subplot 1: CGEN_00 ---
        ax1.plot(fortran_nf_axis, ribbon['c_00'], label='RIBBON.FOR',
                color=self.plot_params['colors'][0], marker=self.plot_params['markers'][0], linestyle=self.plot_params['linestyles'][0])
        ax1.plot(python_nf_axis[mask_py], mom['c_00'][mask_py], label='MoM.PY',
                color=self.plot_params['colors'][1], marker=self.plot_params['markers'][1], linestyle=self.plot_params['linestyles'][1], fillstyle='none')
        ax1.plot(python_nf_axis[mask_py], mom_so['c_00'][mask_py], label='MoM-SO.PY',
                color=self.plot_params['colors'][2], marker=self.plot_params['markers'][2], linestyle=self.plot_params['linestyles'][2], fillstyle='none')

        # --- Subplot 2: CGEN_12 ---
        ax2.plot(fortran_nf_axis, ribbon['c_01'], label='RIBBON.FOR',
                color=self.plot_params['colors'][0], marker=self.plot_params['markers'][0], linestyle=self.plot_params['linestyles'][0])
        ax2.plot(python_nf_axis[mask_py], mom['c_01'][mask_py], label='MoM.PY',
                color=self.plot_params['colors'][1], marker=self.plot_params['markers'][1], linestyle=self.plot_params['linestyles'][1], fillstyle='none')
        ax2.plot(python_nf_axis[mask_py], mom_so['c_01'][mask_py], label='MoM-SO.PY',
                color=self.plot_params['colors'][2], marker=self.plot_params['markers'][2], linestyle=self.plot_params['linestyles'][2], fillstyle='none')

        # Axis configuration for both subplots
        for ax in [ax1, ax2]:
            ax.set_xlabel('Number of Fourier Coefficients (NF)', fontsize=11)
            ax.set_xlim(0.8, max_nf_fortran + 0.2)
            ax.set_xticks(np.arange(1, max_nf_fortran + 1, 1))
            ax.tick_params(top=True, right=True, direction='in', which='both')
            ax.legend(loc='lower right', fontsize=10)
            ax.grid(False)

        # Subplot-specific settings
        ax1.set_ylabel('Generalized Capacitance Matrix, $CGEN$ (pF/m)', fontsize=11)
        ax1.set_title('Auto-Capacitance Term $CGEN_{00}$')
        ax2.set_title('Mutual Capacitance Term $CGEN_{01}$')
        save_figure(fig, self.results_dir, base_filename='generalized_capacitance_convergence')
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def plot_free_space_capacitance_convergence(self):
        """ Generates the free-space capacitance convergence plot, adapting to the number of conductors in the system. """
        if self.results_df is None:
            print("Run the simulations first with 'run_study()'.")
            return

        assert self.results_df.index.max() == self.sum_max - 1, \
            f"The simulation did not run up to the expected maximum value of k={self.sum_max - 1}"

        # Programmatic data extraction
        ribbon, mom, mom_so = {}, {}, {}
        try:
            for i in range(self.N-1):
                for j in range(self.N-1):
                    ribbon[f'c_{i}{j}'] =   self._extract_matrix_element('C0 (RIBBON.FOR)', row=i, col=j)
                    mom_so[f'c_{i}{j}'] =   self._extract_matrix_element('C0 (MOM-SO.PY)', row=i, col=j)
                    mom[f'c_{i}{j}'] =      self._extract_matrix_element('C0 (BARE-WIRE.PY)', row=i, col=j)

        except (TypeError, IndexError, ValueError) as e:
            print(f"Error while extracting elements of the capacitance matrix: {e}")
            print("Check that the simulations have been run and that the 'C' matrix was populated.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)
        fortran_nf_axis = self.results_df.index + 1
        python_nf_axis = 2 * self.results_df.index + 1
        max_nf_fortran = fortran_nf_axis.max()
        mask_py = python_nf_axis <= max_nf_fortran
        
        legend_items_count = 0
        for i in range(self.N-1):
            legend_items_count += 3
            idx = i % len(self.plot_params['markers'])
            
            # Plot of the main diagonal
            ax.plot(fortran_nf_axis, ribbon[f'c_{i}{i}'], label='RIBBON.FOR',
                    color=self.plot_params['colors'][0], marker=self.plot_params['markers'][idx], linestyle=self.plot_params['linestyles'][0])

            ax.plot(python_nf_axis[mask_py], mom[f'c_{i}{i}'][mask_py], label='MoM.PY',
                    color=self.plot_params['colors'][1], marker=self.plot_params['markers'][idx], linestyle=self.plot_params['linestyles'][1], fillstyle='none')

            ax.plot(python_nf_axis[mask_py], mom_so[f'c_{i}{i}'][mask_py], label='MoM-SO.PY',
                    color=self.plot_params['colors'][2], marker=self.plot_params['markers'][idx], linestyle=self.plot_params['linestyles'][2], fillstyle='none')

            # Plot of the off-diagonal elements
            for j in range(i + 1, self.N-1):
                legend_items_count += 3
                ijdx = (i+j) % len(self.plot_params['markers'])
                ax.plot(fortran_nf_axis, ribbon[f'c_{i}{j}'], label=f'$C_{{{i+1}{j+1}}} (RIBBON.FOR)$',
                        color=self.plot_params['colors'][3], marker=self.plot_params['markers'][ijdx], linestyle=self.plot_params['linestyles'][0])

                ax.plot(python_nf_axis[mask_py], mom[f'c_{i}{j}'][mask_py], label=f'$C_{{{i+1}{j+1}}} (MoM.PY)$',
                        color=self.plot_params['colors'][4], marker=self.plot_params['markers'][ijdx], linestyle=self.plot_params['linestyles'][1], fillstyle='none')

                ax.plot(python_nf_axis[mask_py], mom_so[f'c_{i}{j}'][mask_py], label=f'$C_{{{i+1}{j+1}}} (MoM-SO.PY)$',
                        color=self.plot_params['colors'][5], marker=self.plot_params['markers'][ijdx], linestyle=self.plot_params['linestyles'][2], fillstyle='none')

        ax.set_title('Auto-Capacitance Term $C0_{11}$')
        ax.set_xlabel('Number of Fourier Coefficients (NF)', fontsize=12)
        ax.set_ylabel('Free-Space Capacitance Matrix, $C0$ (pF/m)', fontsize=12)
        ax.set_xlim(0.8, max_nf_fortran + 0.2)
        ax.set_xticks(np.arange(1, max_nf_fortran + 1, 1))
        ax.tick_params(top=True, right=True, direction='in', which='both')
        # ax.legend(loc='upper center', fontsize=9, ncol=legend_items_count, bbox_to_anchor=(0.5, 1.1), fancybox=True)
        ax.legend(loc='lower right', fontsize=10)
        ax.grid(False)
        save_figure(fig, self.results_dir, base_filename='free_space_capacitance_convergence')
        plt.tight_layout()
    