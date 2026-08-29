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
from tulip.py_tulip import PyTulip
from analytical_forms.isolated_wires import WiresHomogeneousMedia
from mom.coated_wire_systems import MulticonductorCoatedWireSystems

class RibbonCoatedCablesPULParameters:
    """
    Encapsulates the logic to run and analyze the capacitance convergence
    study, comparing Python MoM and Fortran.
    """
    def __init__(self, project_root: Path, case_name: str, mtl: Dict[str, Any], SUM_MAX: int = 10):
        """
        Initializes the convergence analyzer.

        Args:
            mtl_config (Dict[str, Any]): Dictionary with the MTL model configuration.
            nf_max (int): Maximum number of coefficients / harmonic order to test.
        """

        assert len([key for key in mtl.keys() if isinstance(key, int)]) > 1, "The line must contain more than one conductor."

        self.project_root = project_root
        self.case_name = case_name
        self.script_file_path = Path(__file__).parent / f'{case_name}.py'
        self.mtl_copy = copy.deepcopy(mtl)
        self.freq_range = {'ana': np.logspace(0, 6, num=200), 'mom': np.logspace(0, 6, num=30)}
        self.srw_ratios = {'ana': np.linspace(4.0, 10.0, num=300), 'mom': np.linspace(4.0, 10.0, num=40)}

        # Data parameters
        self.srw_data = {}
        self.srw_mum_data = {}
        self.analytical_data = {}
        self.mom_data = {}
        self.mom_so_data = {}
        self.ribbon_data = {}
        self.tulip_runner = None

        # Additional parameters
        self.c_factor = 1e12  # F/m to nF/km
        self.l_factor = 1e6   # H/m to mH/km
        self.r_factor = 1e3   # Ohm/m to Ohm/km
        self.figsize = (12, 5)
        self.pt1 = 63
        self.pt2 = 10

        self.N = len([key for key in mtl.keys() if isinstance(key, int)])
        self.sum_max = SUM_MAX
        self.results_df = None
        
        # Plotting parameters
        self.plot_params = {
            'linestyles': [':', '-.', '--', '-', ':', '-.', '--'],
            'markers': ['o', 'd', 's', '^', 'v', '<', '>'],
            'colors': ['black', 'gray', 'lightgray', 'gainsboro']
        }

        self.R = self.mtl_copy[0]['radius'][1]
        self.D = np.linalg.norm(np.array(self.mtl_copy[0]['center_point']) - np.array(self.mtl_copy[1]['center_point']))
        self.swr_ratio = self.D / self.R

    def _prepare_fortran_runner(self, mtl):
        """Prepares the parameters and the runner for the Fortran simulation."""

        # Get the 'insulation' dictionary safely.
        #    If 'insulation' does not exist or is None, use an empty dictionary {} as a fallback.
        refIdx = mtl['idx_ref_conductor']
        insulation_dict = mtl[refIdx].get('insulation') or {}

        self.fortran_base_params = {
            'N':    len([key for key in mtl.keys() if isinstance(key, int)]),
            'NF':   mtl[refIdx]['fourier_order'] + 1,
            'IREF': refIdx + 1,
            'RW':   mtl[refIdx]['radius'][1],
            'TD':   insulation_dict.get('thickness', 0.0),
            'ER':   insulation_dict.get('relative_permittivity', 1.0),
            'S':    np.linalg.norm(np.array(mtl[0]['center_point']) - np.array(mtl[1]['center_point'])),
        }

        fortran_exe_path = self.project_root / 'mtl_paul' / 'RIBBON' / 'RIBBON.EXE'
        self.runner = FortranRunner(exe_path=str(fortran_exe_path), silent=True)

    def _prepare_tulip_runner(self, case_name: str = "three_wires_ribbon"):
        """Prepares the runner for the Tulip simulation."""
        EXE_PATH = r"C:\git\tulip\pulmtln-build\rls\bin\Release\pulmtln.exe"
        print(f"--- Preparing to run simulation for case: '{case_name}' ---")
        
        try:
            self.tulip_runner = PyTulip(
                executable_path=EXE_PATH,
                case_name=case_name,
                silent=False,
                overwrite_mesh=False
            )

        except (FileNotFoundError, NotADirectoryError) as e:
            print(f"\nConfiguration Error: {e}")
        except Exception as e:
            print(f"\nAn unexpected error occurred: {e}")
        
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
        print("="*self.pt2 + " BIFILAR COATED-WIRE RIBBON CABLE SIMULATION " + "="*self.pt2)
        print(f"Project: {self.project_root}")
        print(f"Model: {self.mtl_copy['name']}")
        print(f"D/R = {self.swr_ratio:.3f}. Fourier Order (k) = {self.mtl_copy[0]['fourier_order']}.")
        print(f"RIBBON Fourier Coef./cond. (NF) = {self.mtl_copy[0]['fourier_order']+1}.")
        print(f"PYTHON Fourier Coef./cond. (NF) = {2*self.mtl_copy[0]['fourier_order']+1}.")
        print("="*self.pt1)

    def run_single_fortran(self, mtl: dict = None, displayTerminal: bool = True):
        """Runs a single simulation for a specific value of k."""

        if mtl is None:
            self._prepare_fortran_runner(self.mtl_copy)
        else:
            self._prepare_fortran_runner(mtl)

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

    def run_fortran(self, mtl: dict = None, displayTerminal: bool = True):
        """Runs a single simulation for a specific value of k."""

        if mtl is None:
            self._prepare_fortran_runner(self.mtl_copy)
        else:
            self._prepare_fortran_runner(mtl)

        self._prepare_fortran_runner(self.mtl_copy)
        self.runner.run_fortran(self.fortran_base_params)

        self.ribbon_data = {
            freq: {'c': self.c_factor * self.runner.CAP0_matrix[0, 0], 
                   'le': self.l_factor * self.runner.IND_matrix[0, 0]} for freq in self.freq_range['mom']}

    def run_single_tulip(self):
        """
        Runs the simulation using PyTulip and processes the results.
        """
        print("\n============  pyTulip Simulation  ===========")
        if self.tulip_runner is None:
            self._prepare_tulip_runner()

        if self.tulip_runner is None:
            print("Warning: PyTulip not available (executable not found). Step skipped.")
            return

        self.tulip_runner.run()
        
        if not self.tulip_runner.output_filepath.is_file():
            print(f"Error: Output file '{self.tulip_runner.output_filepath.name}' was not generated.")
            return
        
        print("\n--- Accessing Results ---")
        if self.tulip_runner.C_matrix is not None:
            matrix_viewer(self.tulip_runner.C_matrix, "Capacitance Matrix (C)")

        if self.tulip_runner.L_matrix is not None:
            matrix_viewer(self.tulip_runner.L_matrix, "Inductance Matrix (L)")

    def run_mom_methods(self):
        """
        Runs the classical Method of Moments (MoM) simulation for the bifilar transmission line.

        Returns:
            BifilarMoM: Instance of the configured BifilarMoM object.
        """
        print("\n============  pyMoM TwoCoatedWireSystem (Bare-Wire)  ===========")
        bare_wire_mtl = copy.deepcopy(self.mtl_copy)
        bare_wire_mtl['type'] = 'bare_wires'
        for key in bare_wire_mtl.keys():
            if isinstance(key, int):
                bare_wire_mtl[key]['insulation'] = None

        bare_wires_model = MulticonductorTransmissionLine(bare_wire_mtl)
        bare_wires = MulticonductorCoatedWireSystems(bare_wires_model)
        bare_wires.run_simulation()
        bare_wires.print_results()

        print("\n==============      pyMoM TwoCoatedWireSystem      =============")
        coated_wires_model = MulticonductorTransmissionLine(self.mtl_copy)
        IsolatedMTLRepresentation(self.script_file_path, coated_wires_model, units='millimeter').system_schematic()
        coated_wires = MulticonductorCoatedWireSystems(coated_wires_model)
        coated_wires.run_simulation()
        coated_wires.print_results()
        coated_wires.plot_collocation_points()

        self.mom_data = {
            freq: {'c_bare_wire': self.c_factor * bare_wires.C_maxwellian,
                    'c_coated_wire': self.c_factor * coated_wires.C_maxwellian} for freq in self.freq_range['mom']}

    def srw_rates(self):
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

        # Build the bifilar line configuration dynamically for each ratio.
        for ratio in self.srw_ratios['ana']:
            temp_mtl = copy.deepcopy(self.mtl_copy)
            separation = ratio * temp_mtl[0]['radius'][1]
            temp_mtl[1]['center_point'] = (separation, 0.0)

            wires = WiresHomogeneousMedia(temp_mtl)
            le_wires = wires.n_wires_external_inductance_matrix()
            pul_bifilar = wires.bifilar_pul_inductance_and_capacitance()

            self.srw_data[ratio] = {
                'le_wires':     self.l_factor * le_wires,
                'le_exact':     self.l_factor * pul_bifilar['inductance']['exact'],
                'le_bifilar':   self.l_factor * pul_bifilar['inductance']['approximate'],
                'c_exact':      self.c_factor * pul_bifilar['capacitance']['exact'],
                'c_approx':     self.c_factor * pul_bifilar['capacitance']['approximate'],
                'c_wires':      self.c_factor * wires.n_wires_capacitance_matrix(le_wires),
            }

        for ratio in self.srw_ratios['mom']:
            temp_mtl = copy.deepcopy(self.mtl_copy)
            separation = ratio * temp_mtl[0]['radius'][1]
            temp_mtl[1]['center_point'] = (separation, 0.0)

            # === Fortran RIBBON Instance ===            
            self._prepare_fortran_runner(temp_mtl)
            self.runner.run_fortran(self.fortran_base_params)

            # === MoM TwoCoatedWireSystem Instance ===
            mom_coated_model = MulticonductorTransmissionLine(temp_mtl)
            mom_coated = MulticonductorCoatedWireSystems(mom_coated_model)
            mom_coated.run_simulation()
            
            temp_mtl['type'] = 'bare_wires'
            for key in temp_mtl.keys():
                if isinstance(key, int):
                    temp_mtl[key]['insulation'] = None

            mom_bare_model = MulticonductorTransmissionLine(temp_mtl)
            mom_bare = MulticonductorCoatedWireSystems(mom_bare_model)
            mom_bare.run_simulation()

            self.srw_mum_data[ratio] = {
                'ribbon-c':     self.c_factor * self.runner.CAP_matrix.item(),
                'ribbon-c0':    self.c_factor * self.runner.CAP0_matrix.item(),
                'mom-c':        self.c_factor * mom_coated.C_maxwellian.item(),
                'mom-c0':       self.c_factor * mom_bare.C_maxwellian.item(),
            }

    def run_convergence(self):
        """Runs the convergence loop for both simulations and stores the results."""
        print(f"\nRunning Convergence Rate until k = {self.sum_max}!")

        if self.tulip_runner is None:
            self._prepare_tulip_runner()

        if self.tulip_runner is None:
            print("Warning: PyTulip not available. Tulip results will be omitted.")
        else:
            self.tulip_runner.run()

        results = []
        for k in range(0, self.sum_max):
            temp_mtl = copy.deepcopy(self.mtl_copy)
            for key in temp_mtl.keys():
                if isinstance(key, int):
                    temp_mtl[key]['fourier_order'] = k  
                    temp_mtl[key]['insulation']['fourier_order'] = k

            # === Fortran RIBBON Instance ===
            self._prepare_fortran_runner(temp_mtl)
            self.runner.run_fortran(self.fortran_base_params)

            # === MoM TwoCoatedWireSystem Instance ===
            mom_coated_model = MulticonductorTransmissionLine(temp_mtl)
            mom_coated = MulticonductorCoatedWireSystems(mom_coated_model)
            mom_coated.run_simulation()

            results.append({
                'k': k,
                'L (RIBBON.FOR)':    self.runner.IND_matrix if self.fortran_base_params is not None else np.nan,
                'C (RIBBON.FOR)':    self.runner.CAP_matrix if self.fortran_base_params is not None else np.nan,
                'C0 (RIBBON.FOR)':   self.runner.CAP0_matrix if self.fortran_base_params is not None else np.nan,
                'CGEN (RIBBON.FOR)': self.runner.CGEN0_matrix if self.fortran_base_params is not None else np.nan,
                'C (MoM.PY)':        mom_coated.C_maxwellian if mom_coated.C_maxwellian is not None else np.nan,
                'C (SEMBA-TULIP)':   self.tulip_runner.C_matrix if self.tulip_runner is not None and self.tulip_runner.C_matrix is not None else np.nan,
                'L (SEMBA-TULIP)':   self.tulip_runner.L_matrix if self.tulip_runner is not None and self.tulip_runner.L_matrix is not None else np.nan,
            })
            print(f"  Complete for k = {k}.")

        self.results_df = pd.DataFrame(results).set_index('k')

    def plot_inductance_data(self):
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
            'wires':        {'data': (self.freq_range.get('ana'), [data['le_wires']     for data in self.srw_data.values()]), 'label': r'$\ell_{e,n+1 \; wires}$'},
            'exactly':      {'data': (self.freq_range.get('ana'), [data['le_exact']     for data in self.srw_data.values()]), 'label': r'$\ell_{e,Bifilar (Exactly)}$'},
            'bifilar':      {'data': (self.freq_range.get('ana'), [data['le_bifilar']   for data in self.srw_data.values()]), 'label': r'$\ell_{e,Bifilar (Approx.)}$'},
            'analytical':   {'data': (self.freq_range.get('ana'), [data['ls']           for data in self.srw_data.values()]), 'label': '$\ell_s$'},
        }

        fig, ax = plt.subplots(figsize=self.figsize)
        self._configure_plot_appearance(ax, 'Series Inductance p.u.l. (mH/km)', inductance_data, yscale='linear')
        plt.tight_layout()

    def plot_capacitance_data(self):
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
            'mom':      {'data': (self.freq_range.get('mom'), [data['c']        for data in self.mom_data.values()]),        'label': 'MoM'},
            'ribbon':   {'data': (self.freq_range.get('mom'), [data['c']        for data in self.ribbon_data.values()]),     'label': 'RIBBON.FOR'},
            'mom-so':   {'data': (self.freq_range.get('mom'), [data['c']        for data in self.mom_so_data.values()]),     'label': 'MoM-SO'},
            'wires':    {'data': (self.freq_range.get('ana'), [data['c_wires']  for data in self.srw_data.values()]), 'label': 'n+1 wires'},
            'exactly':  {'data': (self.freq_range.get('ana'), [data['c_exact']  for data in self.srw_data.values()]), 'label': 'Bifilar (Exactly)'},
            'bifilar':  {'data': (self.freq_range.get('ana'), [data['c_approx'] for data in self.srw_data.values()]), 'label': 'Bifilar (Approx.)'},
        }

        fig, ax = plt.subplots(figsize=self.figsize)
        self._configure_plot_appearance(ax, 'Capacitance p.u.l. (nF/km)', capacitante_data, yscale='linear')
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
            'ribbon-c':     {'data': (self.srw_ratios.get('mom'), [data['ribbon-c']     for data in self.srw_mum_data.values()]), 'label': '$C$ (RIBBON.FOR)'},
            'mom-c':        {'data': (self.srw_ratios.get('mom'), [data['mom-c']        for data in self.srw_mum_data.values()]), 'label': '$C$ (MoM.PY)'},
            'ribbon-c0':    {'data': (self.srw_ratios.get('mom'), [data['ribbon-c0']    for data in self.srw_mum_data.values()]), 'label': '$C_0$ (RIBBON.FOR)'},
            'mom-c0':       {'data': (self.srw_ratios.get('mom'), [data['mom-c0']       for data in self.srw_mum_data.values()]), 'label': '$C_0$ (MoM.PY)'},
            'exactly':      {'data': (self.srw_ratios.get('ana'), [data['c_exact']      for data in self.srw_data.values()]),     'label': 'Exactly'},
            'approx':       {'data': (self.srw_ratios.get('ana'), [data['c_approx']     for data in self.srw_data.values()]),     'label': 'Approx.'},
        }

        fig, ax = plt.subplots(figsize=self.figsize)
        for key, data in capacitante_data.items():
            freq, value = data['data']
            label = data['label']
            if freq is not None and value is not None:
                if key == 'ribbon-c':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='s', markersize=5)
                elif key == 'ribbon-c0':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='o', markersize=4)
                elif key == 'mom-c':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='s', markersize=9, fillstyle='none', markeredgecolor='k', zorder=2)
                elif key == 'mom-c0':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='o', markersize=9, fillstyle='none', markeredgecolor='k', zorder=2)
                elif key == 'approx':
                    ax.plot(freq, value, label=label, color='k', linestyle='--', linewidth=1.0, zorder=1)
                elif key == 'exactly':
                    ax.plot(freq, value, label=label, color='k', linestyle=':', linewidth=1.0, zorder=1)

        ax.set_xlabel('Ratio of separation to wire radius, s/r$_w$')
        ax.set_ylabel('Per-unit-length capacitance (pF/m)')
        ax.set_xlim(4, 8)
        ax.legend()
        ax.grid(True, linestyle='--', linewidth=0.5)
        plt.tight_layout()
        
    def plot_paul_fig514a(self):
        """
        Generates the inductance convergence plot from the stored results,
        replicating the reference figure.
        """
        if self.results_df is None:
            print("Run the simulations first with 'run_convergence()'.")
            return

        # Ensure the simulation ran to completion
        assert self.results_df.index.max() == self.sum_max - 1, \
            f"The simulation did not run up to the expected maximum value of k={self.sum_max - 1}"

        # Extract the components of the 'L' inductance matrix and convert to uH/m.
        try:
            # Pass the correct conversion factor (self.l_factor)
            l11_ribbon = self._extract_matrix_element('L (RIBBON.FOR)', row=0, col=0)
            l22_ribbon = self._extract_matrix_element('L (RIBBON.FOR)', row=1, col=1)
            l12_ribbon = self._extract_matrix_element('L (RIBBON.FOR)', row=0, col=1)
            l11_tulip  = self._extract_matrix_element('L (SEMBA-TULIP)', row=0, col=0)
            l22_tulip  = self._extract_matrix_element('L (SEMBA-TULIP)', row=1, col=1)
            l12_tulip  = self._extract_matrix_element('L (SEMBA-TULIP)', row=0, col=1)
        except (TypeError, IndexError) as e:
            # Adjust the error messages for the inductance context
            print(f"Error while extracting elements of the inductance matrix: {e}")
            print("Check that the simulations have been run and that the 'L (RIBBON.FOR)' matrix was populated.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)

        for_nf_axis = self.results_df.index + 1
        ax.plot(for_nf_axis, l22_ribbon, color='k', marker='o', linestyle=':', label='$L_{22}$ (RIBBON.FOR)')
        ax.plot(for_nf_axis, l11_ribbon, color='k', marker='o', linestyle=':', label='$L_{11}$ (RIBBON.FOR)', markerfacecolor='white', markeredgecolor='k')
        ax.plot(for_nf_axis, l12_ribbon, color='k', marker='s', linestyle=':', label='$L_{12}$ (RIBBON.FOR)')
        ax.plot(for_nf_axis, l22_tulip, color='gray', linestyle='-.', label='(SEMBA-TULIP)')
        ax.plot(for_nf_axis, l11_tulip, color='gray', linestyle='-.')
        ax.plot(for_nf_axis, l12_tulip, color='gray', linestyle='-.')

        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Inductance (µH/m)', fontsize=12)
        ax.set_xlim(0.8, self.sum_max + 0.2)
        ax.set_xticks(np.arange(1, self.sum_max + 1, 1))        
        ax.set_yticks(np.arange(0.2, 1.2, 0.1))        
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper center', fontsize=9, markerscale=1.0, ncol=4, bbox_to_anchor=(0.5, 1.12), fancybox=True)
        plt.tight_layout()

    def plot_paul_fig514b(self):
        """
        Generates the convergence plot from the stored results,
        replicating the reference figure.
        """
        if self.results_df is None:
            print("Run the simulations first with 'run_convergence()'.")
            return

        # Ensure the simulation ran to completion
        assert self.results_df.index.max() == self.sum_max - 1, \
            f"The simulation did not run up to the expected maximum value of k={self.sum_max - 1}"

        # Extract the components of the 'L' inductance matrix and convert to uH/m.
        try:
            # Pass the correct conversion factor (self.l_factor)
            c11_ribbon = self._extract_matrix_element('C (RIBBON.FOR)', row=0, col=0)
            c22_ribbon = self._extract_matrix_element('C (RIBBON.FOR)', row=1, col=1)
            c12_ribbon = self._extract_matrix_element('C (RIBBON.FOR)', row=0, col=1)
            c11_tulip  = self._extract_matrix_element('C (SEMBA-TULIP)', row=0, col=0)
            c22_tulip  = self._extract_matrix_element('C (SEMBA-TULIP)', row=1, col=1)
            c12_tulip  = self._extract_matrix_element('C (SEMBA-TULIP)', row=0, col=1)
        except (TypeError, IndexError) as e:
            # Adjust the error messages for the inductance context
            print(f"Error while extracting elements of the inductance matrix: {e}")
            print("Check that the simulations have been run and that the 'L (RIBBON.FOR)' matrix was populated.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)

        for_nf_axis = self.results_df.index + 1
        ax.plot(for_nf_axis, c11_ribbon, color='k', marker='o', linestyle=':', label='$C_{11}$ (RIBBON.FOR)')
        ax.plot(for_nf_axis, c22_ribbon, color='k', marker='o', linestyle=':', label='$C_{22}$ (RIBBON.FOR)', markerfacecolor='white', markeredgecolor='k')
        ax.plot(for_nf_axis, c12_ribbon, color='k', marker='s', linestyle=':', label='$C_{12}$ (RIBBON.FOR)')
        ax.plot(for_nf_axis, c11_tulip, color='gray', linestyle='-.',  label='(SEMBA-TULIP)')
        ax.plot(for_nf_axis, c22_tulip, color='gray', linestyle='-.')
        ax.plot(for_nf_axis, c12_tulip, color='gray', linestyle='-.')

        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.sum_max + 1, 1))
        ax.set_yticks(np.arange(15, 40, 5))
        ax.set_xlim(0.8, self.sum_max + 0.2)
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper center', fontsize=9, markerscale=1.0, ncol=4, bbox_to_anchor=(0.5, 1.12), fancybox=True)
        plt.tight_layout()

    def plot_paul_fig514c(self):
        """
        Generates the convergence plot from the stored results,
        replicating the reference figure.
        """
        if self.results_df is None:
            print("Run the simulations first with 'run_convergence()'.")
            return

        # Ensure the simulation ran to completion
        assert self.results_df.index.max() == self.sum_max - 1, \
            f"The simulation did not run up to the expected maximum value of k={self.sum_max - 1}"

        # Extract the components of the 'L' inductance matrix and convert to uH/m.
        try:
            # Pass the correct conversion factor (self.l_factor)
            c0_11 = self._extract_matrix_element('C0 (RIBBON.FOR)', row=0, col=0)
            c0_22 = self._extract_matrix_element('C0 (RIBBON.FOR)', row=1, col=1)
            c0_12 = self._extract_matrix_element('C0 (RIBBON.FOR)', row=0, col=1)
        except (TypeError, IndexError) as e:
            # Adjust the error messages for the inductance context
            print(f"Error while extracting elements of the inductance matrix: {e}")
            print("Check that the simulations have been run and that the 'L (RIBBON.FOR)' matrix was populated.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)

        # Plot the data with the reference figure style
        fortran_nf_axis = self.results_df.index + 1
        ax.plot(fortran_nf_axis, c0_11, color='k', marker='o', linestyle='-',  label='$C0_{11}$')
        ax.plot(fortran_nf_axis, c0_22, color='k', marker='o', linestyle=':',  label='$C0_{22}$', markerfacecolor='white', markeredgecolor='k')
        ax.plot(fortran_nf_axis, c0_12, color='k', marker='s', linestyle='-.', label='$C0_{12}$')

        # Axis configuration to match the reference image
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Bare-Wire Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.sum_max + 1, 1))
        ax.set_xlim(0.8, self.sum_max + 0.2)
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper center', fontsize=9, markerscale=1.0, ncol=3, bbox_to_anchor=(0.5, 1.12), fancybox=True)
        plt.tight_layout()

    def plot_paul_fig514d(self):
        """ Generates the convergence plot from the stored results, replicating the reference figure. """
        if self.results_df is None:
            print("Run the simulations first with 'run_convergence()'.")
            return

        # Ensure the simulation ran to completion
        assert self.results_df.index.max() == self.sum_max - 1, \
            f"The simulation did not run up to the expected maximum value of k={self.sum_max - 1}"

        # Extract the components of the 'L' inductance matrix and convert to uH/m.
        try:
            # Pass the correct conversion factor (self.l_factor)
            cgen_00 = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=0, col=0)
            cgen_11 = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=1, col=1)
            cgen_22 = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=2, col=2)
            cgen_01 = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=0, col=1)
            cgen_02 = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=0, col=2)
            cgen_12 = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=1, col=2)

        except (TypeError, IndexError) as e:
            # Adjust the error messages for the inductance context
            print(f"Error while extracting elements of the inductance matrix: {e}")
            print("Check that the simulations have been run and that the 'L (RIBBON.FOR)' matrix was populated.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)

        # Plot the data with the reference figure style
        fortran_nf_axis = self.results_df.index + 1
        ax.plot(fortran_nf_axis, cgen_00, color='k', marker='o', label='$CGEN_{00}$', linewidth=1.0, zorder=0, markersize=3, linestyle=':')
        ax.plot(fortran_nf_axis, cgen_11, color='k', marker='o', label='$CGEN_{11}$', linewidth=1.0, zorder=1, markersize=3, linestyle='-')
        ax.plot(fortran_nf_axis, cgen_22, color='k', marker='o', label='$CGEN_{22}$', linewidth=1.0, zorder=1, markersize=8, linestyle=':', fillstyle='none')
        ax.plot(fortran_nf_axis, cgen_01, color='k', marker='s', label='$CGEN_{01}$', linewidth=1.0, zorder=0, markersize=4, linestyle='--')
        ax.plot(fortran_nf_axis, cgen_12, color='k', marker='s', label='$CGEN_{12}$', linewidth=1.0, zorder=1, markersize=8, linestyle='--', fillstyle='none')
        ax.plot(fortran_nf_axis, cgen_02, color='k', marker='d', label='$CGEN_{02}$', linewidth=1.0, zorder=0, markersize=6, linestyle='-.')

        # Axis configuration to match the reference image
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Generalized Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.sum_max + 1, 1))
        ax.set_xlim(0.8, self.sum_max + 0.2)
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper center', fontsize=9, markerscale=1.0, ncol=6, bbox_to_anchor=(0.5, 1.12), fancybox=True)
        plt.tight_layout()

    def plot_generalized_capacitance_convergence(self):
        """
        Generates the convergence plot of the generalized capacitance, with subplots
        separated for the CGEN_00 and CGEN_01 terms.
        """
        if self.results_df is None:
            print("Run the simulations first with 'run_convergence()'.")
            return

        assert self.results_df.index.max() == self.sum_max - 1, \
            f"The simulation did not run up to the expected maximum value of k={self.sum_max - 1}"

        # Programmatic data extraction.
        # This part stays flexible to extract every element,
        # even if only some are plotted.
        ribbon, mom, mom_so = {}, {}, {}
        try:
            for i in range(self.N):
                for j in range(self.N):
                    ribbon[f'c_{i}{j}'] =   self._extract_matrix_element('CGEN (RIBBON.FOR)', row=i, col=j)
        except (TypeError, IndexError, ValueError, KeyError) as e:
            print(f"Error while extracting elements of the capacitance matrix: {e}")
            print("Check that the simulations have been run and that the 'CGEN' matrices were populated with the correct dimensions.")
            return

        plt.style.use('default')
        # Create a figure with two subplots (1 row, 2 columns)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=True)
        # fig.suptitle('')

        fortran_nf_axis = self.results_df.index + 1
        python_nf_axis = 2 * self.results_df.index + 1
        max_nf_fortran = fortran_nf_axis.max()
        mask_py = python_nf_axis <= max_nf_fortran
        
        # --- Subplot 1: CGEN_00 ---
        ax1.plot(fortran_nf_axis, ribbon['c_00'], label='RIBBON.FOR',
                color=self.plot_params['colors'][0], marker=self.plot_params['markers'][0], linestyle=self.plot_params['linestyles'][0])

        # --- Subplot 2: CGEN_12 ---
        ax2.plot(fortran_nf_axis, ribbon['c_01'], label='RIBBON.FOR',
                color=self.plot_params['colors'][0], marker=self.plot_params['markers'][0], linestyle=self.plot_params['linestyles'][0])

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
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def plot_dielectric_coated_capacitance_convergence(self):
        """ Generates the free-space capacitance convergence plot, adapting to the number of conductors in the system. """
        if self.results_df is None:
            print("Run the simulations first with 'run_study()'.")
            return

        assert self.results_df.index.max() == self.sum_max - 1, \
            f"The simulation did not run up to the expected maximum value of k={self.sum_max - 1}"

        # Programmatic data extraction
        ribbon, mom =  {}, {}

        try:
            for i in range(self.N-1):
                for j in range(self.N-1):
                    ribbon[f'c_{i}{j}'] = self._extract_matrix_element('C (RIBBON.FOR)', row=i, col=j)
                    mom[f'c_{i}{j}'] = self._extract_matrix_element('C (MoM.PY)', row=i, col=j)

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
        
        for i in range(self.N-1):
            for j in range(i, self.N-1):
                ax.plot(fortran_nf_axis, ribbon[f'c_{i}{j}'], markersize=4, zorder = 1,
                        label=f'$C_{{{i+1}{j+1}}}$ (RIBBON.FOR)', 
                        color=self.plot_params['colors'][i+j],
                        marker=self.plot_params['markers'][i+j],
                        linestyle=self.plot_params['linestyles'][i+j])

                ax.plot(python_nf_axis[mask_py], mom[f'c_{i}{j}'][mask_py], markersize=10, fillstyle='none', zorder = 2,
                        label=f'$C_{{{i+1}{j+1}}}$ (MoM.PY)', 
                        color=self.plot_params['colors'][i+j],
                        marker=self.plot_params['markers'][i+j],
                        linestyle=self.plot_params['linestyles'][i+j])
            
        ax.set_title('Three Ribbon Cable Capacitance Convergence')
        ax.set_xlabel('Number of Fourier Coefficients (NF)', fontsize=12)
        ax.set_ylabel('Capacitance (pF/m)', fontsize=12)
        ax.set_xlim(0.8, max_nf_fortran + 0.2)
        ax.set_xticks(np.arange(1, max_nf_fortran + 1, 1))
        ax.tick_params(top=True, right=True, direction='in', which='both')
        ax.legend(loc='best', fontsize=10)
        ax.grid(False)
        plt.tight_layout()
    