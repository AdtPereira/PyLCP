import numpy as np
import matplotlib.pyplot as plt
from utils.case_utils import save_figure
from .models_base import BasePlotter


class SingleCoreCableModels(BasePlotter):
    """
    Plotter for the Xue-style coaxial single-core cable (SCC) case studies
    (core + sheath, multiple analytical scenarios overlaid): series
    impedance/admittance matrices, their decomposition into
    internal/earth-return/series/shunt/potential terms, and potential
    coefficients, compared against the Magalhães/Xue integral formulation,
    De Conti/Vance approximations, and (optionally) COMSOL FEM results (Js
    method, core/sheath excitation). Used by ``scc_132kV_xue.py`` and
    ``scc_34kV_andreata.py``.

    Config-driven per the project's BasePlotter convention (see
    plotter/models_base.py): ``plot_config`` is injected (not built inside
    the class) — see ``testData/scc_132kV_xue/plot_config.py`` for the
    expected shape. Two rendering paths are used depending on graph shape:

    - Plain matrix/scenario plots (series_impedance_matrix, shunt_admittance_matrix,
      earth_return_impedance, earth_return_admittance, earth_return_potential,
      internal_impedance, internal_admittance, internal_potential) are rendered
      by the fully generic ``BasePlotter._plot_matricial_parameters``, driven
      entirely by ``plot_config[graph_key]['data_series']``.
    - "Composition" plots (series_impedance_composition, shunt_admittance_composition,
      potential_coefficients_composition — one (p, q) matrix element decomposed
      into internal + earth-return + series/shunt/potential + assembled-total
      components) are rendered by ``self._plot_composition``, since the "total"
      component is *computed* here (block-diagonal earth-return expansion, or
      parallel admittance combination), not a plain data lookup.

    Expected pul_data:
        pul_data['frequencies']           ndarray (n_freq,)
        pul_data['analytical'] = {'frequencies': ..., 'scenarios': pul_data['scenarios']}
            (alias required by the generic engine; see plot_config.py header comment
            for the synthetic 'internal' scenario needed by internal_* graphs)
        pul_data['comsol']                None, or {'frequencies': ..., 'scenarios':
                                           {'measured': {'impedance_matrix': ndarray (n_freq,2,2)}}}
                                           as returned by ComsolPostProcessor.get_scc_internal_impedance_matrix()
        pul_data['internal_matrices']      from InternalPerUnitParameters.matrices()
        pul_data['scenarios'][scenario_key] = {
            'mtl': MulticonductorTransmissionLine,
            'earth_return_parameters': {'impedance_matrix', 'admittance_matrix', 'potential_coefficient'},
            'quasi_tem_matrices': {'series_impedance_matrix', 'shunt_admittance_matrix', 'potential_coefficient'},
        }

    COMSOL overlay is only physically meaningful for the assembled self/mutual
    impedance matrix (internal_impedance, and the 'internal_term' curve of
    series_impedance_composition) — it cannot resolve decomposed surface
    components, so it is not offered elsewhere.
    """

    def __init__(self, file_path: str, pul_data: dict, plot_config: dict, autoSave: bool = True):
        """
        Args:
            file_path (str): Path of the calling script (typically ``__file__``).
            pul_data (dict): Per-unit-length data — see the class docstring's
                "Expected pul_data" section.
            plot_config (dict): PLOT_CONFIG dict, e.g. from
                ``testData/scc_132kV_xue/plot_config.py``.
            autoSave (bool): If True (default), each plot is saved to Results/.
        """
        super().__init__(file_path, pul_data, plot_config, autoSave=autoSave)
        self.f = pul_data['frequencies']
        self.w = 2 * np.pi * self.f

    # ------------------------------------------------------------------
    # Plain matrix/scenario plots — thin wrappers over the generic engine
    # ------------------------------------------------------------------

    def series_impedance_matrix(self):
        self._plot_matrix_parameters('series_impedance_matrix')

    def shunt_admittance_matrix(self):
        self._plot_matrix_parameters('shunt_admittance_matrix')

    def series_impedance_earth_return(self):
        self._plot_matrix_parameters('earth_return_impedance')

    def shunt_admittance_earth_return(self):
        self._plot_matrix_parameters('earth_return_admittance')

    def potential_coefficients_earth_return(self):
        self._plot_matrix_parameters('earth_return_potential')

    def internal_impedance_matrix(self):
        self._plot_matrix_parameters('internal_impedance')

    def shunt_admittance_internal(self):
        self._plot_matrix_parameters('internal_admittance')

    def potential_coefficients_internal(self):
        self._plot_matrix_parameters('internal_potential')

    # ------------------------------------------------------------------
    # "Composition" plots — dedicated engine (needs to assemble a computed
    # "total" component that BasePlotter._get_data_from_source can't provide)
    # ------------------------------------------------------------------

    def series_impedance_composition(self, conductor='core'):
        self._plot_composition(f'series_impedance_composition_{conductor}')

    def shunt_admittance_composition(self, conductor='core'):
        self._plot_composition(f'shunt_admittance_composition_{conductor}')

    def potential_coefficients_composition(self, conductor='core'):
        self._plot_composition(f'potential_coefficients_composition_{conductor}')

    def _plot_composition(self, graph_key):
        """
        Renders the 'composition' family: one (p, q) matrix element decomposed
        into internal / earth-return / series-or-shunt-or-potential / assembled-
        total components, for a single analytical scenario (``cfg['scenario_key']``).
        """
        cfg = self.plot_config[graph_key]
        p, q = cfg['p'], cfg['q']
        kind = cfg['kind']
        styles = cfg['styles']
        scenario = self.pul_data['scenarios'][cfg['scenario_key']]
        model = scenario['mtl']
        M = model.num_conductors_per_scc
        N = len(self.f)

        if kind == 'series_impedance':
            Zi = self.pul_data['internal_matrices']['impedance_matrix']
            Zg = scenario['earth_return_parameters']['impedance_matrix']
            Zs = scenario['quasi_tem_matrices']['series_impedance_matrix']
            Z0 = np.zeros_like(Zi, dtype=complex)
            for i in range(N):
                Z0[i, :, :] = np.kron(Zg[i, :, :], np.ones((M, M)))
            Zt = Zi + Z0
            components = {
                'internal_term': (Zi[:, p, q], styles['internal_term']),
                'earth_return_term': (Zg[:, 0, 0], styles['earth_return_term']),
                'series_term': (Zs[:, p, q], styles['series_term']),
                'series_composition': (Zt[:, p, q], styles['series_composition']),
            }
        elif kind == 'shunt_admittance':
            Yi = self.pul_data['internal_matrices']['shunt_admittance_matrix']
            Yg = scenario['earth_return_parameters']['admittance_matrix']
            Ysh = scenario['quasi_tem_matrices']['shunt_admittance_matrix']
            Yt = np.zeros_like(Yi, dtype=complex)
            for i in range(N):
                inv_Yi = np.linalg.inv(Yi[i, :, :])
                inv_Yg = np.linalg.inv(Yg[i, :, :]) * np.ones((M, M))
                Yt[i, :, :] = np.linalg.inv(inv_Yi + inv_Yg)
            components = {
                'internal_term': (Yi[:, p, q], styles['internal_term']),
                'earth_return_term': (Yg[:, 0, 0], styles['earth_return_term']),
                'shunt_term': (Ysh[:, p, q], styles['shunt_term']),
                'series_composition': (Yt[:, p, q], styles['series_composition']),
            }
        elif kind == 'potential_coefficient':
            Pi = self.pul_data['internal_matrices']['potential_coefficient_matrix']
            Pg = scenario['earth_return_parameters']['potential_coefficient']
            Psh = scenario['quasi_tem_matrices']['potential_coefficient']
            Pt = Pi + Pg
            components = {
                'internal_term': (np.full(N, Pi[p, q]), styles['internal_term']),
                'earth_return_term': (Pg[:, 0, 0], styles['earth_return_term']),
                'potential_term': (Psh[:, p, q], styles['potential_term']),
                'series_composition': (Pt[:, p, q], styles['series_composition']),
            }
        else:
            raise ValueError(f"Unknown composition kind '{kind}' for graph_key '{graph_key}'.")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg, right_cfg = cfg['left_plot'], cfg['right_plot']

        for data, style in components.values():
            ax1.plot(self.f, self._calculate_plot_data(data, self.w, left_cfg), **style)
            ax2.plot(self.f, self._calculate_plot_data(data, self.w, right_cfg), **style)

        # COMSOL (Js method, core/sheath excitation) only measures the assembled
        # internal impedance matrix, so the overlay is opt-in per graph_key
        # (cfg['comsol']) and restricted to the 'internal_term' component.
        if cfg.get('comsol') and self.cmsl is not None:
            comsol_scn = self.cmsl.get('scenarios', {}).get('measured')
            if comsol_scn is not None:
                Zi_cmsl = comsol_scn['impedance_matrix'][:, p, q]
                f_cmsl = self.cmsl['frequencies']
                w_cmsl = 2 * np.pi * f_cmsl
                scatter_kwargs = dict(marker='o', s=8, facecolors='none', edgecolors='blue', zorder=10, label='COMSOL (Internal)')
                ax1.scatter(f_cmsl, self._calculate_plot_data(Zi_cmsl, w_cmsl, left_cfg), **scatter_kwargs)
                ax2.scatter(f_cmsl, self._calculate_plot_data(Zi_cmsl, w_cmsl, right_cfg), **scatter_kwargs)

        self._format_axis(ax1, self.xlim, left_cfg)
        self._format_axis(ax2, self.xlim, right_cfg)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=graph_key)


class PrysmianCableModels(BasePlotter):
    """
    Plotter for the Prysmian 138 kV core+sheath cable case study: compares the
    exact closed-form Bessel decomposition of the internal impedance
    [Ametani, 2015] against COMSOL (Js method, core/sheath excitation) for the
    assembled self/mutual impedance matrix, plus earth-return impedance across
    several literature formulations. Used by ``scc_138kV_prysmian.py``.

    Config-driven per the project's BasePlotter convention — see
    ``testData/scc_138kV_prysmian/plot_config.py``.

    Expected pul_data:
        pul_data['frequencies']            ndarray (n_freq,)
        pul_data['analytical'] = {'frequencies': ..., 'scenarios': pul_data['scenarios']}
            (alias required by ground_return_impedance(), which uses the
            generic BasePlotter._plot_matricial_parameters engine)
        pul_data['comsol']                 None, or the shape returned by
                                            ComsolPostProcessor.get_scc_internal_impedance_matrix()
        pul_data['internal_parameters']    from InternalPerUnitParameters.parameters_hybrid()
                                            (or .parameters_by_bessel()):
            {'zcs': {'z11', 'z12', 'z2i'}, 'zs3': {'z20', 'z23'}, 'z2m': ...}
        pul_data['scenarios'][scenario_key]['earth_return_parameters']['impedance_matrix']
            (only for ground_return_impedance())

    Note on z2m/z3m oscillation: the closed-form Bessel formula for a tube's
    inner/outer mutual impedance (z2m, z3m) genuinely decays through a rapidly-
    shrinking sign oscillation at high frequency (confirmed against
    exponentially-scaled Bessel functions and arbitrary-precision arithmetic —
    not a numerical artifact). ``internal_impedance_parameters`` auto-floors the
    resistance y-axis for any graph_key whose series include z2m/z3m, just below
    the smallest of the other (non-oscillating) curves, to hide that decorative
    tail without clipping real data.
    """

    def __init__(self, file_path: str, pul_data: dict, plot_config: dict, autoSave: bool = True):
        """
        Args:
            file_path (str): Path of the calling script (typically ``__file__``).
            pul_data (dict): Per-unit-length data — see the class docstring's
                "Expected pul_data" section.
            plot_config (dict): PLOT_CONFIG dict, e.g. from
                ``testData/scc_138kV_prysmian/plot_config.py``.
            autoSave (bool): If True (default), each plot is saved to Results/.
        """
        super().__init__(file_path, pul_data, plot_config, autoSave=autoSave)
        self.f = pul_data['frequencies']
        self.w = 2 * np.pi * self.f

    def ground_return_impedance(self):
        self._plot_matrix_parameters('ground_return_impedance')

    def internal_impedance_parameters(self, graph_key):
        """
        Assembles the named z-component dict for ``graph_key`` ('core', 'sheath',
        'core_sheath', 'internal_parameters') from pul_data['internal_parameters']
        and renders it via ``self._internal_impedance_subplots``.
        """
        data = self.pul_data['internal_parameters']
        impedance_data = {}

        if graph_key == 'core':
            z11 = data['zcs']['z11']
            z12 = data['zcs']['z12']
            z2i = data['zcs']['z2i']
            impedance_data = {'zcs': z11 + z12 + z2i, 'z11': z11, 'z12': z12, 'z2i': z2i}

        elif graph_key == 'sheath':
            z20 = data['zs3']['z20']
            z23 = data['zs3']['z23']
            z2m = data['z2m']
            impedance_data = {'zs3': z20 + z23, 'z20': z20, 'z23': z23, 'z2m': z2m}

        elif graph_key == 'core_sheath':
            zcs = data['zcs']['z11'] + data['zcs']['z12'] + data['zcs']['z2i']
            zs3 = data['zs3']['z20'] + data['zs3']['z23']
            z2m = data['z2m']
            impedance_data = {
                'Zcc': zcs + zs3 - 2 * z2m,
                'Zss': zs3,
                'Zcs': zs3 - z2m,
            }

        elif graph_key == 'internal_parameters':
            impedance_data = {
                'z11': data['zcs']['z11'],
                'z2m': data['z2m'],
                'z2i': data['zcs']['z2i'],
                'z20': data['zs3']['z20'],
            }

        self._internal_impedance_subplots(graph_key, impedance_data)

    def _internal_impedance_subplots(self, graph_key, impedance_data):
        """Generic plotting function for the named z-component impedance data."""
        config = self.plot_config[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series_config in config['series_to_plot']:
            z = impedance_data[series_config['key']]
            style = {k: series_config[k] for k in ['label', 'color', 'linestyle', 'linewidth']}
            ax1.plot(self.f, np.real(z) * 1e3, **style)
            ax2.plot(self.f, np.imag(z) / self.w * 1e6, **style)

        # COMSOL (Js method, core/sheath excitation) only measures the direct
        # self/mutual impedance matrix [Zi] — i.e. Zcc/Zss/Zcs in the 'core_sheath'
        # plot — not the decomposed surface components (z11/z12/z2i/z20/z23/z2m)
        # shown in the 'core'/'sheath'/'internal_parameters' plots, so the overlay
        # is restricted to 'core_sheath'.
        if graph_key == 'core_sheath' and self.cmsl is not None:
            comsol_scn = self.cmsl.get('scenarios', {}).get('measured')
            if comsol_scn is not None:
                Zi_cmsl = comsol_scn['impedance_matrix']
                f_cmsl = self.cmsl['frequencies']
                w_cmsl = 2 * np.pi * f_cmsl
                comsol_index = {'Zcc': (0, 0), 'Zss': (1, 1), 'Zcs': (0, 1)}
                for i, series_config in enumerate(config['series_to_plot']):
                    idx_i, idx_j = comsol_index[series_config['key']]
                    label = 'COMSOL' if i == 0 else None
                    scatter_kwargs = dict(marker='o', s=8, facecolors='none', edgecolors=series_config['color'], zorder=10, label=label)
                    ax1.scatter(f_cmsl, Zi_cmsl[:, idx_i, idx_j].real * 1e3, **scatter_kwargs)
                    ax2.scatter(f_cmsl, Zi_cmsl[:, idx_i, idx_j].imag / w_cmsl * 1e6, **scatter_kwargs)

        # Configure axes
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(self.xlim)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'Resistance $(\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.legend(fontsize='small')
        ax1.set_title(config['resistance_title'])

        # z2m/z3m (mutual impedance between a tube's inner/outer surfaces) decays
        # through a genuine, rapidly shrinking sign-oscillation at high frequency
        # (confirmed against exponentially-scaled Bessel functions and arbitrary-
        # precision arithmetic — not a numerical artifact). On a log axis this
        # produces a long tail of ever-smaller "arcs" that clutters the plot without
        # being physically relevant. Floor the axis just below the smallest of the
        # other (non-oscillating) curves to hide that tail without clipping real data.
        oscillating_keys = {'z2m', 'z3m'}
        plotted_keys = {series_config['key'] for series_config in config['series_to_plot']}
        if plotted_keys & oscillating_keys:
            stable_values = np.concatenate([
                np.abs(np.real(impedance_data[series_config['key']]))
                for series_config in config['series_to_plot']
                if series_config['key'] not in oscillating_keys
            ]) * 1e3
            stable_values = stable_values[stable_values > 0]
            if stable_values.size:
                ax1.set_ylim(bottom=stable_values.min() * 1e-1)

        # Explicit resistance floor override (e.g. 'core_sheath', which doesn't
        # oscillate but still benefits from trimming the very small tail below
        # the physically relevant range).
        if 'resistance_ylim_bottom' in config:
            ax1.set_ylim(bottom=config['resistance_ylim_bottom'])

        ax2.set_xscale('log')
        ax2.set_xlim(self.xlim)
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'Inductance $(mH/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            base_filename = f'impedance_parameters_{graph_key}'
            if 'parameters' in graph_key:
                base_filename = f'impedance_{graph_key}'
            save_figure(fig, self.results_dir, base_filename=base_filename)
