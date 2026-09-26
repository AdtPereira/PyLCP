"""Modal-domain propagation plots -- Chapter 5 of Andreata's thesis.

Reproduces, for **Configuration 1** (three directly-buried single-core
cables), the three per-configuration figures:

* Fig. 5.5  -- modal attenuation constant  alpha_m(f)      [Np/m]
* Fig. 5.6  -- modal phase velocity         v_m(f)          [m/s]
* Fig. 5.7  -- modal characteristic impedance |Z_cm(f)|     [Ohm]

Input is the dict returned by
``analytical_forms.modal_analysis.ModalDecomposition.modal_parameters``.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

from utils.case_utils import save_figure

__all__ = ["ModalPropagationPlotter", "MODE_STYLE", "print_modal_comparison_report"]

# Fixed colour / dash per mode label, kept consistent across the three figures.
MODE_STYLE = {
    "ground":         {"color": "black",      "linestyle": "-",  "label": "Ground mode"},
    "inter_sheath_1": {"color": "tab:blue",   "linestyle": "--", "label": "Inter-sheath mode 1"},
    "inter_sheath_2": {"color": "tab:blue",   "linestyle": "-.", "label": "Inter-sheath mode 2"},
    "coaxial_1":      {"color": "tab:red",    "linestyle": "--", "label": "Coaxial mode 1"},
    "coaxial_2":      {"color": "tab:green",  "linestyle": "-",  "label": "Coaxial mode 2"},
    "coaxial_3":      {"color": "tab:orange", "linestyle": "-.", "label": "Coaxial mode 3"},
    "ecc":            {"color": "tab:purple", "linestyle": ":",  "label": "ECC mode"},
}


class ModalPropagationPlotter:
    """
    Parameters
    ----------
    file_path : str
        ``__file__`` of the calling case script (used to locate ``Results/``).
    modal : dict
        Output of ``ModalDecomposition.modal_parameters`` -- must contain
        ``frequencies``, ``alpha``, ``vphase``, ``Zcm`` (each ``(Nf, N)``) and
        ``mode_labels`` (tuple of ``N`` strings, or ``None``).
    config_name : str
        Short label for the configuration (goes in the sub-title).
    alpha_unit : {'Np/m', 'Np/km'}
    autoSave : bool
    matlab_modal : dict, optional
        MATLAB reference, as returned by
        ``utils.matlab_data.MatlabDataReader.get_modal_scenario_data`` --
        ``frequencies``, ``mode_labels`` and per-quantity ``(Nf, 6)`` arrays
        (``alpha``, ``vphase``, ``Zcm``, ...). When given, each figure
        overlays the matching mode (matched by label, not by column index --
        the two decompositions do not share a frequency grid or an
        eigenvector ordering) as open-circle markers. ``None`` (default)
        disables the overlay.
    """

    def __init__(self, file_path, modal, config_name="Configuration 1",
                 alpha_unit="Np/m", autoSave=True, matlab_modal=None):
        self.script_path = Path(file_path)
        self.modal = modal
        self.config_name = config_name
        self.alpha_unit = alpha_unit
        self.autoSave = autoSave
        self.matlab_modal = matlab_modal

        self.f = np.asarray(modal["frequencies"], dtype=float)
        self.xlim = (self.f.min(), self.f.max())
        self.n = modal["alpha"].shape[1]
        self.labels = modal.get("mode_labels") or tuple(f"mode_{j}" for j in range(self.n))
        # plotting (and legend) order follows MODE_STYLE, not the tracked column order
        style_order = list(MODE_STYLE)
        self.order = sorted(range(self.n), key=lambda j: (
            style_order.index(self.labels[j]) if self.labels[j] in style_order else len(style_order), j))

        self.results_dir = os.path.join("testData", self.script_path.stem, "Results")
        os.makedirs(self.results_dir, exist_ok=True)

    # ------------------------------------------------------------------ #
    def _style(self, label, j):
        base = MODE_STYLE.get(label, {"linestyle": "-", "label": label})
        sty = {"linewidth": 1.4}
        sty.update({k: v for k, v in base.items() if k != "label"})
        sty["label"] = base.get("label", label)
        if "color" not in sty:
            sty["color"] = f"C{j}"
        return sty

    def _new_axis(self, title):
        fig, ax = plt.subplots(figsize=(7.0, 5.0))
        fig.suptitle(f"{title}\n{self.config_name}", fontsize=11)
        ax.set_xscale("log")
        ax.set_xlim(self.xlim)
        ax.set_xlabel("Frequency (Hz)")
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        return fig, ax

    def _finish(self, fig, ax, base_filename):
        handles, labels = ax.get_legend_handles_labels()
        if self.matlab_modal is not None:
            proxy = Line2D([], [], marker='o', linestyle='None', markersize=6,
                            markerfacecolor='none', markeredgecolor='black',
                            label='MATLAB (Andreata)')
            handles, labels = handles + [proxy], labels + ['MATLAB (Andreata)']
        ax.legend(handles, labels, fontsize="small")
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=base_filename)
        return fig

    def _overlay_matlab(self, ax, label, j, key, transform=None):
        """Scatters the MATLAB reference column matching mode `label` (looked
        up by name in ``self.matlab_modal['mode_labels']`` -- the MATLAB
        modal decomposition uses its own eigenvector ordering, so column `j`
        of `self.modal` has no reason to line up with column `j` of the
        MATLAB arrays) on top of the just-plotted pyLCP curve. No-op if
        `self.matlab_modal` is unset or does not carry `key`/`label`."""
        mtlb = self.matlab_modal
        if mtlb is None or key not in mtlb:
            return
        mtlb_labels = mtlb.get("mode_labels") or ()
        if label not in mtlb_labels:
            return
        idx = mtlb_labels.index(label)
        values = mtlb[key][:, idx]
        if transform is not None:
            values = transform(values)
        color = self._style(label, j).get("color", f"C{j}")
        ax.scatter(mtlb["frequencies"], values, marker='o', s=10,
                   facecolors='none', edgecolors=color, linewidths=0.9, zorder=9)

    # ------------------------------------------------------------------ #
    def modal_attenuation(self):
        """Fig. 5.5 -- alpha_m(f)."""
        alpha = np.asarray(self.modal["alpha"], dtype=float).copy()
        unit_scale = 1e3 if self.alpha_unit == "Np/km" else 1.0
        alpha *= unit_scale
        fig, ax = self._new_axis(r"Modal attenuation constant $\alpha_m$")
        for j in self.order:
            ax.plot(self.f, alpha[:, j], **self._style(self.labels[j], j))
            self._overlay_matlab(ax, self.labels[j], j, "alpha",
                                 transform=lambda v: v * unit_scale)
        ax.set_yscale("log")
        ax.set_ylabel(rf"$\alpha_m$ ({self.alpha_unit})")
        return self._finish(fig, ax, "modal_attenuation")

    def modal_phase_velocity(self):
        """Fig. 5.6 -- v_m(f)."""
        v = np.abs(np.asarray(self.modal["vphase"], dtype=float))
        fig, ax = self._new_axis(r"Modal phase velocity $v_m$")
        for j in self.order:
            ax.plot(self.f, v[:, j], **self._style(self.labels[j], j))
            self._overlay_matlab(ax, self.labels[j], j, "vphase", transform=np.abs)
        ax.set_ylabel(r"$v_m$ (m/s)")
        v_top = np.nanpercentile(v, 99)
        if self.matlab_modal is not None and "vphase" in self.matlab_modal:
            v_top = max(v_top, np.nanpercentile(np.abs(self.matlab_modal["vphase"]), 99))
        ax.set_ylim(0, min(3.2e8, v_top * 1.1))
        return self._finish(fig, ax, "modal_phase_velocity")

    def modal_char_impedance(self):
        """Fig. 5.7 -- |Z_cm(f)|."""
        zc = np.abs(np.asarray(self.modal["Zcm"], dtype=complex))
        fig, ax = self._new_axis(r"Modal characteristic impedance $|Z_{cm}|$")
        for j in self.order:
            ax.plot(self.f, zc[:, j], **self._style(self.labels[j], j))
            self._overlay_matlab(ax, self.labels[j], j, "Zcm", transform=np.abs)
        ax.set_yscale("log")
        ax.set_ylabel(r"$|Z_{cm}|$ ($\Omega$)")
        return self._finish(fig, ax, "modal_char_impedance")

    def plot_all(self):
        return [
            self.modal_attenuation(),
            self.modal_phase_velocity(),
            self.modal_char_impedance(),
        ]


def print_modal_comparison_report(modal, matlab_modal, target_freqs_hz=(1e2, 1e4, 1e6),
                                  title="Modal parameters -- pyLCP vs. MATLAB (Andreata)"):
    """Prints, per mode and per target frequency, the relative error of
    pyLCP's modal parameters against the external MATLAB reference.

    The two decompositions do not share a frequency grid (pyLCP: the case's
    own sweep; MATLAB: the "standard" 91-point grid, see
    ``MatlabDataReader.get_modal_scenario_data``), so each target frequency is
    resolved independently on each grid to its nearest sample -- adequate for
    a handful of spot-check frequencies on a 10-points/decade log grid, not a
    substitute for the visual overlay in ``ModalPropagationPlotter``.

    Parameters
    ----------
    modal : dict
        Output of ``ModalDecomposition.modal_parameters`` (needs
        ``frequencies``, ``mode_labels``, ``alpha``, ``vphase``, ``Zcm``).
    matlab_modal : dict or None
        Output of ``MatlabDataReader.get_modal_scenario_data``. If ``None``,
        prints a one-line notice and returns.
    target_freqs_hz : tuple of float
        Frequencies (Hz) at which the comparison is reported.
    """
    line = "=" * 88
    print("\n" + line)
    print(f"  {title}")
    print(line)
    if matlab_modal is None:
        print("  no MATLAB modal reference available -- skipping.")
        print(line + "\n")
        return

    f_pylcp = np.asarray(modal["frequencies"], dtype=float)
    f_mtlb = np.asarray(matlab_modal["frequencies"], dtype=float)
    labels = modal.get("mode_labels")
    if labels is None:
        print("  pyLCP modes are unlabelled (mode_labels=None) -- skipping.")
        print(line + "\n")
        return

    quantities = (
        ("alpha_m", "alpha", "alpha", np.asarray, "Np/m"),
        ("v_m", "vphase", "vphase", np.abs, "m/s"),
        ("|Z_cm|", "Zcm", "Zcm", np.abs, "Ohm"),
    )

    mtlb_labels = matlab_modal.get("mode_labels") or ()
    for label in labels:
        if label not in mtlb_labels:
            print(f"  {label:16s}  -- no MATLAB reference for this mode.")
            continue
        j = labels.index(label)
        idx_m = mtlb_labels.index(label)
        print(f"  -- mode: {label} --")
        for disp_name, pylcp_key, mtlb_key, transform, unit in quantities:
            row = []
            for f_target in target_freqs_hz:
                kp = int(np.argmin(np.abs(f_pylcp - f_target)))
                km = int(np.argmin(np.abs(f_mtlb - f_target)))
                v_pylcp = transform(modal[pylcp_key][kp, j])
                v_mtlb = transform(matlab_modal[mtlb_key][km, idx_m])
                denom = abs(v_mtlb) if abs(v_mtlb) > 0 else 1.0
                rel = abs(v_pylcp - v_mtlb) / denom
                row.append(f"{f_target:8.0e} Hz: pyLCP={v_pylcp:.4g} "
                          f"MATLAB={v_mtlb:.4g} ({rel:.2%})")
            print(f"    {disp_name:8s} [{unit:5s}]  " + "  |  ".join(row))
    print(line + "\n")
