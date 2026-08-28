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
from pathlib import Path

from utils.case_utils import save_figure

__all__ = ["ModalPropagationPlotter", "MODE_STYLE"]

# Fixed colour / dash per mode label, kept consistent across the three figures.
MODE_STYLE = {
    "ground":         {"color": "black",      "linestyle": "-",  "label": "Modo terra"},
    "inter_sheath_1": {"color": "tab:blue",   "linestyle": "--", "label": "Modo entre blindagens 1"},
    "inter_sheath_2": {"color": "tab:blue",   "linestyle": "-.", "label": "Modo entre blindagens 2"},
    "coaxial_1":      {"color": "tab:red",    "linestyle": "--", "label": "Modo coaxial 1"},
    "coaxial_2":      {"color": "tab:green",  "linestyle": "-",  "label": "Modo coaxial 2"},
    "coaxial_3":      {"color": "tab:orange", "linestyle": "-.", "label": "Modo coaxial 3"},
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
    """

    def __init__(self, file_path, modal, config_name="Configuracao 1",
                 alpha_unit="Np/m", autoSave=True):
        self.script_path = Path(file_path)
        self.modal = modal
        self.config_name = config_name
        self.alpha_unit = alpha_unit
        self.autoSave = autoSave

        self.f = np.asarray(modal["frequencies"], dtype=float)
        self.xlim = (self.f.min(), self.f.max())
        self.n = modal["alpha"].shape[1]
        self.labels = modal.get("mode_labels") or tuple(f"mode_{j}" for j in range(self.n))

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
        ax.set_xlabel("Frequencia (Hz)")
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        return fig, ax

    def _finish(self, fig, ax, base_filename):
        ax.legend(fontsize="small")
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=base_filename)
        return fig

    # ------------------------------------------------------------------ #
    def modal_attenuation(self):
        """Fig. 5.5 -- alpha_m(f)."""
        alpha = np.asarray(self.modal["alpha"], dtype=float).copy()
        if self.alpha_unit == "Np/km":
            alpha *= 1e3
        fig, ax = self._new_axis(r"Constante de atenuacao modal $\alpha_m$")
        for j in range(self.n):
            ax.plot(self.f, alpha[:, j], **self._style(self.labels[j], j))
        ax.set_yscale("log")
        ax.set_ylabel(rf"$\alpha_m$ ({self.alpha_unit})")
        return self._finish(fig, ax, "modal_attenuation")

    def modal_phase_velocity(self):
        """Fig. 5.6 -- v_m(f)."""
        v = np.abs(np.asarray(self.modal["vphase"], dtype=float))
        fig, ax = self._new_axis(r"Velocidade de fase modal $v_m$")
        for j in range(self.n):
            ax.plot(self.f, v[:, j], **self._style(self.labels[j], j))
        ax.set_ylabel(r"$v_m$ (m/s)")
        ax.set_ylim(0, min(3.2e8, np.nanpercentile(v, 99) * 1.1))
        return self._finish(fig, ax, "modal_phase_velocity")

    def modal_char_impedance(self):
        """Fig. 5.7 -- |Z_cm(f)|."""
        zc = np.abs(np.asarray(self.modal["Zcm"], dtype=complex))
        fig, ax = self._new_axis(r"Impedancia caracteristica modal $|Z_{cm}|$")
        for j in range(self.n):
            ax.plot(self.f, zc[:, j], **self._style(self.labels[j], j))
        ax.set_yscale("log")
        ax.set_ylabel(r"$|Z_{cm}|$ ($\Omega$)")
        return self._finish(fig, ax, "modal_char_impedance")

    def plot_all(self):
        return [
            self.modal_attenuation(),
            self.modal_phase_velocity(),
            self.modal_char_impedance(),
        ]
