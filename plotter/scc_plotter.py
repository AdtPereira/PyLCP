from .models_base import BasePlotter

class SCCPlotter(BasePlotter):
    def __init__(self, file_path: str, pul_data: dict, plot_config: dict, autoSave: bool = True):
        super().__init__(file_path, pul_data, plot_config, autoSave=autoSave)

    def scc_series_impedance_matrix(self, graph_key_list):
        for key in graph_key_list:
            self._plot_matricial_input(key)

    def scc_shunt_admittance_matrix(self, graph_key_list):
        for key in graph_key_list:
            self._plot_matricial_input(key)

    def scc_earth_return_impedance_matrix(self):
        keys = ['earth_return_impedance_self',
                'earth_return_impedance_mutual_ab',
                'earth_return_impedance_mutual_ac']
        for key in keys:
            self._plot_matricial_input(key)

    def scc_earth_return_admittance_matrix(self):
        keys = ['earth_return_admittance_self',
                'earth_return_admittance_mutual_ab',
                'earth_return_admittance_mutual_ac']
        for key in keys:
            self._plot_matricial_input(key)

    def scc_earth_propagation_constant(self):
        self._plot_non_matricial_parameter('earth_propagation_constant')

class HDPEPlotter(BasePlotter):
    def __init__(self, file_path: str, pul_data: dict, plot_config: dict, autoSave: bool = True):
        super().__init__(file_path, pul_data, plot_config, autoSave=autoSave)

    def hdpe_internal_impedance_matrix(self):
        for key in ['internal_impedance_matrix']:
            self._plot_matricial_upper_triangular(key)

    def hdpe_internal_impedance_elements(self):
        for key in ['coaxial_cable', 'internal_impedance_elements']:
            self._plot_non_matricial_list_parameter(key)
            
