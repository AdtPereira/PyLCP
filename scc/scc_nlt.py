import numpy as np
from scc_parameters import TransmissionLineParameters
from scc_functions import VectorialOperation  

class MonoNetworkTopology():

    def __init__(self, s, syst, RS, LX):
        # Create an instance of the TransmissionLineParameters class
        self.s = s
        self.syst = syst
        self.RS = RS
        self.LX = LX
        self.tl = TransmissionLineParameters(s, syst)

    def StepUnitSource(self, ksi, T):
        self.VS = (1/self.s) * np.exp(-ksi*T * self.s)

    def TerminalVoltages(self, VS):
        # Calls the methods of the TransmissionLineParameters class 
        # to calculate the quadrupole parameters
        self.tl.SeriesImpedance()
        self.tl.ShuntAdmittance()
        self.tl.PropagationFunction()

        # Calls the QuadripoleParameters method to calculate Ybus
        self.Ybus = self.tl.QuadripoleParameters(self.RS, self.LX)
            
        # V = np.zeros((self.syst.nc * NUM_BARS, 1))        
        self.V = np.linalg.solve(self.Ybus, np.array([[VS/self.RS], [0]]))

class NumericalLaplaceTransform():

    def __init__(self, N, T):
        self.N = N
        self.dt = T/N
        self.c = np.log(N)/T
        self.m = np.arange(N)
        self.t = self.dt*self.m       
        self.dw = 2 * np.pi/T
        self.s = self.c + 1j* self.m * self.dw 
        self.f = np.imag(self.s) / 2 / np.pi 

    def Main_NLT(self, Vk, Vm):
        self.Cn = np.exp(self.c * self.m * self.dt) / self.dt
        self.H = 0.5 * (1 + np.cos(2 * np.pi * self.m / self.N))
        
        Vk = [x * y for x, y in zip(Vk, self.H)]
        Vm = [x * y for x, y in zip(Vm, self.H)]
        
        Vktd = [x * y for x, y in zip(self.Cn, np.real(np.fft.ifft(Vk)))] 
        Vmtd = [x * y for x, y in zip(self.Cn, np.real(np.fft.ifft(Vm)))] 

        return Vktd, Vmtd