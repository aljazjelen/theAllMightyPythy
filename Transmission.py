import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


class Transmission():
    def __init__(self,RPM,inertia,K,C):
        self.TrRPM = RPM
        self.TransWFreq = self.TrRPM/60*2*np.pi
        self.inertia = inertia
        self.TrK = K
        self.TrC = C

    def phiTr(self,t):
        return self.TransWFreq*t