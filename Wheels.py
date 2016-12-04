import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


class Wheels():
    def __init__(self,RPM,inertia,K,C):
        self.WheRPM = RPM
        self.WheWFreq = self.WheRPM/60*2*np.pi
        self.inertia = inertia
        self.WheK = K
        self.WheC = C

    def phiWhe(self,t):
        return self.WheWFreq*t