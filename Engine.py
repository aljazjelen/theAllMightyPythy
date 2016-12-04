import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

class sineEngine():
    def __init__(self,RPM,T,T_fluctuations,nmbrCylinders):
        self.rpm = RPM
        self.engineWFreq = self.rpm/60*2*np.pi
        self.torque = T
        self.torque_fluctuations = T_fluctuations
        self.cylinders = nmbrCylinders
        self.drivingWFreq = self.engineWFreq * self.cylinders / 2

    def EngineTorque(self, t):
        return self.torque + self.torque_fluctuations * np.sin(self.drivingWFreq * t)