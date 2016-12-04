import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

class DMF():
    def __init__(self,K1,K2,C1,C2,J1,J2):
        self.k1 = K1
        self.k2 = K2
        self.c1 = C1
        self.c2 = C2
        self.j1 = J1
        self.j2 = J2

