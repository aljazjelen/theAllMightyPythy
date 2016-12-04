"""***********************************************************
    Author:     Aljaž Jelen
    Title:
    Date:       December 2016
***********************************************************"""

""" Systems import """
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


class DMF():
    """ Class which represents a physical model of wheels   """
                                            #TODO create time-dynamic class
    def __init__(self,K1,K2,C1,C2,J1,J2):
        """ Method for instantiating DualMassFlywheel object
        *_Parameters_*
            :param K1:      Stiffness of first plate    [N/rad]
            :param K2:      Stiffness of second plate   [N/rad]
            :param C1:      Damping of first plate      [N/rad*s]
            :param C2:      Damping of second plate     [N/rad*s]
            :param J1:      Inertia of first plate      [m^2*kg]
            :param J2:      Inertia of second plate     [m^2*kg]
        """
        self.k1 = K1
        self.k2 = K2
        self.c1 = C1
        self.c2 = C2
        self.j1 = J1
        self.j2 = J2

