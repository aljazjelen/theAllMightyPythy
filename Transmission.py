"""***********************************************************
    Author:     Aljaž Jelen
    Title:
    Date:       December 2016
***********************************************************"""

""" Systems import """
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


class Transmission():
    """ Class which represents a physical model of transmission """
                                            #TODO create time-dynamic class
    def __init__(self,RPM,inertia,K,C):
        """ Method for instantiating Transmission object
        *_Parameters_*
            :param RPM:     Expected input revolutions per minute       [RPM]
            :param inertia: Approximated inertia of spurs and shafts    [m^2*kg]
            :param K:       Approximated stiffness                      [N/rad]
            :param C:       Approximated damping coefficient            [N/rad*s]
        *_Members_*
            :member TransWFreq:    Angular frequency of transmission    [rad/s]
        """
        self.TrRPM = RPM
        self.TransWFreq = self.TrRPM/60*2*np.pi
        self.inertia = inertia
        self.TrK = K
        self.TrC = C

    def phiTr(self,t):
        """ Method for calculating angle done by transmission input shaft
        *_Parameters_*
            :param t:       Time at which angle is calculated   [s]
        *_Returns_*
            :return:        Angle at specified time             [rad]
        """
        return self.TransWFreq*t