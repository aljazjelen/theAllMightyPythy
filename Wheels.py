"""***********************************************************
    Author:     Aljaž Jelen
    Title:
    Date:       December 2016
***********************************************************"""

""" Systems import """
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


class Wheels():
    """ Class which represents a physical model of wheels   """
                                            #TODO create time-dynamic class
    def __init__(self,RPM,inertia,K,C):
        """ Method for instantiating Wheels object
        *_Parameters_*
            :param RPM:     Expected input revolutions per minute           [RPM]
            :param inertia: Approximated inertia of all wheels              [m^2*kg]
            :param K:       Approximated stiffness of all wheels            [N/rad]
            :param C:       Approximated damping coefficient of all wheels  [N/rad*s]
        *_Members_*
            :member WheWFreq:    Angular frequency of wheels                [rad/s]
        """
        self.WheRPM = RPM
        self.WheWFreq = self.WheRPM/60*2*np.pi
        self.inertia = inertia
        self.WheK = K
        self.WheC = C

    def phiWhe(self,t):
        """ Method for calculating angle done by wheels
        *_Parameters_*
            :param t:       Time at which angle is calculated   [s]
        *_Returns_*
            :return:        Angle at specified time             [rad]
        """
        return self.WheWFreq*t