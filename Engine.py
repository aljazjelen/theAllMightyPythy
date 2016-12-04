"""***********************************************************
    Author:     Aljaž Jelen
    Title:
    Date:       December 2016
***********************************************************"""

""" Systems import """
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


class sineEngine():
    """ Class which represents a physical model of engine with sine torque distribution
        This simply means that output torque is fluctuating around certain predetermined,
        static torque amplitude.
        Well described by equation:
                M(t) = M0 + M1*sin(wt)
    """
                                            #TODO create time-dynamic class
    def __init__(self,RPM,T,T_fluctuations,nmbrCylinders):
        """ Method for instantiating sineEngine object
        *_Parameters_*
            :param RPM:             Engine revolutions per minute   [RPM]
            :param T:               Static (mean) torque amplitude  [Nm]
            :param T_fluctuations:  Torque fluctuations (sine)      [Nm]
            :param nmbrCylinders:   Number of cylinders in engine   [/]
        *_Members_*
            :member engineWFreq:    Angular frequency of engine     [rad/s]
            :member drivingWFreq:   Angular frequency of expansions [rad/s]
        """
        self.rpm = RPM
        self.engineWFreq = self.rpm/60*2*np.pi
        self.torque = T
        self.torque_fluctuations = T_fluctuations
        self.cylinders = nmbrCylinders
        self.drivingWFreq = self.engineWFreq * self.cylinders / 2

    def EngineTorque(self, t):
        """ Method for calculating sineEngine output torque
        *_Parameters_*
            :param t:   Time at which total torque is calculated    [s]
        *_Returns_*
           :return:     Total torque of sineEngine                  [Nm]
        """
        return self.torque + self.torque_fluctuations * np.sin(self.drivingWFreq * t)


class realEngine():
    """ Class which represents a physical model of engine with real torque distribution
        Data of output torque is collected from real-life model instead of approximation
        with sine wave.
    """
                #TODO Finish class - maybe add option to read Torque from desired file!
    def __init__(self):
        """

        """