"""***********************************************************
    Author:     Aljaž Jelen
    Title:
    Date:       December 2016
***********************************************************"""

""" Systems import """
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.linalg import solve
import cmath

""" My custom imports
    needed to assemble Powertrain model
"""
from DMF import DMF
from Engine import sineEngine
from Transmission import Transmission
from Wheels import Wheels

class Powertrain():
    """ Class which represents a physical model of entire powertrain, which is
        assembled from engine, double-mass-flywheel, transmission and wheels
        Class is constructed in such a way, that removing the latest part is as
        easy as assigning inertia of latter to value 0.
        Removing parts in the middle of chain is NOT possible (yet)

        This class is used to recreate real assembly and return various 'rigid body
        dynamics' information of it.
    """
                                        #TODO expand it, a lot more functions and flexibility
    def __init__(self,engine,DMF,transmission,wheels):
        """ Method for instantiating Powertrain object
        *_Parameters_*
            :param engine:          Engine object with all physical data and methods
            :param DMF:             DualMassFlyWheel object with all physical data and methods
            :param transmission:    Transmission object with all physical data and methods
            :param wheels:          Wheels object with all physical data and methods
        *_Members_*
            :member initData:       Initial conditions array of vibrations  [position, velocity]
            :member M:              Inertia matrix                          [m^2*kg]
            :member K:              Stiffness matrix                        [N/rad]
            :member C:              Damping matrix                          [N/rad*s]
            :member Mr:             Matrix for evaluating frequency domain  [N/rad*s]
            :member timeProperties: Simulation properties of time domain    [startTime, endTime, timeStep]
            :member timeDomain:     Time domain array with specified step   [s]
            :member derivSolution:  Solution of ODEs in time domain (odeint)[angle1, velocity1,..., angleN, velocityN]
            :member freqSolution:   Solution of ODEs in frequency domain    [s]
            :member calcMatrices:   Sets M, K and C member matrices
            :member setInitials:    Extracts and sets initial conditions
            :member solveTimeDomain:Solves system of equation in time domain with 'odeint'
        """
        self.dmf = DMF
        self.engine = engine
        self.trans = transmission
        self.wheels = wheels

        self.initData = np.array(())    # Contains initial data according to objects used - #TODO check setInitials()
        self.M = np.array(())
        self.K = np.array(())
        self.C = np.array(())

        self.Mfr = np.array(())
        self.Kfr = np.array(())
        self.Cfr = np.array(())
        self.Ffr = np.array(())

        self.timeProperties = {'startTime': 0.0, 'endTime':10,'timeStep':0.0001}
        self.timeDomain = np.arange(self.timeProperties['startTime'], self.timeProperties['endTime'], self.timeProperties['timeStep'])

        self.derivSolution = []
        self.freqSolution = []
        self.calcMatrices()
        self.setInitials()
        self.solveTimeDomain()


    def calcMatrices(self):
        """ Method to set M, K and C matrices for calculations and more robust ODE handling. Matrices coefficients are
            set according to properties of Engine, Transmission, DMF, Wheels. Matrices are a size of N,
            where N is number of elements in chain of mathematical model. Each sub-part of model (fe Wheels)
            is easily neglected by setting its inertia to 0. Neglected can only be models at the very end of the chain.
        *_Parameters_*
            NONE
        *_Returns_*
            :return:    NONE
        """
        if (self.trans.inertia is 0):
            print("Transmission has no inertia, therefore will not be taken into account...")
            self.M = np.array(((self.dmf.j1, 0), (0, self.dmf.j2)))           # This is N x N inertia matrix
            self.K = np.array(((self.dmf.k1, -self.dmf.k1), (-self.dmf.k1, self.dmf.k1 + self.dmf.k2)))  # This is N x N stiffness matrix
            self.C = np.array(((self.dmf.c1, -self.dmf.c1), (-self.dmf.c1, self.dmf.c1 + self.dmf.c2)))  # This is N x N damping matrix
            #self.F = np.array(((0), (0), (engine.EngineTorque(t)), (c2 * trans.TransWFreq + k2 * trans.phiTr(t))))
        else:
            self.M = np.array(((self.dmf.j1, 0, 0), (0, self.dmf.j2,0), (0, 0, self.trans.inertia)))  # This is N x N inertia matrix
            self.K = np.array(((self.dmf.k1, -self.dmf.k1, 0),(-self.dmf.k1, self.dmf.k1 + self.dmf.k2, -self.dmf.k2),(0, -self.dmf.k2, self.dmf.k2 + self.trans.TrK)))  # This is N x N stiffness matrix
            self.C = np.array(((self.dmf.c1, -self.dmf.c1, 0),(-self.dmf.c1, self.dmf.c1 + self.dmf.c2, -self.dmf.c2),(0, -self.dmf.c2, self.dmf.c2 + self.trans.TrC)))  # This is N x N damping matrix
            if (self.wheels.inertia is 0):
                print("Wheels have no inertia, therefore will not be taken into account...")
                print("Not yet implemented, derive equations first!")
                #self.M = np.array(((self.dmf.j1, 0), (0, self.dmf.j2)))  # This is N x N inertia matrix
                #self.K = np.array(((self.dmf.k1, -self.dmf.k1),(-self.dmf.k1, self.dmf.k1 + self.dmf.k2)))  # This is N x N stiffness matrix
                #self.C = np.array(((self.dmf.c1, -self.dmf.c1),(-self.dmf.c1, self.dmf.c1 + self.dmf.c2)))  # This is N x N damping matrix
                # self.F = np.array(((0), (0), (engine.EngineTorque(t)), (c2 * trans.TransWFreq + k2 * trans.phiTr(t))))
            else:
                print("Not yet implemented, derive equations first!")

    def calcFrequencyMatrices(self,omega):
        """ Method to set matrix Mfr for calculating behaviour in frequency domain
            according to properties of Engine, Transmission, DMF, Wheels
        *_Parameters_*
            :param omega:       Frequency at which solution in frequency domain is desired  [rad/s]
        *_Returns_*
        :return:    NONE
        """
        print("Calculating Matrices for frequency: " + str(omega))
        self.Mfr = np.array(((-self.dmf.j1*omega*omega+self.dmf.k1+omega*self.dmf.c1*1j, -self.dmf.k1-1j*omega*self.dmf.c1), (-self.dmf.k1-1j*omega*self.dmf.c1, -self.dmf.j2+self.dmf.k1+self.dmf.k2+1j*omega*self.dmf.c1+1j*omega*self.dmf.c2)))  # This is N x N inertia matrix
        self.Kfr = np.array(
            ((self.dmf.k1, -self.dmf.k1), (-self.dmf.k1, self.dmf.k1 + self.dmf.k2)))  # This is N x N stiffness matrix
        self.Cfr = np.array(
            ((self.dmf.c1, -self.dmf.c1), (-self.dmf.c1, self.dmf.c1 + self.dmf.c2)))  # This is N x N damping matrix
        #self.Ffr = np.array(((self.engine.torque),(self.dmf.k2+1j*self.dmf.c2*omega)))
        #self.Ffr = np.array(((self.engine.torque_fluctuations), (0)))
        self.Ffr = np.array(((1j*self.engine.torque_fluctuations), (0)))



    def calcForFrequency(self,omega):
        """ Method to calculate complex constants for each degree of freedom - \theta which describes
            model's behaviour at given frequency and time
        *_Parameters_*
            :param omega:       Frequency at which solution in frequency domain is desired  [rad/s]
        *_Returns_*
            :return:    NONE
        """
        self.calcFrequencyMatrices(omega)
        print("[Powertrain] First element in frequency matrix is: " + str(self.Mfr[0,0]))
        self.freqSolution = solve(self.Mfr, self.Ffr)

    def setInitials(self,phi1=0,phi2=0,phiTr=0):
        """ Method for setting initial conditions of simulation
            Each additional element in chain (with inertia greater than 0) introduces 2 new initial conditions
            due to order of ODE - 2nd order
        *_Parameters_*
            :param phi1:    Initial angle of first plate in DualMassFlywheel model  [rad]
            :param phi2:    Initial angle of second plate in DualMassFlywheel model [rad]
            :param phiTr:   Initial angle of transmission model                     [rad]
        *_Members_*
            :member initData:       Initial conditions array of vibrations          [position1, velocity1,... positionN]
        :return:
        """
                                            #TODO - make a use out of arguments! for now they have no purpose
        phi1_init = phi1
        phi2_init = phi2
        dphi1_init = self.engine.engineWFreq
        dphi2_init = self.trans.TransWFreq
        self.initData = np.array(((phi1_init), (phi2_init), (dphi1_init), (dphi2_init)))

        # Checks if transmission is present and add it to initial conditions if so
        if (self.trans.inertia is not 0):
            print("Initializing: Transmission inertia is not set to 0, therefore Wheels frequency is taken as last!")
            dphiTr_init = self.wheels.WheWFreq
            phiTr_init = phiTr
            self.initData = np.array(((phi1_init), (phi2_init), (phiTr_init), (dphi1_init), (dphi2_init), (dphiTr_init)))
            if (self.wheels.inertia is not 0):  # TODO add floor variables (car comparing to floor)
                print("Initializing: Wheels inertia is not set to 0, therefore Floor-Wheel-Car frequency is taken as last!")
                print("Not yet implemented, derive equations first!")



    def solveTimeDomain(self):
        """ Method to solve a given system of equations and generate solutions which are stored in member derivSolution
        *_Parameters_*
            NONE
        :return:    NONE
        """
        self.derivSolution = odeint(self.deriv, self.initData, self.timeDomain, rtol=1e-6,
                                    atol=1e-6)  # TODO make some simulation method?

    def deriv(self,u, t):
        """ Method which consistently combines M, C and K matrices into system of A and B with external force F. It is
            broken down into 2 sets which together form a linear system of equations. Each set contains N number of
            equations, where N is number of objects present in mathematical model. For more information please check
            chapter "Programming model".
            This system is then pushed into solving stage where 'odeint' method handles it through specified time
            domain and given tolerances.


            We are solving system of equations: y' = A*y(t) + B*f(t)
            where we, by 'odeint' obtain vector u, which contains vector y and vector y'

            SciPy official example of refracting 2nd order non-homogenous differential equations into 2 subsets:
            https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.integrate.odeint.html

        *_Parameters_*
            :param u:   Previous solution pushed back into method 'odeint' for solving and obtaining next solution in t
            :param t:   Point of time in time domain where 'odeint' is currently solving
        :return:        Solution vector at given time t, with size 2xtime (2 since we've got angle and velocity)
        """
        invM = np.linalg.inv(self.M)            # create inverse of Inertia matrix of less notation later on

        matrix_size = len(self.M)               # To check how many masses are we dealing with
        a_ident = np.identity(matrix_size)
        a_zeros = np.zeros((matrix_size, matrix_size))
        a_upper = np.hstack((a_zeros, a_ident))
        a_lower = np.hstack((np.dot(-invM, self.K), np.dot(-invM, self.C)))
        A = np.vstack((a_upper, a_lower))

        # b_upper = np.array(((0, 0), (0, 0)))
        # b_lower = invM
        # b_upper = np.array(((0, 0, 0, 0), (0, 0, 0, 0)))
        b_upper = np.hstack((a_zeros, a_zeros))
        b_lower = np.hstack((a_zeros, invM))
        B = np.vstack((b_upper, b_lower))

        " Obtain stiffness and damping of last object, to be able to calculate torque in vector of external Forces "
        end_stiffness = self.K[matrix_size - 1, matrix_size - 2] + self.K[matrix_size - 1, matrix_size - 1]
        end_damping = self.C[matrix_size - 1, matrix_size - 2] + self.C[matrix_size - 1, matrix_size - 1]

        " External Forces Vector - known variables (External forces, end Torque)"
        # F = np.array(((engine.EngineTorque(t)), (c2 * trans.TransWFreq + k2 * trans.phiTr(t))))
        F_dphi = np.zeros(matrix_size)
        F_ddphi = np.zeros(matrix_size)
        F_ddphi[0] = self.engine.EngineTorque(t)
        endWfreq = self.trans.TransWFreq
        endPhi = self.trans.phiTr(t)

        if (self.trans.inertia is not 0):
            endWfreq = self.wheels.WheWFreq
            endPhi = self.wheels.phiWhe(t)
            if (self.wheels.inertia is not 0):           # TODO add floor variables (car comparing to floor)
                ndWfreq = self.wheels.WheWFreq
                endPhi = self.wheels.phiWhe(t)

        F_ddphi[matrix_size - 1] = end_damping * endWfreq + end_stiffness * endPhi
        F = np.hstack((F_dphi, F_ddphi))
        F.transpose()

        solutionVector = np.dot(A, u) + np.dot(B, F)
        " Vector of solutions stores free coordinates and their first derivatives "
        " in order:     solutionVector[0 to i]      : angles (non derivatives)"
        "               solutionVector[i+1 to end]  : angular velocity (1st derivative)"
        " This vector is continuously being pushed into OdeInt method to integrate throughout"
        " entire time domain "
        return solutionVector
