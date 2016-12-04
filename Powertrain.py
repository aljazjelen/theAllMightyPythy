" System imports "
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.linalg import solve
import cmath

" My imports "
from DMF import *
from Engine import sineEngine
from Transmission import Transmission
from Wheels import Wheels

class Powertrain():
    def __init__(self,engine,DMF,transmission,wheels):
        self.dmf = DMF
        self.engine = engine
        self.trans = transmission
        self.wheels = wheels
        self.initData = np.array(())    # Contains initial data according to objects used - check setInitials()
        self.M = np.array(())
        self.K = np.array(())
        self.C = np.array(())

        self.Mfr = np.array(())
        self.Kfr = np.array(())
        self.Cfr = np.array(())
        self.Ffr = np.array(())

        """ TIME DOMAIN """
        self.timeProperties = {'startTime': 0.0, 'endTime':10,'timeStep':0.0001}
        self.timeDomain = np.arange(self.timeProperties['startTime'], self.timeProperties['endTime'], self.timeProperties['timeStep'])

        self.derivSolution = []
        self.freqSolution = []
        self.calcMatrices()             # Calculate stiffness, damping, mass matrices
        self.setInitials()              # Sets most ordinary initial conditions


    def calcMatrices(self):
        """ MATRICES OF COEFFICIENTS """
        " Change matrices to expand the model "
        " Say we have N number of masses "
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
        self.calcFrequencyMatrices(omega)
        print("[Powertrain] First element in frequency matrix is: " + str(self.Mfr[0,0]))
        self.freqSolution = solve(self.Mfr, self.Ffr)
        x = solve(self.Mfr, self.Ffr)
        #print(x)

        X = cmath.polar(x[1])
        #print(cmath.polar(x[1]))

        #print(x[1]*cmath.exp(1j*omega*1))


    def setInitials(self,phi1=0,phi2=0,phiTr=0):
        """ INITIAL CONDITIONS """
        " Are subject to refracting 2nd order ODE to 1st order ODEs and assigning initial conditions "
        " to all substitutes of 2nd derivative (angular/linear acceleration)"
        phi1_init = 0
        phi2_init = 0

        dphi1_init = self.engine.engineWFreq
        dphi2_init = self.trans.TransWFreq

        self.initData = np.array(((phi1_init), (phi2_init), (dphi1_init), (dphi2_init)))

        if (self.trans.inertia is not 0):
            print("Initializing: Transmission inertia is not set to 0, therefore Wheels frequency is taken as last!")
            dphiTr_init = self.wheels.WheWFreq
            phiTr_init = 0
            self.initData = np.array(((phi1_init), (phi2_init), (phiTr_init), (dphi1_init), (dphi2_init), (dphiTr_init)))
            if (self.wheels.inertia is not 0):  # TODO add floor variables (car comparing to floor)
                print("Initializing: Wheels inertia is not set to 0, therefore Floor-Wheel-Car frequency is taken as last!")
                print("Not yet implemented, derive equations first!")

        #time = np.arange(self.timeDomain['startTime'], self.timeDomain['endTime'], self.timeDomain['timeStep'])
        print(self.M)
        self.derivSolution = odeint(self.deriv, self.initData, self.timeDomain, rtol=1e-6, atol=1e-6)  #TODO make some simulation method?


    def deriv(self,u, t):
        invM = np.linalg.inv(self.M)

        " We are solving system of equations: y' = A*y(t) + B*f(t) "
        " where we, by odeint obtain vector u, which contains vector y and vector y' "
        matrix_size = len(self.M)  # To check how many masses are we dealing with
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

        " Obtain stiffness and damping of last object, to be able to calculate torque in "
        " Forces vector"
        end_stiffness = self.K[matrix_size - 1, matrix_size - 2] + self.K[matrix_size - 1, matrix_size - 1]
        end_damping = self.C[matrix_size - 1, matrix_size - 2] + self.C[matrix_size - 1, matrix_size - 1]

        " Forces Vector - known variables (External forces, end Torque"
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
        # F = np.vstack((np.transpose(F_dphi),np.transpose(F_ddphi)))
        # F = np.array(((0), (0), (engine.EngineTorque(t)), (end_damping * trans.TransWFreq + end_stiffness * trans.phiTr(t))))
        " Vector of solutions stores free coordinates and their first derivatives "
        " in order:     solutionVector[0 to i]      : angles (non derivatives)"
        "               solutionVector[i+1 to end]  : angular velocity (1st derivative)"
        """
        print(self.M)
        print(self.K)
        print(self.C)
        print(A)
        print(B)
        print(F)
        """
        solutionVector = np.dot(A, u) + np.dot(B, F)
        " This vector is continuously being pushed into OdeInt method to integrate throughout"
        " entire time domain "
        return solutionVector
