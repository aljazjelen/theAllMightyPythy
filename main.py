""" System imports """
import matplotlib.pyplot as plt
import numpy as np

""" Mathematical representation of real models imports """
from Powertrain import Powertrain
from Powertrain import sineEngine
from Powertrain import Transmission
from Powertrain import DMF
from Powertrain import Wheels


"""
    As programming classes are consistently created to represent their real/physical counterpart,
    creation of power-train assembly is as easiest as it can get.
    Just set all properties of each model consistently and create Powertrain object with previously defined
    Engine, DMF, Transmission and Wheels.
    All models are with Powertrain class connected together, forming a chain of mathematical models.

    ______________ Powertrain _______________
    |Engine|--|DMF|--|Transmission|--|Wheels|

    In case you want to exclude any of objects in the chain you can do it by assigning value 0 to inertia. As many
    elements can be removed as long as they are last at the chain. Breaking of a chain is not possible!
"""


""" PHYSICAL PROPERTIES """
" Engine "
RPM = 900
meanTorque = 300
fluctTorque = 500
cyl = 6
engine = sineEngine(RPM, meanTorque, fluctTorque, cyl)

" DoubleMassFlywheel "
j1 = 1.8
j2 = 0.6
k1 = 2e4
k2 = 1.1e4
c1 = 30
c2 = 1
dmf = DMF(k1, k2, c1, c2, j1, j2)

" Transmission "
jTr = 0
kTr = 10000000000000
cTr = 0
trans = Transmission(RPM,jTr,kTr,cTr)

jTr = 0.5
kTr = 1000000
cTr = 0
trans2 = Transmission(RPM,jTr,kTr,cTr)

jTr = 0.005
kTr = 1000000
cTr = 0
trans3 = Transmission(RPM,jTr,kTr,cTr)

" Wheels "
jWhe = 0
kWhe = 100000000000
cWhe = 0
wheels = Wheels(RPM,jWhe,kWhe,cWhe)

"""     POWERTRAIN ASSEMBLY     """
pwtr = Powertrain(engine,dmf,trans,wheels)
pwtr2 = Powertrain(engine,dmf,trans2,wheels)
pwtr3 = Powertrain(engine,dmf,trans3,wheels)



""" FREQUENCY RESPONSE """      #TODO    // // // UNDER CONSTRUCTION // // //

" Engine "
RPM = 900
meanTorque = 0
fluctTorque = 500
cyl = 6
engineMean = sineEngine(RPM, meanTorque, fluctTorque, cyl)

" No frequency, for comparison "
wheels_freq = Wheels(RPM,jWhe,kWhe,cWhe)
trans_freq = Transmission(RPM,jTr,kTr,cTr)

#pwtr = Powertrain(engineMean,dmf,trans_freq,wheels_freq)
pwtr = Powertrain(engineMean,dmf,trans,wheels)

