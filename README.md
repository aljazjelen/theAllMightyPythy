# theAllMightyPythy is a repository of all Python projects of mine.


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


PHYSICAL PROPERTIES
Engine:
    engine = sineEngine(RPM, meanTorque, fluctTorque, cyl)

" DoubleMassFlywheel
dmf = DMF(k1, k2, c1, c2, j1, j2)

" Transmission "
trans = Transmission(RPM,jTr,kTr,cTr)

" Wheels "
wheels = Wheels(RPM,jWhe,kWhe,cWhe)

"""     POWERTRAIN ASSEMBLY     """
pwtr = Powertrain(engine,dmf,trans,wheels)
