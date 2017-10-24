# here is a file with configurations: initial conditions and functions, which are mentioned in manuals for my variant

import numpy

isDebug = True
DEFAULT_EPSILON = 1e-12

# differential
DIFFERENCIATION_FUNCTION = lambda x: numpy.cos(x**4 + 5*x**3 + 3*x + 11)
DIFFERENCIATION_POINT = (numpy.sqrt(3) - 9) / 5

# integration
INTEGRATION_FUNCTION = lambda x: numpy.sin(1 / (x**2 + 1))
INTEGRATION_BOUNDARIES = {"a": 0.5,
                          "b": 0.8}
# improper integration
IMPROPER_INTEGRATION_FUNCTION = lambda x: x/((numpy.abs(1 - x**2)) ** (0.2))
IMPROPER_INTEGRATION_BOUNDARIES = {"a": 0.3,
                                   "b": 1}

# system of differential equation
#DIFFERENTIAL_EQUATION_CONDITIONS = {"a" : 2.09, "b" : 2.45, "c": 1.65, "d" : 1.25, "x0" : 3.73, "y0" : 2.55, "t0" : 21, "t" : 51}

DIFFERENTIAL_EQUATION_CONDITIONS = {"a" : 1.2, "b" : 0.6, "c": 0.8, "d" : 0.3, "x0" : 2, "y0" : 1, "t0" : 0, "t" : 30}