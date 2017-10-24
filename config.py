import numpy

isDebug = False

DIFFERENCIATION_FUNCTION = lambda x: numpy.cos(x**4 + 5*x**3 + 3*x + 11)
DIFFERENCIATION_POINT = (numpy.sqrt(3) - 9) / 5

INTEGRATION_FUNCTION = lambda x: numpy.sin(1 / (x**2 + 1))
INTEGRATION_BOUNDARIES = {"a": 0.5,
                          "b": 0.8}

IMPROPER_INTEGRATION_FUNCTION = lambda x: x/((numpy.abs(1 - x**2)) ** (0.2))
IMPROPER_INTEGRATION_BOUNDARIES = {"a": 0.3,
                                   "b": 1}

DEFAULT_EPSILON = 1e-12