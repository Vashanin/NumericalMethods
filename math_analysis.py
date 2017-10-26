import numpy
import matplotlib.pyplot as plt
import scipy.optimize
import scipy.misc
from config import isDebug, DIFFERENCIATION_FUNCTION, INTEGRATION_FUNCTION, INTEGRATION_BOUNDARIES,\
    IMPROPER_INTEGRATION_BOUNDARIES, IMPROPER_INTEGRATION_FUNCTION, DEFAULT_EPSILON
from interpolation.lagrange_method import Lagrange_Interpolation

class Derivation:
    def __init__(self, accuracy_order=4):
        try:
            self.initial_function = DIFFERENCIATION_FUNCTION
            self.__epsilon__ = DEFAULT_EPSILON

            if (accuracy_order == 2 or accuracy_order == 4):
                self.accuracy_order = accuracy_order
            else:
                raise Exception("Invalid accuracy order!")

        except Exception as e:
            print("Differentiation constructor error: " + str(e.args))

    @staticmethod
    def maximize_the_derivative(function, order_of_the_derivation, lower, upper):
        derivative = lambda x: (-1) * scipy.misc.derivative(function, x, n=order_of_the_derivation,
                                                            order=order_of_the_derivation + 1)

        x0 = numpy.asarray([0.1])
        res = scipy.optimize.minimize(derivative, x0, method='TNC', bounds=((lower, upper),),
                                      options={'xtol': 1e-10, 'disp': False})
        if isDebug:
            print(str(res) + "\n")

        return (-1) * res.fun[0]

    def get_max_valuation(self):
        return 1e12
    def get_step(self):
        M = self.get_max_valuation()
        e = self.get_epsilon()
        return ((25*e)/(4*M))**(1/5)
    def set_epsilon(self, epsilon: float):
        self.__epsilon__ = epsilon
    def get_epsilon(self):
        return self.__epsilon__
    def set_initial_function(self, initial_function):
        self.initial_function = initial_function
    def set_accuracy_order(self, accuracy_order):
        if (accuracy_order == 2 or accuracy_order == 4):
            self.accuracy_order = accuracy_order

    def count(self, initial_point, disp=False):
        result = 0.0

        if disp:
            print("DIFFERENTIATION")

        def D(h):
            f = self.initial_function
            x = initial_point

            return (f(x+h) - f(x-h)) / (2*h)

        h = self.get_step()

        if (self.accuracy_order == 2):
            result = D(h)

        if (self.accuracy_order == 4):
            result = (4*D(h) - D(2*h)) / 3

        if isDebug:
            print("\th: {}\n"
                  "\tD(h): {}\n"
                  "\tD(2h): {}\n".format(h, D(h), D(2*h)))

        error_estimation = (3 * self.__epsilon__) / (2 * h) + (self.get_max_valuation() * h ** 4) / 30

        if disp:
            print("\tPoint: {}\n"
                  "\tDifferential value: {}\n"
                  "\tError estimation: {}\n".format(initial_point, result, error_estimation))

        return result
    def draw(self, lower, upper):
        atoms = numpy.arange(lower, upper, 0.001)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.figure(1).subplots_adjust(bottom=0.03, left=0.03, right=0.97, top=0.96, hspace=0.1)

        plt.subplot(211)
        plt.title("Differentiation")

        plt.plot(atoms, self.initial_function(atoms), color="blue", label="f(x)")
        plt.legend()
        plt.grid()

        plt.subplot(212)

        plt.plot(atoms, self.count(atoms), color="red", label="df(x)/dx")
        plt.legend()
        plt.grid()

        plt.show()

class Integration:
    def __init__(self, initial_function=INTEGRATION_FUNCTION, boundaries=INTEGRATION_BOUNDARIES):
        self.initial_function = initial_function
        self.boundaries = boundaries
        self.epsilon = DEFAULT_EPSILON

    def set_epsilon(self, epsilon):
        self.epsilon = epsilon
    def set_boundaries(self, boundaries):
        self.boundaries = boundaries
    def set_initial_function(self, initial_function):
        self.initial_function = initial_function

    def count(self, method, disp=True, show_error_estimation=True):
        result = None

        if (method == "trapezoid"):
            result = self.trapezoid_method(disp=disp, show_error_estimetion=show_error_estimation)
        elif (method == "simpson"):
            result = self.simpson_method(disp=disp, show_error_estimetion=show_error_estimation)
        elif (method == "gaussian"):
            result = self.gaussian_method(disp=disp, show_error_estimetion=show_error_estimation)
        else:
            print("Undefined numerical method: {}".format(method))

        return result

    def trapezoid_method(self, disp, show_error_estimetion):
        try:
            if disp:
                print("TRAPEZOID METHOD")

            def get_step(M2):
                a = self.boundaries["a"]
                b = self.boundaries["b"]

                max_h = numpy.sqrt(12 * self.epsilon /
                                   (M2 * numpy.abs(b - a)))
                min_n = numpy.abs(b-a) / max_h

                INCREASE_ACCURACY = 3
                n = numpy.round(min_n) + INCREASE_ACCURACY
                h = numpy.abs(b-a) / n

                return {"amount": int(n), "step": float(h)}

            a = self.boundaries["a"]
            b = self.boundaries["b"]

            M2 = Derivation.maximize_the_derivative(self.initial_function,
                                                    order_of_the_derivation=4,
                                                    lower=a,
                                                    upper=b)
            M2 = abs(M2)
            optimal = get_step(M2)
            n = optimal["amount"]
            h = optimal["step"]

            f = self.initial_function
            x = list(numpy.arange(a, b+h, h))

            integral = f(x[0]) + f(x[n])

            for i in range(1, n):
                integral += 2 * f(x[i])

            integral *= h/2

            error_estimation = None

            if show_error_estimetion:
                M2 = Derivation.maximize_the_derivative(self.initial_function,
                                                        order_of_the_derivation=2,
                                                        lower=a,
                                                        upper=b)
                M2 = abs(M2)
                error_estimation = M2 * numpy.abs(b - a) * (h ** 2) / 12

            if isDebug:
                print("\tM2: {}\n"
                      "\tn: {}\n"
                      "\th: {}\n"
                      "\tx: {}\n".format(str(M2), str(n), str(h), str(numpy.asarray(x))))

            if disp:
                print("\tIntegral value: {}\n"
                      "\tError estimation: {}\n".format(str(integral), str(error_estimation)))

            return {"integral_value": integral,
                    "error_estimation": error_estimation}
        except Exception as e:
            if isDebug:
                print("Exception occurs in Integration().__trapezoid_method__: " + str(e.args))
    def simpson_method(self, disp, show_error_estimetion):
        try:
            if disp:
                print("SIMPSON METHOD")

            def get_step(amount_of_segments):
                return (self.boundaries["b"] - self.boundaries["a"]) / (2 * amount_of_segments)

            f = self.initial_function

            m = 0
            iteration = []
            while True:
                m += 1
                h = get_step(m)

                x = list(numpy.arange(self.boundaries["a"], self.boundaries["b"] + h, h))

                integral = 0.0
                integral += f(x[0])

                for i in range(1, m+1):
                    integral += 4*f(x[2*i-1])

                for i in range(1, m):
                    integral += 2*f(x[2*i])

                integral += f(x[2 * m])
                integral /= (3/h)

                if isDebug:
                    print("\t[Iteration #{}]\n"
                          "\tm: {}\n"
                          "\tIntegral value: {}\n".format(len(iteration), m, integral))

                iteration.append(integral)
                iteration_amount = len(iteration)
                if (iteration_amount > 1):
                    if (numpy.abs(iteration[iteration_amount - 1] - iteration[iteration_amount-2]) < self.epsilon):
                        break

            integral = iteration[len(iteration)-1]
            a = self.boundaries["a"]
            b = self.boundaries["b"]

            error_estimation = None
            if show_error_estimetion:
                M4 = Derivation.maximize_the_derivative(self.initial_function,
                                                        order_of_the_derivation=4,
                                                        lower=a,
                                                        upper=b)
                error_estimation = M4 * numpy.abs(b - a) * (get_step(m) ** 4) / 180

            if isDebug:
                print("\tM4: {}\n"
                      "\th: {}\n"
                      "\tm: {}".format(str(M4), str(get_step(m)), str(str(m))))

            if disp:
                print("\tAmount of iterations: {}\n"
                      "\tIntegral value: {}\n"
                      "\tError estimation: {}\n".format(str(len(iteration)), str(integral), str(error_estimation)))

            return {"integral_value": integral,
                    "error_estimation": error_estimation,
                    "amount_of_iterations": len(iteration)}

        except Exception as e:
            if isDebug:
                print("Exception occured in Integration().__simpson_method__: " + str(e.args))
    def gaussian_method(self, disp, show_error_estimetion):
        try:
            if disp:
                print("GAUSSIAN METHOD")

            a = self.boundaries["a"]
            b = self.boundaries["b"]

            f = lambda t: self.initial_function((a + b) / 2 + t * (b - a) / 2)

            Xi = [-0.90617985, -0.53846931, 0, 0.53846931, 0.90617985]
            q = [0.23692688, 0.47862868, 0.56888889, 0.47862868, 0.23692688]
            f_Xi = [f(Xi[i]) for i in range(len(Xi))]

            integral = ((b-a)/2) * (q[0]*f(Xi[0])
                            + q[1]*f(Xi[1])
                            + q[2]*f(Xi[2])
                            + q[3]*f(Xi[3])
                            + q[4]*f(Xi[4]))
            if isDebug:
                print("\tXi: {}\n"
                      "\tf(Xi): {}\n"
                      "\tq: {}\n".format(str(Xi), str(f_Xi), str(q)))

            if disp:
                print("\tIntegral value: {}\n".format(integral))

            return {"integral_value": integral}
        except Exception as e:
            print("Exception occured in Integration().__gaussian_method__: " + str(e.args))

class Improper_Integral(Integration):
    def __init__(self, initial_function=IMPROPER_INTEGRATION_FUNCTION, boundaries=IMPROPER_INTEGRATION_BOUNDARIES):
        Integration.__init__(self, initial_function, boundaries)

        self.initial_function = initial_function
        self.boundaries = boundaries

    def count(self, method, disp=True, show_error_estimation=True):
        result = None

        if method == "cutting":
            result = self.cutting_method(disp=disp)
        elif method == "quadrature-formulas":
            result = self.quadrature_formulas_method(disp=disp)
        elif method == "monte-carlo":
            result = self.monte_carlo_method(disp=disp)
        else:
            print("Undefined numerical method: {}".format(method))

        return result

    def cutting_method(self, disp):
        try:
            if disp:
                print("IMPROPER INTEGRAL")

            a = self.boundaries["a"]
            b = self.boundaries["b"]

            integral = Integration()
            integral.set_initial_function(self.initial_function)
            iterations = []

            delta = 1e-15

            while True:
                integral.set_boundaries({"a" : a,
                                         "b": b - delta})

                integral.set_epsilon(epsilon=1e-10)

                value = integral.count(method="trapezoid", disp=isDebug, show_error_estimation=False)

                iterations.append(value["integral_value"])

                if (len(iterations) > 1):
                    if (numpy.abs(iterations[len(iterations) - 1] - iterations[len(iterations) - 2]) < (self.epsilon)**2):
                        break

                delta -= delta/10

            result = iterations[len(iterations) - 1]
            if disp:
                print("\tImproper integral value: {}\n".format(result))

            return {"integral_value": result}
        except Exception as e:
            if isDebug:
                print("Exception occured in Improper_Integral().__cutting_method__: " + str(e.args))
    def quadrature_formulas_method(self, disp):
        try:
            if disp:
                print("IMPROPER INTEGRAL")

            a = self.boundaries["a"]
            b = self.boundaries["b"]

            nodes_amount = 10

            optimal_nodes = [((b-a)/2) * numpy.cos((2*k + 1) * numpy.pi / (2 * nodes_amount)) + ((b+a)/2) for k in range(nodes_amount)]

            lagrange_interpolation = Lagrange_Interpolation(function=self.initial_function,
                                                            nodes=optimal_nodes,
                                                            boundaries_of_interpolation=self.boundaries)

            interpolant = lambda x: lagrange_interpolation.interpolant(x)

            if isDebug:
                print("\tNodes amount: {}\n"
                      "\tInterpolation nodes: {}\n".format(nodes_amount, optimal_nodes))
                lagrange_interpolation.draw_graphic()

            integral = Integration(initial_function=interpolant, boundaries=self.boundaries)
            integral.set_epsilon(epsilon=1e-4)
            result = integral.count(method="simpson", disp=isDebug, show_error_estimation=False)

            if disp:
                print("\tImproper integral value: {}\n".format(result["integral_value"]))

            return result
        except Exception as e:
            if isDebug:
                print("Exception occured in Improper_Integral().__quadrature_formulas_method__: " + str(e.args))
    def monte_carlo_method(self, disp):
        try:
            if disp:
                print("IMPROPER INTEGRAL")

            a = self.boundaries["a"]
            b = self.boundaries["b"]
            f = self.initial_function

            amount_of_simulations = 1000000
            simulations = numpy.random.uniform(low=a, high=b, size=amount_of_simulations)

            integral_value = 0.0
            for x in simulations:
                integral_value += f(x) * (b - a)

            integral_value /= amount_of_simulations

            if disp:
                print("\tImproper integral value: {}\n".format(integral_value))

            return {"integral_value": integral_value}
        except Exception as e:
            if isDebug:
                print("Exception occured at monte_carlo_method: " + str(e.args))