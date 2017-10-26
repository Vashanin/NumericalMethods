import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.misc

class Lagrange_Interpolation:
    def __init__(self, function, nodes, boundaries_of_interpolation=None):
        self.__function__ = function

        nodes.sort()
        self.nodes = nodes
        self.function_values = [function(nodes[i]) for i in range(len(nodes))]
        self.boundaries_of_interpolation = boundaries_of_interpolation \
            if (boundaries_of_interpolation != None) \
            else {"a": nodes[0], "b": nodes[len(nodes) - 1]}

    def interpolant(self, x):
        sum = 0.0
        n = len(self.function_values)
        for i in range(n):
            numerator = 1.0
            denominator = 1.0
            for j in range(n):
                if (j != i):
                    numerator *= (x - self.nodes[j])
                    denominator *= (self.nodes[i] - self.nodes[j])
            sum += self.function_values[i] * numerator / denominator

        return sum

    def draw_graphic(self):
        x = np.arange(self.boundaries_of_interpolation["a"], self.boundaries_of_interpolation["b"], 0.01)

        plt.xlabel('X')
        plt.ylabel('Y')

        plt.figure(1).subplots_adjust(bottom=0.05, left=0.05, right=0.95, top=0.95, hspace=0.1)

        plt.plot(x, self.interpolant(x), color="red", label="Interpolant function")
        plt.plot(self.nodes, self.function_values, 'ro', label="Interpolant nodes")
        plt.plot(x, self.__function__(x), color="blue", label="Initial function")

        plt.title("Lagrange interpolation")
        plt.legend()
        plt.grid()

        plt.show()


    def run(self, disp=False):
        try:
            dimension = len(self.nodes)

            while True:
                user_input = input("\nDo you want to check error at some point(y/n) >> ")
                if (user_input == "n"):
                    break
                if (user_input == "y"):
                    user_point = float(input("Enter some point(float value) >> "))
                    print("Your point: " + str(user_point))
                    interpolant_value_at_user_point = self.interpolant(user_point)
                    print("Interpolant function value in the point: " + str(interpolant_value_at_user_point))
                    function_value_at_user_point = self.__function__(user_point)
                    print("Initial function value in the point: " + str(function_value_at_user_point))

                    inc = 2
                    if (dimension % 2 == 0):
                        inc = 3

                    error = scipy.misc.derivative(self.__function__, user_point, n=(dimension+1), order=(dimension+inc)) / math.factorial(dimension + 1)
                    for i in range(dimension):
                        error *= user_point - self.nodes[i]

                    print("Error evaluation: " + str(error))
            if disp:
                print("Drawing graphic...")
                self.draw_graphic()
        except Exception as e:
            print("Sorry, something went wrong in Polynomial method: " + e.args[0])