import numpy as np
import math
import matplotlib.pyplot as plt

def interpolant(argument, coefficients):
    sum = 0.0
    for i in range(len(coefficients)):
        sum += coefficients[i] * np.power(argument, i)

    return sum

def draw_graphic(nodes, function_values, function, coefficients, lower_bound_of_interpolation, upper_bound_of_interpolation):
    x = np.arange(lower_bound_of_interpolation, upper_bound_of_interpolation, 0.001)

    plt.title('Polynomial interpolation')
    plt.xlabel('X')
    plt.ylabel('Y')

    plt.plot(x, interpolant(x, coefficients), color="red", label="Interpolant function")
    plt.plot(nodes, function_values, 'ro', label="Interpolant nodes")
    plt.plot(x, function(x), color="blue", label="Initial function")

    plt.title("Polynomial interpolation")
    plt.legend()
    plt.grid()

    plt.show()

def run(nodes, function):
    try:
        dimension = len(nodes)
        function_values = [function(nodes[i]) for i in range(dimension)]
        matrix_of_coefficients = [[math.pow(nodes[i], j) for j in range(dimension)] for i in range(dimension)]
        coefficients = np.linalg.solve(np.array(matrix_of_coefficients), np.array(function_values))

        print("Coefficients of the equation: " + str([(coefficient) for coefficient in coefficients]))

        while True:
            user_input = input("\nDo you want to check error at some point(y/n) >> ")
            if (user_input == "n"):
                break
            if (user_input == "y"):
                user_point = float(input("Enter some point(float value) >> "))
                print("| Your point: " + str(user_point))
                interpolant_value_at_user_point = interpolant(user_point, coefficients)
                print("| Interpolant function value in the point: " + str(interpolant_value_at_user_point))
                function_value_at_user_point = function(user_point)
                print("| Initial function value in the point: " + str(function_value_at_user_point))
                error = np.abs(interpolant_value_at_user_point - function_value_at_user_point)
                print("| Absolute error: " + str(error))

        print("Drawing graphic...")
        draw_graphic(nodes, function_values, function, coefficients, nodes[0], nodes[len(nodes) - 1])
    except Exception as e:
        print("Sorry, something went wrong in Polynomial method: " + e.args[0])