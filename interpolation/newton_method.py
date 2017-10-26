import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.misc

def interpolant(x, nodes, function_value, divided_difference):
    sum = function_value[0]
    for i in range(len(nodes) - 1):
        item = 1.0
        for j in range(i+1):
            item *= (x - nodes[j])
        item *= divided_difference[i + 1][0]

        sum += item

    return sum

def draw_graphic(nodes, function_values, function, divided_difference, lower_bound_of_interpolation, upper_bound_of_interpolation):
    x = np.arange(lower_bound_of_interpolation, upper_bound_of_interpolation, 0.0001)

    plt.xlabel('X')
    plt.ylabel('Y')

    plt.figure(1).subplots_adjust(bottom=0.05, left=0.05, right=0.95, top=0.95, hspace=0.1)
    plt.ylim(ymax=20)
    plt.ylim(ymin=0)

    plt.plot(x, interpolant(x, nodes, function_values, divided_difference), color="red", label="Interpolant function")
    plt.plot(nodes, function_values, 'ro', label="Interpolant nodes")
    plt.plot(x, function(x), color="blue", label="Initial function")

    plt.title("Newton interpolation")
    plt.legend()
    plt.grid()

    plt.show()

def run(nodes, function):
    try:

        dimension = len(nodes)
        values = [function(nodes[i]) for i in range(dimension)]

        divided_differences = []
        divided_differences.append(values)

        final_divided_difference = 0
        for i in range(dimension - 1):
            column_of_divided_difference = []
            for j in range((dimension-1) - i):
                item = (divided_differences[i][j+1] - divided_differences[i][j])/(nodes[(j+1) + i] - nodes[j])
                column_of_divided_difference.append(item)
            divided_differences.append(column_of_divided_difference)

            isReached = True
            for item in column_of_divided_difference:
                if (abs(item) > 1e-200):
                    isReached = False
            if not isReached:
                continue
            final_divided_difference = i + 1
            break

        print("Drawing graphic...")

        draw_graphic(nodes, values, function, divided_differences, nodes[0], nodes[len(nodes)-1])
    except Exception as e:
        print("Sorry, something went wrong in Newton method: " + e.args[0])