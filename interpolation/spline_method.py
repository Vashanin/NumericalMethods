import numpy as np
import matplotlib.pyplot as plt

S1 = -1
S2 = 1

def get_m(h, d):
    matrix = []
    coefficients = []
    dimension = len(h) + 1

    row = [0 for k in range(dimension)]
    row[0] = 1
    row[1] = 0.5
    matrix.append(row)
    coefficients.append((3/h[0])*(d[0] - S1))

    for i in range(1, dimension-1):
        row = [0 for k in range(dimension)]

        row[i-1] = h[i-1]
        row[i] = 2*(h[i-1] + h[i])
        row[i+1] = h[i]

        matrix.append(row)
        coefficients.append(6*(d[i]-d[i-1]))

    row = [0 for k in range(dimension)]
    N = len(h)
    row[N - 1] = 0.5
    row[N] = 1
    matrix.append(row)
    coefficients.append((3/h[N-1])*(S2 - d[N-1]))

    return list(np.linalg.solve(matrix, coefficients))

def linear_interpolant(x, nodes, function_value):
    y = []
    for i in range(len(x)):
        for j in range(len(nodes) - 1):
            if (nodes[j] <= x[i] and x[i] <= nodes[j+1]):
                k = (function_value[j+1] - function_value[j]) / (nodes[j+1] - nodes[j])
                b = (function_value[j+1] * nodes[j] - function_value[j] * nodes[j+1]) / (nodes[j] - nodes[j+1])

                y.append(k * x[i] + b)
                break
    return y

def cubic_interpolant(x, nodes, function_value):
    y = []

    h = [(nodes[k+1] - nodes[k]) for k in range(len(nodes) - 1)]
    d = [(function_value[k+1] - function_value[k]) / h[k] for k in range(len(nodes) - 1)]
    m = get_m(h, d)

    x = list(x)

    b0 = [function_value[k] for k in range(len(nodes) - 1)]
    b1 = [(d[k] - h[k] * (2 * m[k] + m[k+1]) / 6) for k in range(len(nodes) - 1)]
    b2 = [(m[k]/2) for k in range(len(nodes) - 1)]
    b3 = [((m[k + 1] - m[k]) / (6 * h[k])) for k in range(len(nodes) - 1)]

    for i in range(len(x)):
        for j in range(len(nodes) - 1):
            if (nodes[j] <= x[i] and x[i] <= nodes[j + 1]):
                y.append(b0[j] + b1[j]*(x[i] - nodes[j]) + b2[j]*(x[i] - nodes[j])*(x[i] - nodes[j]) + b3[j]*(x[i] - nodes[j])*(x[i] - nodes[j])*(x[i] - nodes[j]))
                break

    return np.asarray(y)

def draw_graphic(nodes, function):
    dimension = len(nodes)
    function_values = [function(nodes[i]) for i in range(dimension)]

    x = np.arange(nodes[0], nodes[len(nodes) - 1], 0.01)

    y_linear = linear_interpolant(x, nodes, function_values)
    y_cubic = cubic_interpolant(x, nodes, function_values)

    plt.xlabel('X')
    plt.ylabel('Y')

    plt.figure(1).subplots_adjust(bottom=0.03, left=0.03, right=0.97, top=0.96, hspace=0.1)
    plt.subplot(211)
    plt.plot(x, y_linear, color="red", label="Linear spline")
    plt.plot(nodes, function_values, 'ro', label="Interpolant nodes")
    plt.plot(x, function(x), color="blue", label="Initial function")
    plt.title("Spline interpolation")
    plt.legend()
    plt.grid()

    plt.subplot(212)
    plt.plot(x, y_cubic, color="red", label="Cubic spline")
    plt.plot(nodes, function_values, 'ro', label="Interpolant nodes")
    plt.plot(x, np.log(3.0 * x - 12.0), color="blue", label="Initial function")
    plt.legend()
    plt.grid()

    plt.show()

def run(nodes, function):
    try:
        draw_graphic(nodes, function)
    except Exception as e:
        print("Sorry, something went wrong in Spline method: " + e.args[0])