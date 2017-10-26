import numpy
import scipy.optimize

def run(data):
    print("\nLINEAR APPROXIMATION")

    x = data["x"]
    y = data["y"]

    def function(a):
        sum = 0.0
        for i in range(len(y)):
            sum += (y[i] - (a[0]*x[i] + a[1])) ** 2
        return sum

    x0 = numpy.array([6.1, 1.2])
    res = scipy.optimize.minimize(function, x0, method='nelder-mead', options={'xtol': 1e-10, 'disp': True})

    print("Approximation coefficients: " + str(res.x))

    return (lambda x: res.x[0]*x + res.x[1])