from config import DIFFERENTIAL_EQUATION_CONDITIONS, DEFAULT_EPSILON, isDebug
import numpy
import matplotlib.pyplot as plt

class Differential_Equations_System:
    f1 = None
    f2 = None

    x = None
    y = None

    x0 = None
    y0 = None
    start = None
    end = None

    epsilon = None

    def __init__(self):
        initial_coditions = DIFFERENTIAL_EQUATION_CONDITIONS
        a = initial_coditions["a"]
        b = initial_coditions["b"]
        c = initial_coditions["c"]
        d = initial_coditions["d"]

        self.x0 = initial_coditions["x0"]
        self.y0 = initial_coditions["y0"]

        self.start = initial_coditions["start"]
        self.end = initial_coditions["end"]

        self.f1 = lambda t, y, x: a*x - b*x*y
        self.f2 = lambda t, y, x: -c*y + d*x*y

        self.epsilon = DEFAULT_EPSILON

    def __calculate_for_step__(self, h :float):

        t = numpy.arange(self.start, self.end + h, h)

        x = [self.x0]
        y = [self.y0]

        for i in range(len(t) - 1):
            k1 = h * self.f1(t=t[i], y=y[i], x=x[i])
            q1 = h * self.f2(t=t[i], y=y[i], x=x[i])

            k2 = h * self.f1(t=t[i] + h/2, y=y[i] + k1/2, x=x[i] + q1/2)
            q2 = h * self.f2(t=t[i] + h/2, y=y[i] + k1/2, x=x[i] + q1/2)

            k3 = h * self.f1(t=t[i] + h/2, y=y[i] + k2/2, x=x[i] + q2/2)
            q3 = h * self.f2(t=t[i] + h/2, y=y[i] + k2/2, x=x[i] + q2/2)

            k4 = h * self.f1(t=t[i] + h, y=y[i] + k3, x=x[i] + q3)
            q4 = h * self.f2(t=t[i] + h, y=y[i] + k3, x=x[i] + q3)

            x.append(
                x[i] + (k1 + 2*k2 + 2*k3 + k4)/6
            )
            y.append(
                y[i] + (q1 + 2*q2 + 2*q3 + q4)/6
            )

        return t, x, y

    def __count_error__(self, h, t, x, y):
        try:
            t_test, x_test, y_test = self.__calculate_for_step__(h / 2)

            x_difference = []
            y_difference = []

            for i in range(len(t)):
                if (len(t_test) > 2*i) and (t[i] == t_test[2*i]):
                    x_difference.append(abs(x[i] - x_test[2*i]))
                    y_difference.append(abs(y[i] - y_test[2*i]))

            Rx = max(x_difference)
            Ry = max(y_difference)
            R = max([Rx, Ry])

            return R, Rx, Ry

        except Exception as e:
            if isDebug:
                print("Exception occured at count_error function: " + str(e.args))

    def draw_graphic(self):
        try:
            h = self.epsilon ** 0.25
            print("SYSTEM OF DIFFERENTIAL EQUATIONS")

            i = 0
            while True:
                t, x, y = self.__calculate_for_step__(h)

                R, Rx, Ry = self.__count_error__(h, t, x, y)

                print("[Iteration {}]\n"
                      "Estimated error: Rx = {}\n"
                      "                 Ry = {}\n"
                      "                 R = {}\n"
                      "                 epsilon = {}\n"
                      "                 h = {}".format(i, Rx, Ry, R, self.epsilon, h))

                if (R > self.epsilon):
                    h /= 2
                    i += 1
                else:
                    break

            plt.figure(1).subplots_adjust(bottom=0.05, left=0.05, right=0.95, top=0.95, hspace=0.2)

            plt.subplot(211)
            plt.xlabel('t')

            plt.title("Time series")

            plt.plot(t, x, color="blue", label="x(t) - prey")
            plt.plot(t, y, color="red", label="y(t) - predator")

            plt.legend()
            plt.grid()

            plt.subplot(212)
            plt.xlabel("Prey")
            plt.ylabel("Predator")
            plt.title("Prey-Predator dependency")

            plt.plot(x, y, color="red")
            plt.legend()
            plt.grid()

            plt.show()
        except Exception as e:
            if isDebug:
                print("Exception occures in draw_graphic method: " + str(e.args))