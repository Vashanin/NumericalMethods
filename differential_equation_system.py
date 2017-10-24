from config import DIFFERENTIAL_EQUATION_CONDITIONS, DEFAULT_EPSILON
import numpy
import matplotlib.pyplot as plt

class System_Of_Differential_Equations:
    f1 = None
    f2 = None

    k1 = None
    k2 = None
    k3 = None
    k4 = None

    q1 = None
    q2 = None
    q3 = None
    q4 = None

    x = None
    y = None

    x0 = None
    y0 = None
    t_start = None
    t_end = None

    def __init__(self):
        initial_coditions = DIFFERENTIAL_EQUATION_CONDITIONS
        a = initial_coditions["a"]
        b = initial_coditions["b"]
        c = initial_coditions["c"]
        d = initial_coditions["d"]

        self.x0 = initial_coditions["x0"]
        self.y0 = initial_coditions["y0"]

        self.t_start = initial_coditions["t0"]
        self.t_end = initial_coditions["t"]

        self.f1 = lambda y, x: a*x - b*x*y
        self.f2 = lambda y, x: -c*y + d*x*y

        self.k1 = lambda h, y, x: h * self.f1(y, x)
        self.k2 = lambda h, y, x, k1, q1: h * self.f1(y + k1/2, x + q1/2)
        self.k3 = lambda h, y, x, k2, q2: h * self.f1(y + k2/2, x + q2/2)
        self.k4 = lambda h, y, x, k3, q3: h * self.f1(y + k3, x + q3)

        self.q1 = lambda h, y, x: h * self.f2(y, x)
        self.q2 = lambda h, y, x, k1, q1: h * self.f2(y + k1/2, x + q1/2)
        self.q3 = lambda h, y, x, k2, q2: h * self.f2(y + k2/2, x + q2/2)
        self.q4 = lambda h, y, x, k3, q3: h * self.f2(y + k3, x + q3)

        self.x = lambda x, q1, q2, q3, q4: x + (q1 + 2 * q2 + 2 * q3 + q4)/6
        self.y = lambda y, k1, k2, k3, k4: y + (k1 + 2 * k2 + 2 * k3 + k4)/6

    def __calculate_for_step__(self, h):

        t = numpy.arange(self.t_start, self.t_end + h, h)

        x = [self.x0]
        y = [self.y0]

        for i in range(len(t) - 1):
            k1 = self.k1(h=h, y=y[i], x=x[i])
            q1 = self.q1(h=h, y=y[i], x=x[i])

            k2 = self.k2(h=h, y=y[i], x=x[i], k1=k1, q1=q1)
            q2 = self.q2(h=h, y=y[i], x=x[i], k1=k1, q1=q1)

            k3 = self.k3(h=h, y=y[i], x=x[i], k2=k2, q2=q2)
            q3 = self.q3(h=h, y=y[i], x=x[i], k2=k2, q2=q2)

            k4 = self.k4(h=h, y=y[i], x=x[i], k3=k3, q3=q3)
            q4 = self.q4(h=h, y=y[i], x=x[i], k3=k3, q3=q3)

            x.append(
                self.x(x=x[i], q1=q1, q2=q2, q3=q3, q4=q4)
            )
            y.append(
                self.y(y=y[i], k1=k1, k2=k2, k3=k3, k4=k4)
            )

        return {"t": t, "x": x, "y": y}

    def draw_graphic(self):
        data = self.__calculate_for_step__(0.0001)

        t = data["t"]
        print(len(t))
        x = data["x"]
        print(len(x))
        y = data["y"]
        print(len(y))

        plt.figure(1).subplots_adjust(bottom=0.05, left=0.05, right=0.95, top=0.95, hspace=0.2)

        plt.subplot(211)
        plt.xlabel('t')

        plt.title("Time series")

        plt.ylim(-10, 10)
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