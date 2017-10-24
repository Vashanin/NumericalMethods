from math_analysis import Derivation, Integration, Improper_Integral
from config import isDebug, DIFFERENCIATION_POINT
from differential_equation_system import System_Of_Differential_Equations

class Manager:
    __instance = None

    __diff__ = None
    __integral__ = None
    __improper_integral__ = None

    @staticmethod
    def inst():
        if Manager.__instance == None:
            Manager.__instance = Manager()
        return Manager.__instance

    def run(self):
        while True:
            commands = input("> ")
            if (commands == "quit"):
                break

            command = commands.split(" ")
            try:
                if (command[0] == "diff"):
                    if (command[1] == "init"):
                        if (len(command) == 2):
                            command.append(str(4))
                        self.__diff__ = Derivation(accuracy_order=int(command[2]))

                    elif (command[1] == "count"):
                        if (self.__diff__ == None):
                            print("Firstly, you have to init diff object.")
                        if (len(command) == 2):
                            command.append(str(DIFFERENCIATION_POINT))
                        self.__diff__.count(float(command[2]), disp=True)

                    elif (command[1] == "set-epsilon"):
                        if (self.__diff__ == None):
                            print("Firstly, you have to init diff object.")
                        self.__diff__.set_epsilon(float(command[2]))

                    elif (command[1] == "draw-graphic"):
                        if (self.__diff__ == None):
                            print("Firstly, you have to init diff object.")
                        self.__diff__.draw(float(command[2]), float(command[3]))
                    else:
                        print("Invalid command: " + str(command[1]))
                        print("Here are next arguments for integrate:\n"
                              "\tinit - initiation(first of all)\n"
                              "\tcount [value] - count differencial(if value isn't mentioned, program use default value)\n"
                              "\tset-epsilon [value] - set computation mistake(change default value)\n"
                              "\tdraw-graphic [lower] [upper] - draw graphic of first derivation, within mentioned boundaries\n")
                elif (command[0] == "integral"):
                    if (command[1] == "init"):
                        self.__integral__ = Integration()

                    elif (command[1] == "count"):
                        if (len(command) < 3):
                            command.append("no-value")

                        if (command[2] == "simpson-method"):
                            self.__integral__.count(method="simpson")

                        elif (command[2] == "trapezoid-method"):
                            self.__integral__.count(method="trapezoid")

                        elif (command[2] == "gaussian-method"):
                            self.__integral__.count(method="gaussian")
                        else:
                            print("You have to mention numerical method. Here are next integration methods:\n"
                                  "\tsimpson-method\n"
                                  "\ttrapezoid-method\n"
                                  "\tgaussian-method")

                    elif (command[1] == "set-epsilon"):
                        self.__integral__.set_epsilon(float(command[2]))

                    elif (command[1] == "set-boundaries"):
                        self.__integral__.set_boundaries(
                            {"a" : float(command[2]), "b" : float(command[3])}
                        )

                    else:
                        print("Invalid command: " + str(command[1]))
                        print("Here are next arguments for integrate:\n"
                              "\tinit - initiation(first of all)\n"
                              "\tcount - count integral\n"
                              "\tset-epsilon [value] - set computation mistake(change default value)\n"
                              "\tset-boundaries [lower] [upper] - set boundaries of integration(change default value)\n")
                elif (command[0] == "improper-integral"):

                    if (command[1] == "init"):
                        self.__improper_integral__ = Improper_Integral()

                    elif (command[1] == "count"):
                        if (len(command) < 3):
                            command.append("no-value")

                        if (command[2] == "cutting-method"):
                            self.__improper_integral__.count(method="cutting")

                        elif (command[2] == "quadrature-method"):
                            self.__improper_integral__.count(method="quadrature-formulas")

                        else:
                            print("You have to mention numerical method. Here are next integration methods:\n"
                                  "\tсutting-method\n"
                                  "\tquadrature-method\n")

                    elif (command[1] == "set-epsilon"):
                        self.__improper_integral__.set_epsilon(float(command[2]))

                    elif (command[1] == "set-boundaries"):
                        self.__improper_integral__.set_boundaries(
                            {"a": float(command[2]), "b": float(command[3])}
                        )
                    else:
                        print("Invalid command: " + str(command[1]))
                        print("Here are next arguments for integrate:\n"
                              "\tinit - initiation(first of all)\n"
                              "\tcount - count integral\n"
                              "\tset-epsilon [value] - set computation mistake(change default value)\n"
                              "\tset-boundaries [lower] [upper] - set boundaries of integration(change default value)\n")
                elif (command[0] == "help"):
                    print('List of сommands:\n'
                          '\tdiff init [accuracy_order] - create Differentiation instance\n'
                          '\tdiff count [point] - count differencial at point with mentioned value(change default value)\n'
                          '\tdiff set-epsilon [value] - set computation error(change default value)\n'
                          '\tdiff draw-graphic [from] [to] - draw graphic in mentioned interval\n\n'
                          '\tintegral init - create Integration instance\n'
                          '\tintegral count [method] - count integral, using mentioned numerical method\n'
                          '\tintegral set-epsilon [value] - set computation error(change default value)\n'
                          '\tintegral set-boundaries [lower] [upper] - set boundaries of integration(change default value)\n\n'
                          '\timproper-integral init - create Improper_Integral instance\n'
                          '\timproper-integral count [method] - count integral, using mentioned numerical method\n'
                          '\timproper-integral set-epsilon [value] - set computation error(change default value)\n'
                          '\timproper-integral set-boundaries [lower] [upper] - set boundaries of integration(change default value)\n'
                          '\nWARNING: Before counting differential or integral, you have to init objects!')
                else:
                    print("Invalid command: " + str(command[0]))

            except Exception as e:
                print("Something went wrong. Use a help.")
                if isDebug:
                    print("Exception: " + str(e.args))

if isDebug:
    system_of_differential_equations = System_Of_Differential_Equations()
    system_of_differential_equations.draw_graphic()
else:
    manager = Manager.inst()
    manager.run()