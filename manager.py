from math_analysis import Derivation, Integration, Improper_Integral
from config import isDebug, DIFFERENCIATION_POINT
from equations_system import Differential_Equations_System

class Manager:
    __instance = None

    is_initiate = False

    diff = None
    integral = None
    improper_integral = None
    equations_system = None

    @staticmethod
    def inst():
        if Manager.__instance == None:
            Manager.__instance = Manager()

            if not Manager.is_initiate:
                Manager.diff = Derivation(accuracy_order=4)
                Manager.integral = Integration()
                Manager.improper_integral = Improper_Integral()
                Manager.equations_system = Differential_Equations_System()
                Manager.is_initiate = True

        return Manager.__instance

    def run(self):
        while True:
            commands = input("> ")
            if (commands == "quit"):
                break

            command = commands.split(" ")
            try:
                if (command[0] == "diff"):
                    if (command[1] == "set-accuracy-order"):
                        accuracy_order = int(command[2])
                        Manager.diff.set_accuracy_order(accuracy_order)

                    elif (command[1] == "count"):
                        if (len(command) == 2):
                            command.append(str(DIFFERENCIATION_POINT))
                        Manager.diff.count(float(command[2]), disp=True)

                    elif (command[1] == "set-epsilon"):
                        Manager.diff.set_epsilon(float(command[2]))

                    elif (command[1] == "draw-graphic"):
                        Manager.diff.draw(float(command[2]), float(command[3]))
                    else:
                        print("Invalid command: " + str(command[1]))
                        print("Here are next arguments for integrate:\n"
                              "\tset-accuracy-order - set accuracy order of derivation (2 or 4)\n"
                              "\tcount [value] - count differencial(if value isn't mentioned, program use default value)\n"
                              "\tset-epsilon [value] - set computation mistake(change default value)\n"
                              "\tdraw-graphic [lower] [upper] - draw graphic of first derivation, within mentioned boundaries\n")
                elif (command[0] == "integral"):
                    if (command[1] == "count"):
                        if (len(command) < 3):
                            command.append("no-value")

                        if (command[2] == "simpson-method"):
                            Manager.integral.count(method="simpson")

                        elif (command[2] == "trapezoid-method"):
                            Manager.integral.count(method="trapezoid")

                        elif (command[2] == "gaussian-method"):
                            Manager.integral.count(method="gaussian")
                        else:
                            print("You have to mention numerical method. Here are next integration methods:\n"
                                  "\tsimpson-method\n"
                                  "\ttrapezoid-method\n"
                                  "\tgaussian-method")

                    elif (command[1] == "set-epsilon"):
                        Manager.integral.set_epsilon(float(command[2]))

                    elif (command[1] == "set-boundaries"):
                        Manager.integral.set_boundaries(
                            {"a" : float(command[2]), "b" : float(command[3])}
                        )

                    else:
                        print("Invalid command: " + str(command[1]))
                        print("Here are next arguments for integrate:\n"
                              "\tcount - count integral\n"
                              "\tset-epsilon [value] - set computation mistake(change default value)\n"
                              "\tset-boundaries [lower] [upper] - set boundaries of integration(change default value)\n")
                elif (command[0] == "improper-integral"):
                    if (command[1] == "count"):
                        if (len(command) < 3):
                            command.append("no-value")

                        if (command[2] == "cutting-method"):
                            Manager.improper_integral.count(method="cutting")
                        elif (command[2] == "quadrature-method"):
                            Manager.improper_integral.count(method="quadrature-formulas")
                        elif (command[2] == "monte-carlo-method"):
                            Manager.improper_integral.count(method="monte-carlo")
                        else:
                            print("You have to mention numerical method. Here are next integration methods:\n"
                                  "\tсutting-method\n"
                                  "\tquadrature-method\n"
                                  "\tmonte-carlo-method\n")

                    elif (command[1] == "set-epsilon"):
                        Manager.improper_integral.set_epsilon(float(command[2]))

                    elif (command[1] == "set-boundaries"):
                        Manager.improper_integral.set_boundaries(
                            {"a": float(command[2]), "b": float(command[3])}
                        )
                    else:
                        print("Invalid command: " + str(command[1]))
                        print("Here are next arguments for integrate:\n"
                              "\tinit - initiation(first of all)\n"
                              "\tcount - count integral\n"
                              "\tset-epsilon [value] - set computation mistake(change default value)\n"
                              "\tset-boundaries [lower] [upper] - set boundaries of integration(change default value)\n")
                elif (command[0] == "equations-system"):
                    if (command[1] == "count"):
                        Manager.equations_system.draw_graphic()
                    else:
                        print("Invalid command: " + str(command[1]))
                        print("Here are next arguments for integrate:\n"
                              "\tcount - compute results and show as graphic\n")
                elif (command[0] == "help"):
                    print('List of сommands:\n'
                          '\tdiff set-accuracy_order - set accuracy order (2 or 4)\n'
                          '\tdiff count [point] - count differencial at point with mentioned value(change default value)\n'
                          '\tdiff set-epsilon [value] - set computation error(change default value)\n'
                          '\tdiff draw-graphic [from] [to] - draw graphic in mentioned interval\n\n'
                          '\tintegral count [method] - count integral, using mentioned numerical method\n'
                          '\tintegral set-epsilon [value] - set computation error(change default value)\n'
                          '\tintegral set-boundaries [lower] [upper] - set boundaries of integration(change default value)\n\n'
                          '\timproper-integral count [method] - count integral, using mentioned numerical method\n'
                          '\timproper-integral set-epsilon [value] - set computation error(change default value)\n'
                          '\timproper-integral set-boundaries [lower] [upper] - set boundaries of integration(change default value)\n\n'
                          '\tequations-system count - compute result and show as graphic\n')
                else:
                    print("Invalid command: " + str(command[0]))

            except Exception as e:
                print("Something went wrong. Use a help.")
                if isDebug:
                    print("Exception: " + str(e.args))


manager = Manager.inst()
manager.run()