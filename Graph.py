class Target:
    def __init__(self, i, j):
        self.i = i
        self.j = j

    def __repr__(self):
        return "(" + str(self.i) + "," + str(self.j) + ")"

    def __hash__(self):
        return (self.i, self.j).__hash__()

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __sub__(self, other):
        return ((self.i - other.i) ** 2 + (self.j - other.j) ** 2) ** (1 / 2)


class Graph:

    def __init__(self, filename: str, n: int = 10, m: int = -1):
        self.name = filename.removesuffix(".dat")
        print("Parsing", filename)
        _ = open("data/" + filename, "r")
        if m == -1:
            print("Assuming m = n")
            m = n
        data = _.read()
        _.close()
        lines = data.split(sep="\n")
        deleted_points = set()
        self.points = set()
        for line_index in range(2, len(lines) - 1):
            try:
                line = lines[line_index]
                coord = line[line.find(":") + 1:]
                x, y = coord.split(sep=",")
                x = int(x.replace(" ", "").removeprefix("("))
                y = int(y.replace(" ", "").removesuffix(")"))
                if x >= n:
                    print("n too low, assuming bigger")
                    n = x + 1
                if y >= m:
                    print("m too low, assuming bigger")
                    m = y + 1
                deleted_points.add(Target(x, y))
            except ValueError:
                pass
        if len(deleted_points) == 0:
            print("Parsing " + filename + " as format 1 failed, trying format 2")
            for line_index in range(1, len(lines) - 1):
                line = lines[line_index]
                coord = line[line.find(":") + 1:]
                _, _, _, x, y = coord.split(sep=" ")
                x = float(x)
                y = float(y.replace(";", ""))
                if x >= n:
                    print("n too low, assuming bigger")
                    n = int(x + 1)
                if y >= m:
                    print("m too low, assuming bigger")
                    m = int(y + 1)
                self.points.add(Target(x, y))
        else:
            self.points = {Target(x, y) if Target(x, y) not in deleted_points else None for x in range(n) for y in
                           range(m)}
            self.points.discard(None)
        if len(self.points) == 0:
            raise Exception("Parsing " + filename + " failed")
        self.points.add(Target(0, 0))
        self.n = n
        self.m = m

    def __repr__(self):
        return self.points.__repr__()

    def plot(self):
        import matplotlib.pyplot as plt
        plt.plot([point.i for point in self.points], [point.j for point in self.points], 'r.')
        if not plt.gca().yaxis.get_inverted():
            plt.gca().invert_yaxis()
        plt.show()


class instance_graph:
    def __init__(self, graph: Graph, k, Rcapt, Rcom):
        self.points = graph.points.copy()
        self.points_list = list(self.points)
        self.n = graph.n
        self.m = graph.m
        self.k = k
        self.Rcom = Rcom
        self.Rcapt = Rcapt
        self.name = graph.name
        assert Rcom >= Rcapt

        self.neighborhood_dict = {point: set() for point in self.points}
        visited = set()
        for target in self.points:
            visited.add(target)
            for other in self.points - visited:
                d = target - other
                if d < Rcom:
                    self.neighborhood_dict[target].add(other)
                    self.neighborhood_dict[other].add(target)

    def __repr__(self):
        return self.points.__repr__()

    def plot(self):
        import matplotlib.pyplot as plt
        plt.plot([point.i for point in self.points], [point.j for point in self.points], 'r.')
        if not plt.gca().yaxis.get_inverted():
            plt.gca().invert_yaxis()
        plt.show()

    def evaluate(self, sensor_encoding, verbose=False):
        import matplotlib.pyplot as plt

        sensor_dict = {target: target in sensor_encoding for target in self.points}

        sensors = {point if sensor_dict[point] else None for point in self.points}
        sensors.discard(None)
        k_dict = {target: (1 if target in sensors else 0) for target in self.points}
        k_dict[Target(0, 0)] = self.k
        visited = set()
        connexes = {(point,) for point in sensors}.union({(Target(0, 0),)})
        for target in self.points:
            visited.add(target)
            for other in (self.points - visited).intersection(self.neighborhood_dict[target]):
                d = target - other
                target_connexe = None
                other_connexe = None
                if target in sensors.union({Target(0, 0)}) and other in sensors.union({Target(0, 0)}):
                    for connexe in connexes:
                        if target in connexe:
                            target_connexe = connexe
                        if other in connexe:
                            other_connexe = connexe
                    if target_connexe != other_connexe:
                        connexes.remove(target_connexe)
                        connexes.remove(other_connexe)
                        connexes.add(target_connexe + other_connexe)
                        if verbose:
                            plt.plot([target.i, other.i], [target.j, other.j], 'g')
                if d < self.Rcapt:
                    if other in sensors:
                        k_dict[target] += 1
                    if target in sensors:
                        k_dict[other] += 1

        if verbose:
            sensors = {point if sensor_dict[point] else None for point in self.points}
            sensors.discard(None)
            plt.plot([point.i for point in self.points], [point.j for point in self.points], 'r.')
            plt.plot([point.i for point in sensors], [point.j for point in sensors], 'bo')
            plt.plot([0], [0], 'bd')
            for point in sensors:
                plt.gca().add_patch(plt.Circle((point.i, point.j), self.Rcapt, color='blue', fill=False))
            if not plt.gca().yaxis.get_inverted():
                plt.gca().invert_yaxis()
            plt.savefig(self.name + ".jpg")
            plt.clf()

        if len(connexes) != 1:
            return False, len(sensors)
        for target in self.points:
            if k_dict[target] < self.k:
                return False, len(sensors)
        return True, len(sensors)

    def simulated_annealing(self):
        from math import exp
        from random import random
        n = len(self.points)
        energy_function = lambda feasible, cost: n * (1 - feasible) + cost
        current_temperature = n//2
        temperature_decrease_function = lambda T: T * 0.95 + 0.05/1000
        stopping_criteria = lambda T: T < 0.01
        iterations_per_temperature = n//2

        def next_solution(solution):
            from random import choice
            point = choice(self.points_list)
            if point in solution:
                return solution - {point}
            else:
                return solution.union({point})

        current_solution = set()
        best_solution = current_solution
        current_energy = energy_function(*self.evaluate(current_solution))
        best_score = current_energy
        print("Applying simulated annealing")

        while not stopping_criteria(current_temperature):
            for _ in range(iterations_per_temperature):
                new_solution = next_solution(current_solution)
                new_energy = energy_function(*self.evaluate(new_solution))
                energy_delta = new_energy - current_energy
                if energy_delta < 0 or random() < exp(-energy_delta / current_temperature):
                    current_solution = new_solution
                    current_energy = energy_function(*self.evaluate(current_solution))
                    if current_energy < best_score:
                        best_score = current_energy
                        best_solution = current_solution
                        print("New best, temperature:", current_temperature, " energy:", current_energy)
            print("Temperature:", current_temperature, " energy:", current_energy)
            current_temperature = temperature_decrease_function(current_temperature)
        return best_solution, best_score


if __name__ == "__main__":
    k_ = 1
    Rcapt_ = 2
    Rcom_ = 3

    def full_test(filename):
        g = instance_graph(Graph(filename), k_, Rcapt_, Rcom_)
        sol, _ = g.simulated_annealing()
        print(sol)
        print(g.evaluate(sol, True))

    filenames = ["test.dat",
                 "grille1010_1.dat",
                 "grille1010_2.dat",
                 "captANOR150_7_4.dat"
                 ]
    for file in filenames:
        full_test(file)

