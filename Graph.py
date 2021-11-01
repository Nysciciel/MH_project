import matplotlib.pyplot as plt
from math import sqrt, exp, ceil
from random import random, sample, expovariate, seed, shuffle
from datetime import datetime
from time import time
from os import mkdir
from collections import deque


class Target:
    def __init__(self, i, j):
        self.i = i
        self.j = j

    def __repr__(self):
        return "(" + str(self.i) + "," + str(self.j) + ")"

    def __hash__(self):
        return hash((self.i, self.j))

    def __eq__(self, other):
        return self.i == other.i and self.j == other.j

    def __sub__(self, other):
        return sqrt((self.i - other.i) ** 2 + (self.j - other.j) ** 2)


class Graph:

    def __init__(self, filename: str):
        print("Parsing", filename)
        _ = open("data/" + filename, "r")
        data = _.read()
        _.close()
        lines = data.split(sep="\n")
        deleted_points = set()
        self.points = set()
        if filename.startswith("grille"):
            dim_str = filename[:filename.find("_")].replace("grille", "")
            n = m = int(dim_str[:len(dim_str) // 2])
            self.name = filename.removesuffix(".dat")
            for line_index in range(2, len(lines) - 1):
                try:
                    line = lines[line_index]
                    coord = line[line.find(":") + 1:]
                    x, y = coord.split(sep=",")
                    x = int(x.replace(" ", "").removeprefix("("))
                    y = int(y.replace(" ", "").removesuffix(")"))
                    assert x < n
                    assert y < m
                    deleted_points.add(Target(x, y))
                except ValueError:
                    pass
            self.points = {Target(x, y) if Target(x, y) not in deleted_points else None for x in range(n) for y in
                           range(m)}
            self.points.discard(None)
        elif filename.startswith("captANOR"):
            dim_str = filename[filename.find("_") + 1:].replace(".dat", "")
            m = n = int(dim_str[:dim_str.find("_")])
            self.name = filename.removesuffix(".dat")
            if len(deleted_points) == 0:
                for line_index in range(2, len(lines) - 1):
                    line = lines[line_index]
                    coord = line[line.find(":") + 1:]
                    _, _, _, x, y = coord.split(sep=" ")
                    x = float(x)
                    y = float(y.replace(";", ""))
                    assert x < n
                    assert y < m
                    self.points.add(Target(x, y))
        else:
            raise Exception("Not a known graph format")
        if len(self.points) == 0:
            raise Exception("Parsing " + filename + " failed")
        self.points.add(Target(0, 0))
        self.n = n
        self.m = m

    def __repr__(self):
        return self.points.__repr__()

    def plot(self):
        plt.plot([point.i for point in self.points], [point.j for point in self.points], 'r.')
        if not plt.gca().yaxis.get_inverted():
            plt.gca().invert_yaxis()
        plt.title(self.name)
        plt.show()
        plt.clf()


class instance_graph:
    def __init__(self, graph: Graph, k_, Rcapt_, Rcom_):
        self.points = graph.points.copy()
        self.points_list = list(self.points)
        self.n = graph.n
        self.m = graph.m
        self.k = k_
        self.Rcom = Rcom_
        self.Rcapt = Rcapt_
        self.name = graph.name
        assert self.Rcom >= self.Rcapt

        self.neighborhood_dict = {point: set() for point in self.points}
        visited = set()
        for target in self.points:
            visited.add(target)
            for other in self.points - visited:
                d = target - other
                if d <= self.Rcom:
                    self.neighborhood_dict[target].add(other)
                    self.neighborhood_dict[other].add(target)
        self.eval_dict = {}

    def __repr__(self):
        return self.points.__repr__()

    def plot(self):
        plt.plot([point.i for point in self.points], [point.j for point in self.points], 'r.')
        if not plt.gca().yaxis.get_inverted():
            plt.gca().invert_yaxis()
        plt.show()

    def evaluate(self, sensor_encoding, verbose=False, seed_=None):
        if (tuple(sensor_encoding) in self.eval_dict) and not verbose:
            return self.eval_dict[tuple(sensor_encoding)]
        sensor_dict = {target: target in sensor_encoding for target in self.points}
        sensors = {point if sensor_dict[point] else None for point in self.points}
        sensors.discard(None)
        k_dict = {target: (1 if target in sensors else 0) for target in self.points}
        k_dict[Target(0, 0)] = self.k
        visited = set()
        connexes = {(point,) for point in sensors}.union({(Target(0, 0),)})
        for target in sensors:
            visited.add(target)
            for other in self.neighborhood_dict[target] - visited:
                d = target - other
                if other in sensors or other == Target(0, 0):
                    target_connexe = other_connexe = None
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
                if d <= self.Rcapt:
                    if other in sensors and k_dict[target] < self.k:
                        k_dict[target] += 1
                    if target in sensors and k_dict[other] < self.k:
                        k_dict[other] += 1

        if verbose:
            sensors = {point if sensor_dict[point] else None for point in self.points}
            sensors.discard(None)
            plt.plot([point.i for point in self.points], [point.j for point in self.points], 'r.')
            plt.plot([point.i for point in sensors], [point.j for point in sensors], 'bo')
            plt.plot([0], [0], 'bd')
            if seed_ is not None:
                plt.title("Seed:" + str(seed_))
            for target in self.points:
                if k_dict[target] < self.k:
                    plt.plot([target.i], [target.j], color='black', marker='.')
            for point in sensors:
                plt.gca().add_patch(plt.Circle((point.i, point.j), self.Rcapt, color='blue', fill=False))
            if not plt.gca().yaxis.get_inverted():
                plt.gca().invert_yaxis()
            try:
                plt.savefig(self.name + "_" + str(self.k) + str(self.Rcapt) + str(self.Rcom) + ".jpg")
            except OSError:
                print("Couldn't save " + self.name + ".jpg")
            plt.clf()
            plt.close()

        connected = all({k_dict[target] == self.k for target in self.points})
        degree = sum(k_dict.values()) / (self.k * len(self.points)) if not connected else 1
        connexity = 0 if len(sensors) == 1 else (len(sensors) - len(connexes)) / (len(sensors) - 1)
        feasible = len(connexes) == 1 and connected
        self.eval_dict[tuple(sensor_encoding)] = feasible, len(sensors) / len(self.points), degree, connexity
        return feasible, len(sensors) / len(self.points), degree, connexity

    def simulated_annealing(self, seed_=None):
        seed_ = abs(int(hash(str(datetime.now())))) if seed_ is None else seed_
        seed(seed_)

        n = len(self.points)
        epsilon = 0.0005
        q = 0.90
        energy_function = lambda feasible, cost, degree, connexity: \
            3 * (1 - feasible) + \
            2 * (1 - degree) + \
            2 * (1 - connexity) + \
            cost

        current_temperature = 10
        temperature_decrease_function = lambda T: T * q + (1 - q) * (epsilon / 3)
        stopping_criteria = lambda T: T < epsilon
        iterations_per_temperature = n // 3

        def next_solution(solution):
            n_ = len(self.points)
            lambd = 1
            dist = ceil(expovariate(1 / lambd))
            while dist > n_:
                dist = ceil(expovariate(1 / lambd))
            points = sample(self.points_list, dist)
            for point in points:
                if point in solution:
                    solution = solution - {point}
                else:
                    solution = solution.union({point})
            return solution

        current_solution = self.points.copy()
        best_solution = current_solution
        current_energy = energy_function(*self.evaluate(current_solution))
        best_score = float('inf')
        # print("Applying simulated annealing, k:", self.k,
        #       ",Rcapt:", self.Rcapt, ",Rcom:", self.Rcom, ",seed:", seed_)

        while not stopping_criteria(current_temperature):
            for _ in range(iterations_per_temperature):
                new_solution = next_solution(current_solution)
                new_eval = self.evaluate(new_solution)
                new_energy = energy_function(*new_eval)
                energy_delta = new_energy - current_energy
                if energy_delta < 0 or random() < exp(-energy_delta / current_temperature):
                    current_solution = new_solution
                    current_energy = new_energy
                    if current_energy < best_score and new_eval[0]:
                        best_score = current_energy
                        best_solution = current_solution
                        # print("New best, temperature:", current_temperature,
                        #       " eval:", new_eval,
                        #       " sensors:", len(best_solution))
            # print("Temperature:", current_temperature, " eval:", new_eval, "energy:", new_energy)
            current_temperature = temperature_decrease_function(current_temperature)
        old_best_solution = best_solution.copy()
        for target in old_best_solution:
            new_eval = self.evaluate(best_solution - {target})
            if new_eval[0]:
                # print("late optimization")
                best_solution -= {target}
        return best_solution

    def naive_reduction(self, seed_=None):
        seed_ = abs(int(hash(str(datetime.now())))) if seed_ is None else seed_
        seed(seed_)

        current_solution = self.points.copy()
        optimal = False
        while not optimal:
            optimal = True
            points = list(current_solution.copy())
            shuffle(points)
            for point in points:
                if self.evaluate(current_solution - {point})[0]:
                    current_solution.remove(point)
                    optimal = False
        return current_solution

    def tabou(self, seed_=None):
        seed_ = abs(int(hash(str(datetime.now())))) if seed_ is None else seed_
        seed(seed_)

        # iterations = len(self.points)*3
        duration = 600 * len(self.points) / 1500

        tabou_length = len(self.points) // 3

        cost_function = lambda feasible, cost, degree, connexity: \
            3 * (1 - feasible) + \
            2 * (1 - degree) + \
            2 * (1 - connexity) + \
            cost

        def neighborhood(solution):
            return {tuple(solution - {point}) if point in solution else
                    tuple(solution.union({point})) for point in self.points}

        current_solution = self.points.copy()
        best_solution = current_solution
        best_cost = float('inf')

        tabou = deque([], tabou_length)

        t0 = time()
        it = 0
        while time() - t0 < duration:
            it += 1
            neighbors = neighborhood(current_solution)
            for invalid in tabou:
                neighbors.discard(tuple(invalid))
            if len(neighbors) == 0:
                break
            tabou.append(current_solution)
            best_cost_neighbor = float('inf')
            for neighbor in neighbors:
                neighbor = set(neighbor)
                eval_ = self.evaluate(neighbor)
                new_cost = cost_function(*eval_)
                if new_cost < best_cost_neighbor:
                    current_solution = neighbor
                    best_cost_neighbor = new_cost
                    if new_cost < best_cost:
                        best_solution = current_solution
                        best_cost = new_cost

        old_best_solution = best_solution.copy()
        for target in old_best_solution:
            new_eval = self.evaluate(best_solution - {target})
            if new_eval[0]:
                # print("late optimization")
                best_solution -= {target}
        print("iterations:", it)
        return best_solution


def full_test(graph, seed_=None, solver=instance_graph.simulated_annealing):
    seed_ = abs(int(hash(str(datetime.now())))) if seed_ is None else seed_
    t0 = time()
    print("\nTest for", graph.name, "k=", graph.k, "Rcapt=", graph.Rcapt, "Rcom=", graph.Rcom,
          "algorithm=", solver.__name__, )
    sol = solver(graph, seed_)
    print("result:", len(sol), "time:", time() - t0)
    graph.evaluate(sol, True, seed_)


def make_stats(graph, n, solver=instance_graph.simulated_annealing):
    print("\nStats for", graph.name, "k=", graph.k, "Rcapt=", graph.Rcapt, "Rcom=", graph.Rcom,
          "Times=", n, "algorithm=", solver.__name__)
    t0 = time()
    plt.hist([len(solver(graph)) for _ in range(n)], density=True)
    plt.title("Avg time:" + str((time() - t0) / n))
    try:
        mkdir("results/" + graph.name + "_" + str(graph.k) + str(graph.Rcapt) + str(graph.Rcom))
    except FileExistsError:
        pass
    plt.savefig("results/" + graph.name + "_" + str(graph.k) + str(graph.Rcapt) + str(graph.Rcom) + "/"
                + solver.__name__ + "_stats" + str(n) + ".jpg")
    plt.clf()
    plt.close()


if __name__ == "__main__":

    filenames = [
        "grille1010_1.dat",
        # "grille1010_2.dat",
        # "captANOR150_7_4.dat",
        # "grille1515_1.dat",
        # "grille1515_2.dat",
        # "grille2020_1.dat",
        # "grille2020_2.dat",
        # "grille2525_1.dat",
        # "grille2525_2.dat",
        # "grille3030_1.dat",
        # "grille3030_2.dat",
        # "grille4040_1.dat",
        # "grille4040_2.dat",
        # "captANOR400_7_10_2021.dat",
        # "captANOR900_14_20_2021.dat",
        # "captANOR1600_16_100_2021.dat"
    ]

    for file in filenames:
        basic_graph = Graph(file)
        for k in [1] if file.startswith("grille") else [1, 2, 3]:
            for Rcom in [1, 2, 3]:
                for Rcapt in {1: [1], 2: [1, 2], 3: [2]}[Rcom]:
                    # graph_ = instance_graph(basic_graph, k, Rcapt, Rcom)
                    for algorithm in [
                        instance_graph.tabou]:  # , instance_graph.simulated_annealing, instance_graph.naive_reduction]:
                        graph_ = instance_graph(basic_graph, k, Rcapt, Rcom)
                        full_test(graph_, solver=algorithm)
                        # make_stats(graph, 10, algorithm)
