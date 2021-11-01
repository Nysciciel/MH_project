import matplotlib.pyplot as plt
from math import sqrt, exp, ceil
from random import random, sample, expovariate, shuffle
from time import time
from os import mkdir
from collections import deque
from line_profiler_pycharm import profile


class Target:
    # class for a target on the graph
    # initialised with its coordinates

    def __init__(self, i, j):
        self.i = i
        self.j = j

    def __repr__(self):
        return "(" + str(self.i) + "," + str(self.j) + ")"

    def __hash__(self):
        return hash((self.i, self.j))

    def __eq__(self, other):
        return self.i == other.i and self.j == other.j

    # computes the distance between himself and some other target
    def __sub__(self, other):
        return sqrt((self.i - other.i) ** 2 + (self.j - other.j) ** 2)


class Graph:
    # class for a graph, without the information of k Rcapt or Rcom
    # takes as argument the name of the file of the data, located in data/
    # example: "grille1010_1.dat"

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
    # class for the actual graph with all the information needed for the problem
    # every algorithm to solve the problem is a method of this
    # takes as argument a Graph from earlier as well as k Rcapt and Rcom

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
        plt.clf()
        plt.close()

    def evaluate(self, sensor_encoding, verbose=False):
        # method to evaluate a solution
        # takes a set of Targets, which are the sensors of the solution
        # returns :
        # a boolean, true if the solution is valid, false otherwise
        # the cost, as the number of sensors divided by the total number of points such that it's between 0 and 1
        # and other indicators of the quality of the solution, between 0 and 1, the higher, the better:
        # the degree, corresponding to the k-coverage
        # and the connexity corresponding to the sensor communication
        # also the verbose can be set to true to get a graph representation of the solution

        if (tuple(sensor_encoding) in self.eval_dict) and not verbose:
            return self.eval_dict[tuple(sensor_encoding)]
        sensor_dict = {target: target in sensor_encoding for target in self.points}
        sensors = {point if sensor_dict[point] else None for point in self.points}
        sensors.discard(None)
        k_dict = {target: (1 if target in sensors else 0) for target in self.points}
        k_dict[Target(0, 0)] = self.k

        visited = {Target(0, 0)}
        to_visit = self.neighborhood_dict[Target(0, 0)].intersection(sensors)
        if verbose:
            for other in to_visit:
                plt.plot([0, other.i], [0, other.j], 'g')
        while len(to_visit) != 0:
            target = to_visit.pop()
            visited.add(target)
            for other in self.neighborhood_dict[target]:
                if (self.Rcapt == self.Rcom or target - other <= self.Rcapt) and k_dict[other] < self.k:
                    k_dict[other] += 1
                if other in sensors and other not in visited:
                    to_visit.add(other)
                    if verbose:
                        plt.plot([target.i, other.i], [target.j, other.j], 'g')

        if verbose:
            sensors = {point if sensor_dict[point] else None for point in self.points}
            sensors.discard(None)
            plt.plot([point.i for point in self.points], [point.j for point in self.points], 'r.')
            plt.plot([point.i for point in sensors], [point.j for point in sensors], 'bo')
            plt.plot([0], [0], 'bd')
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

        connected = all({k_dict[target] >= self.k for target in self.points})
        degree = sum(k_dict.values()) / (self.k * len(self.points)) if not connected else 1
        connexity = 1 - len(sensors - visited)
        feasible = (connexity == 1) and connected
        self.eval_dict[tuple(sensor_encoding)] = feasible, len(sensors) / len(self.points), degree, connexity

        return feasible, len(sensors) / len(self.points), degree, connexity

    def simulated_annealing(self):

        epsilon = 0.0001
        q = 0.95
        energy_function = lambda feasible, cost, degree, connexity: \
            3 * (1 - feasible) + \
            2 * (1 - degree) + \
            2 * (1 - connexity) + \
            cost

        current_temperature = 20
        temperature_decrease_function = lambda T: T * q + (1 - q) * (epsilon / 5)
        stopping_criteria = lambda T: T < epsilon
        iterations_per_temperature = len(self.points)
        sample_size = 1

        def next_solution(solution, lambd):
            n_ = len(self.points)
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

        while not stopping_criteria(current_temperature):
            for _ in range(iterations_per_temperature):
                new_solutions = [next_solution(current_solution, ceil(current_temperature / epsilon) / 10) for __ in
                                 range(sample_size)]
                new_evals = {tuple(i): self.evaluate(i) for i in new_solutions}
                new_solution = min(new_evals, key=lambda x: energy_function(*new_evals.get(x)))
                new_eval = new_evals[new_solution]
                new_solution = set(new_solution)
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

    def naive_reduction(self):

        current_solution = self.points.copy()
        points = list(current_solution.copy())
        shuffle(points)
        for point in points:
            # if random() < 1 / 10:
            #     return current_solution
            if self.evaluate(current_solution - {point})[0]:
                current_solution.remove(point)
        return current_solution

    def tabou(self):
        duration = 600 * len(self.points) / 1500

        tabou_length = len(self.points)

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
        old_cost = cost_function(*self.evaluate(current_solution))

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
            best_cost_neighbor = float('inf')
            new_cost = old_cost
            for neighbor in neighbors:
                neighbor = set(neighbor)
                eval_ = self.evaluate(neighbor)
                new_cost = cost_function(*eval_)
                if new_cost < best_cost_neighbor:
                    current_solution = neighbor
                    best_cost_neighbor = new_cost
                    if new_cost < best_cost:
                        it = 0
                        best_solution = current_solution
                        best_cost = new_cost
            if random() < 1 - 1 / (abs(best_cost_neighbor - old_cost) + 1):
                tabou.append(current_solution)
            old_cost = new_cost

        old_best_solution = best_solution.copy()
        for target in old_best_solution:
            new_eval = self.evaluate(best_solution - {target})
            if new_eval[0]:
                # print("late optimization")
                best_solution -= {target}
        return best_solution

    def random_embrace(self):

        duration = (600 * len(self.points) / 1500 * 9) // 10
        cost_function = lambda feasible, cost_, degree, connexity: \
            3 * (1 - feasible) + \
            2 * (1 - degree) + \
            2 * (1 - connexity) + \
            cost_
        t0 = time()
        best_solution = None
        best_cost = float('inf')
        it = 0
        while time() - t0 < duration:
            it += 1
            solution = self.naive_reduction()
            cost = cost_function(*self.evaluate(solution))
            if cost < best_cost:
                best_cost = cost
                best_solution = solution
        return best_solution

    def anneal_tabou(self):
        epsilon = 0.0001
        q = 0.90
        energy_function = lambda feasible, cost, degree, connexity: \
            3 * (1 - feasible) + \
            2 * (1 - degree) + \
            2 * (1 - connexity) + \
            cost

        tabou_length = len(self.points) // 10
        tabou = deque([], tabou_length)
        current_temperature = 20
        temperature_decrease_function = lambda T: T * q + (1 - q) * (epsilon / 4)
        stopping_criteria = lambda T: T < epsilon
        iterations_per_temperature = len(self.points)

        def next_solution(solution, lambd, pool_):
            n_ = len(pool_)
            dist = ceil(expovariate(1 / lambd))
            while dist > n_:
                dist = ceil(expovariate(1 / lambd))
            points = sample(pool_, dist)
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

        while not stopping_criteria(current_temperature):
            for _ in range(iterations_per_temperature):
                new_solution = next_solution(current_solution, ceil(current_temperature / epsilon) / 10,
                                             list(self.points - set(tabou)))
                new_eval = self.evaluate(new_solution)
                new_solution = set(new_solution)
                new_energy = energy_function(*new_eval)
                energy_delta = new_energy - current_energy
                p = random()
                thresh = exp(-energy_delta / current_temperature)
                if energy_delta < 0 or p < thresh:
                    if energy_delta >= 0 and p < thresh:
                        tabou.extend(current_solution.symmetric_difference(new_solution))
                        # points = list(current_solution.copy())
                        # shuffle(points)
                        # for point in points:
                        #     eval_ = self.evaluate(current_solution - {point})
                        #     if eval_[0]:
                        #         current_solution.remove(point)
                        #         tabou.append(current_solution)
                        #         current_energy = energy_function(*eval_)
                        #         if current_energy < best_score and eval_[0]:
                        #             best_score = current_energy
                        #             best_solution = current_solution
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


def full_test(graph: instance_graph, solver=instance_graph.simulated_annealing):
    # runs a test of the solver on the graph
    # and plots out the graph
    # and measures the time taken
    t0 = time()
    print("\nTest for", graph.name, "k=", graph.k, "Rcapt=", graph.Rcapt, "Rcom=", graph.Rcom,
          "algorithm=", solver.__name__, )
    sol = solver(graph)
    print("result:", len(sol), "time:", time() - t0)
    graph.evaluate(sol, True)


def make_stats(graph: instance_graph, n, solver=instance_graph.simulated_annealing):
    # runs n times the algorithm on the graph, and plots an histogram of the results
    # stored in results/ "graph name" / " k, Rcapt, Rcom " / "algorithm name _ stats" n
    print("\nStats for", graph.name, "k=", graph.k, "Rcapt=", graph.Rcapt, "Rcom=", graph.Rcom,
          "Times=", n, "algorithm=", solver.__name__)
    t0 = time()
    plt.hist([len(solver(graph)) for _ in range(n)], density=True)
    plt.title("Avg time:" + str((time() - t0) / n))
    try:
        mkdir("results")
    except FileExistsError:
        pass
    try:
        mkdir("results/" + graph.name)
    except FileExistsError:
        pass
    try:
        mkdir("results/" + graph.name + "/" + str(graph.k) + str(graph.Rcapt) + str(graph.Rcom))
    except FileExistsError:
        pass
    plt.savefig("results/" + graph.name + "/" + str(graph.k) + str(graph.Rcapt) + str(graph.Rcom) + "/"
                + solver.__name__ + "_stats" + str(n) + ".jpg")
    plt.clf()
    plt.close()


if __name__ == "__main__":

    filenames = [
        "grille1010_1.dat",
        # "grille1010_2.dat",
        "grille1515_1.dat",
        # "grille1515_2.dat",
        # "grille2020_1.dat",
        # "grille2020_2.dat",
        # "grille2525_1.dat",
        # "grille2525_2.dat",
        # "grille3030_1.dat",
        # "grille3030_2.dat",
        # "grille4040_1.dat",
        # "grille4040_2.dat",
        "captANOR150_7_4.dat",
        # "captANOR400_7_10_2021.dat",
        # "captANOR900_14_20_2021.dat",
        # "captANOR1600_16_100_2021.dat"
    ]

    for file in filenames:
        print("")
        basic_graph = Graph(file)
        for k in [1] if file.startswith("grille") else [1, 2, 3]:
            for Rcom in [1, 2, 3]:
                for Rcapt in {1: [1], 2: [1, 2], 3: [2]}[Rcom]:
                    print("")
                    graph_ = instance_graph(basic_graph, k, Rcapt, Rcom)
                    for algorithm in [
                        instance_graph.anneal_tabou,
                        instance_graph.tabou,
                        instance_graph.simulated_annealing,
                        instance_graph.random_embrace,
                        instance_graph.naive_reduction
                    ]:
                        full_test(graph_, solver=algorithm)
                        # make_stats(graph_, 10, algorithm)
