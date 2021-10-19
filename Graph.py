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
            self.points.remove(None)
        if len(self.points) == 0:
            raise Exception("Parsing " + filename + " failed")
        self.points.add(Target(0, 0))
        self.n = n
        self.m = m

    def __repr__(self):
        return self.points.__repr__()

    def plot(self):
        for y in range(self.m):
            for x in range(self.n):
                if Target(x, y) in self.points:
                    print("°", end="")
                else:
                    print(" ", end="")
            print("|", end="\n")

    def true_plot(self):
        import matplotlib.pyplot as plt
        plt.plot([point.i for point in self.points], [point.j for point in self.points], 'r.')
        plt.gca().invert_yaxis()
        plt.show()

    def evaluate(self, k, Rcom, Rcapt, sensor_encoding):
        sensor_dict = {target: target in sensor_encoding for target in self.points}

        k_dict = {target: (1 if sensor_dict[target] else 0) for target in self.points}
        k_dict[Target(0, 0)] = k
        visited = set()
        connexes = {(point,) for point in self.points}
        res = 0
        for target in self.points:
            visited.add(target)
            if sensor_dict[target]:
                res += 1
            for other in self.points - visited:
                d = target - other
                if d < Rcom:
                    target_connexe = None
                    other_connexe = None
                    for connexe in connexes:
                        if target in connexe:
                            target_connexe = connexe
                        if other in connexe:
                            other_connexe = connexe
                    if target_connexe != other_connexe:
                        connexes.remove(target_connexe)
                        connexes.remove(other_connexe)
                        connexes.add(target_connexe + other_connexe)
                    if d < Rcapt:
                        if sensor_dict[other]:
                            k_dict[target] += 1
                        if sensor_dict[target]:
                            k_dict[other] += 1

        print(connexes)
        print(k_dict)

        import matplotlib.pyplot as plt
        sensors = {point if sensor_dict[point] else None for point in self.points}
        sensors.remove(None)
        plt.plot([point.i for point in self.points], [point.j for point in self.points], 'r.')
        plt.plot([point.i for point in sensors], [point.j for point in sensors], 'bo')
        plt.gca().invert_yaxis()
        plt.show()

        if len(connexes) != 1:
            return False, res
        for target in self.points:
            if k_dict[target] < k:
                return False, res
        return True, res


if __name__ == "__main__":
    g1 = Graph("grille1010_1.dat")
    g2 = Graph("grille1010_2.dat")
    g3 = Graph("captANOR150_7_4.dat")
    g4 = Graph("test.dat")
    print("")
    print(g4.evaluate(1, 4, 4, [Target(0.01, 2.15), Target(2.86, 5.41)]))
