class Target:
    def __init__(self, i, j, sensor=False):
        self.i = i
        self.j = j
        self.sensor = sensor

    def __repr__(self):
        return "(" + str(self.i) + "," + str(self.j) + ["", " sensor"][self.sensor] + ")"

    def __hash__(self):
        return tuple(self.__dict__.values()).__hash__()

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()


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
        plt.plot([point.i for point in self.points], [point.j for point in self.points], 'ro')
        plt.gca().invert_yaxis()
        plt.show()


if __name__ == "__main__":
    g1 = Graph("grille1010_1.dat")
    g2 = Graph("grille1010_2.dat")
    g3 = Graph("captANOR150_7_4.dat")
    g1.plot()
    print("")
    g2.plot()
    g1.true_plot()
    g2.true_plot()
    g3.true_plot()
