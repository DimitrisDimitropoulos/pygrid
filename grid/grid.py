from scipy.optimize import newton
from border.border import Border
import matplotlib.pyplot as plt


class Grid:
    def __init__(self, border):
        """
        Constructor for the Grid class.
        :param border: A Border object
        :param num_points: Number of points in the grid
        """
        self.border = border

    def _calculate_lambda(self, x, d0, lambda_initial, nodes):
        """
        Private method to calculate lambda for a given x, d0, and initial lambda
        by solving a nonlinear equation using the Newton-Raphson algorithm.
        :param x: The x position
        :param d0: Initial distance or spacing
        :param lambda_initial: Initial guess for lambda
        :return: Calculated value of lambda
        """
        y_len = self.border.get_distance(x)

        # Define the nonlinear equation to solve
        def equation(lambda_):
            return (
                d0 * (lambda_ ** ((nodes * 0.5) - 1))
                - 0.5 * y_len * lambda_
                + 0.5 * y_len
                - d0
            )

        # Use Newton-Raphson method to find the root
        lambda_solution = newton(equation, lambda_initial, maxiter=100000, tol=1e-12)
        return lambda_solution

    def calcX(self):
        """
        Calculate a vector of x-values from 0 to xMax with num_points.
        :return: A numpy array of x-values
        """
        xs = []
        dx0 = 0.5
        lambda_ = 1.1
        # Prepend values using geometric progression until x is greater than 0
        x = self.border.x1
        dx = dx0
        while x > 0:
            xs.insert(0, x)
            dx *= lambda_
            x -= dx
        # Reset x to the initial value
        x = self.border.x1
        # Append values with equal spacing dx0
        while x <= self.border.x4:
            xs.append(x)
            x += dx0
        # Append values using geometric progression until x exceeds xMax
        dx = dx0
        while x <= self.border.xMax:
            xs.append(x)
            dx *= lambda_
            x += dx
        # Append values greater than xMax using the same geometric progression
        while x <= self.border.xMax * lambda_:
            xs.append(x)
            dx *= lambda_
            x += dx
        # make first at 0 and last at xMax
        # NOTE: this is a manipulation of the original code
        xs[0] = 0
        xs[-1] = self.border.xMax
        return xs

    def plotX(self):
        """
        Plot the vector of x-values as points.
        """
        xs = self.calcX()
        ys = [0] * len(xs)  # y-values set to zero for a 1D plot
        plt.figure(figsize=(10, 6))
        plt.scatter(xs, ys, label="x values")
        plt.xlabel("x value")
        plt.title("1D Plot of x values")
        plt.legend()
        plt.grid(True)
        plt.show()

    def calcY(self, x):
        ys = []
        y_len = self.border.get_distance(x)
        y_star = self.border.get_upper(x) - y_len / 2
        y = self.border.get_lower(x)
        dy0 = 10**-4
        lambda_initial = 1.1  # Initial guess for lambda
        # Calculate lambda using the private method
        lambda_ = self._calculate_lambda(x, dy0, lambda_initial, 180)
        ys.append(y)
        dy = dy0
        while y < y_star:
            y += dy
            ys.append(y)
            dy *= lambda_
        # Delete the last value and replace it with y_star
        ys = ys[:-1]
        ys.append(y_star)
        # Calculate the symmetric points against y_star
        # ys = ys[:-1]  # remove the last value
        upper_part = [2 * y_star - y for y in reversed(ys)]
        # Append the upper part nodes
        ys.extend(upper_part)
        # print the length the x and lambda
        # print(f"X: {x}, YLen: {len(ys)}, Lambda: {lambda_}")
        return ys

    def plotY(self, x):
        """
        Plot the vector of y-values as points on the Y axis.
        """
        ys = self.calcY(x)
        xs = [0] * len(ys)  # x-values set to zero for a 1D plot
        plt.figure(figsize=(10, 6))
        plt.scatter(xs, ys, label="y values")
        plt.ylabel("y value")
        plt.title("1D Plot of y values")
        plt.legend()
        plt.grid(True)
        plt.show()

    def xgrid(self, i, j):
        """
        Return the x-value at grid point (i, j).
        Since x-values do not depend on j, the j parameter is not used.
        :param i: Index along the x-axis
        :param j: Index along the y-axis (not used)
        :return: x-value at grid point (i, j)
        """
        xs = self.calcX()
        return xs[i]

    def ygrid(self, i, j):
        """
        Return the y-value at grid point (i, j).
        The y-values depend on x, so we use x index i to get ys at that x.
        :param i: Index along the x-axis
        :param j: Index along the y-axis
        :return: y-value at grid point (i, j)
        """
        x = self.calcX()[i]
        ys = self.calcY(x)
        return ys[j]

    def plotYForAllX(self):
        """
        Plot the vector of y-values for all x-values in the same plot.
        """
        xs = self.calcX()
        plt.figure(figsize=(10, 6))
        for x in xs:
            try:
                ys = self.calcY(x)
                x_values = [x] * len(ys)  # x-values set to the current x for each y
                plt.scatter(x_values, ys, color="black", s=1)
            except ValueError as e:
                print(f"Skipping x={x}: {e}")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Nodes of the grid")
        plt.legend()
        plt.grid(True)
        plt.savefig("./out/grid.svg")
        plt.show()
