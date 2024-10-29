import numpy as np
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
        lambda_ = 1.1
        ys.append(y)
        dy = dy0
        while y < y_star:
            y += dy
            ys.append(y)
            dy *= lambda_
        # New y star is the end of the list
        y_star_upper = ys[-1]
        # Calculate nodes for the upper part using a new y_star_upper
        y = self.border.get_upper(x)
        dy = dy0
        upper_ys = []
        while y > y_star_upper:
            upper_ys.append(y)
            dy *= lambda_
            y -= dy
            # delete the last element if the difference with y_star_upper is less 0.2
            if abs(y - y_star_upper) < 0.05:
                upper_ys.pop()
        # Append the upper part nodes in reverse order
        ys.extend(reversed(upper_ys))
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
        plt.xlabel("x value")
        plt.ylabel("y value")
        plt.title("1D Plot of y values for all x values")
        plt.legend()
        plt.grid(True)
        plt.savefig("grid.svg")
        plt.show()
