# border/border.py
import numpy as np


class Border:
    def __init__(self, x1, x2, x3, x4, xMax, y1, y2, xMin):
        """
        Constructor for the Border class.
        :param x1: x1 coordinate
        :param x2: x2 coordinate
        :param x3: x3 coordinate
        :param x4: x4 coordinate
        :param xMax: maximum x-coordinate
        :param y1: upper y-coordinate
        :param y2: lower y-coordinate
        :param xMin: lower x bound
        """
        self.x1 = x1
        self.x2 = x2
        self.x3 = x3
        self.x4 = x4
        self.xMax = xMax
        self.y1 = y1
        self.y2 = y2
        self.xMin = xMin

    def get_upper(self, x):
        """
        Get the y-coordinate of the upper border at a given x-coordinate.
        :param x: The x-coordinate
        :return: The y-coordinate of the upper border
        """
        return self.y1

    def _calculate_coefficients(self):
        """
        Private method to calculate the coefficients a3, a2, a1, a0 by solving a linear system.
        Constraints:
            P(x2) = 0
            P(x3) = -abs(y2 - y1)
            P'(x2) = 0
            P'(x3) = 0
        """
        x2 = self.x2
        x3 = self.x3
        y1 = self.y1
        y2 = self.y2
        # Calculate delta_y
        delta_y = -abs(y2 - y1)  # P(x3) = -|y2 - y1|
        # Set up the equations based on the constraints
        # Equation 1: P(x2) = 0
        # Equation 2: P(x3) = delta_y
        # Equation 3: P'(x2) = 0
        # Equation 4: P'(x3) = 0
        # Coefficient matrix A
        A = np.array(
            [
                [x2**3, x2**2, x2, 1],  # P(x2) = 0
                [x3**3, x3**2, x3, 1],  # P(x3) = delta_y
                [3 * x2**2, 2 * x2, 1, 0],  # P'(x2) = 0
                [3 * x3**2, 2 * x3, 1, 0],  # P'(x3) = 0
            ]
        )
        # Right-hand side vector b
        b = np.array(
            [
                0,  # P(x2) = 0
                delta_y,  # P(x3) = -abs(y2 - y1)
                0,  # P'(x2) = 0
                0,  # P'(x3) = 0
            ]
        )
        # Solve the linear system A * [a3, a2, a1, a0]^T = b
        coefficients = np.linalg.solve(A, b)
        return coefficients

    def get_lower(self, x):
        """
        Get the y-coordinate of the lower border at a given x-coordinate.
        :param x: The x-coordinate
        :return: The y-coordinate of the lower border
        """
        a3, a2, a1, a0 = self._calculate_coefficients()
        if x < self.x2 and x >= self.xMin:
            return 0
        elif self.x2 <= x <= self.x3:
            return a3 * x**3 + a2 * x**2 + a1 * x + a0
        elif self.x3 < x <= self.xMax:
            return a3 * self.x3**3 + a2 * self.x3**2 + a1 * self.x3 + a0
        else:
            raise ValueError("x is out of the valid range")

    def get_distance(self, x):
        """
        Get the distance between the upper and lower border at a given x-coordinate.
        :param x: The x-coordinate
        :return: The distance between the upper and lower borders
        """
        return abs(self.get_upper(x)) + abs(self.get_lower(x))
