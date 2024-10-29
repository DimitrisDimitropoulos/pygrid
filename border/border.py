# border/border.py
import math


class Border:
    def __init__(self, x1, x2, x3, x4, xMax, y1, y2):
        """
        Constructor for the Border class.
        :param x1: x1 coordinate
        :param x2: x2 coordinate
        :param x3: x3 coordinate
        :param x4: x4 coordinate
        :param xMax: maximum x-coordinate
        :param y1: upper y-coordinate
        :param y2: lower y-coordinate
        """
        self.x1 = x1
        self.x2 = x2
        self.x3 = x3
        self.x4 = x4
        self.xMax = xMax
        self.y1 = y1
        self.y2 = y2

    def get_upper(self, x):
        """
        Get the y-coordinate of the upper border at a given x-coordinate.
        :param x: The x-coordinate
        :return: The y-coordinate of the upper border
        """
        return self.y1

    def get_lower(self, x):
        """
        Get the y-coordinate of the lower border at a given x-coordinate.
        :param x: The x-coordinate
        :return: The y-coordinate of the lower border
        """
        a3 = 8.68055 * 1e-4
        a2 = -0.145834
        a1 = 7.79167
        a0 = -134.445
        if x < self.x2 and x >= 0:
            return 0
        elif self.x2 <= x <= self.x3:
            return a3 * x**3 + a2 * x**2 + a1 * x + a0
        elif self.x3 < x < self.xMax:
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
