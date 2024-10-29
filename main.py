# main.py
from border.border import Border
from grid.grid import Grid

# Example usage
border = Border(30, 44, 68, 70, 100, 4, 10)
x = 65
print("Upper:", border.get_upper(x))
print("Lower:", border.get_lower(x))
print("Distance:", border.get_distance(x))

# Create Grid instance and print calcX results
grid = Grid(border)
x_values = grid.calcX()
print("x_values:", x_values)

# grid.plotX()

# grid.calcY(65)
# print("y_values:", grid.calcY(65))
# grid.plotY(65)

grid.plotYForAllX()
