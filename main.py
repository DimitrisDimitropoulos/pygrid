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

# print("y_values:", grid.calcY(65))
# grid.plotY(65)

# print(grid.xgrid(0, 10))
# print(grid.ygrid(0, 10))

# grid.plotYForAllX()
# grid.plotGrid()

# Generate xgrid and ygrid for the grid
imax: int = len(x_values)
jmax: int = len(grid.calcY(x_values[0]))
xgrid = [[grid.xgrid(i, j) for j in range(jmax)] for i in range(imax)]
ygrid = [[grid.ygrid(i, j) for j in range(jmax)] for i in range(imax)]

# Print the number of elements in ygrid
num_rows = len(ygrid)
print("Number of rows in ygrid:", num_rows)

if num_rows > 0:
    num_columns = len(ygrid[0])
    print("Number of columns in ygrid:", num_columns)
else:
    print("ygrid is empty")


grid.plotGrid()

# Call generate_grid_ascii method
grid.generate_grid_ascii(imax, jmax, xgrid, ygrid)
