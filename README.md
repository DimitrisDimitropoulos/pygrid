# CFD Grid Project

This project generates computational grids for fluid dynamics simulations using
custom-defined border and grid classes. The primary driver of the application
is the `main.py` file, which creates a `Border` and a `Grid` object, computes
the grid points, plots the grid, and finally outputs an ASCII representation of
the grid.

## Key Components

- **Border Class (`border/border.py`):**  
  Defines the geometric boundary of the domain using a set of x and y coordinates. It provides methods to get the upper and lower borders as well as the distance between them at a given x-coordinate.

- **Grid Class (`grid/grid.py`):**  
  Generates x- and y-coordinates for the grid based on the defined border. It includes methods to:

    - Compute and return a vector of x-values (`calcX`).
    - Compute and return a vector of y-values for a given x (`calcY`).
    - Plot x-values, y-values, and the entire grid.
    - Export the grid to an ASCII file.

- **Main Script (`main.py`):**  
  Serves as the entry point to the application. It:
    - Instantiates a `Border` object with appropriate coordinates.
    - Uses the Border object to create a `Grid` object.
    - Prints sample values (upper, lower border, distance).
    - Computes and prints the generated x-values.
    - Generates the grid and plots it.
    - Finally, exports the grid structure to an ASCII file.

## Requirements

- Python 3.x
- [NumPy](https://numpy.org/)
- [SciPy](https://www.scipy.org/)
- [Matplotlib](https://matplotlib.org/)

You can install the required packages using pip:

```bash
pip install numpy scipy matplotlib

```

