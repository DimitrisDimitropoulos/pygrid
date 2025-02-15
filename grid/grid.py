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
        :return: A list of x-values
        """
        xs = []
        seen = set()
        dx0 = 0.5
        # lambda_ = 1.07
        lambda_ = 1.1
        # Prepend values using geometric progression until x is greater than 0
        x = self.border.x1
        dx = dx0
        while x > self.border.xMin:
            if x not in seen:
                xs.insert(0, x)
                seen.add(x)
            dx *= lambda_
            x -= dx
        # Reset x to the initial value
        x = self.border.x1
        # Append values with equal spacing dx0
        while x <= self.border.x4:
            if x not in seen:
                xs.append(x)
                seen.add(x)
            x += dx0
        # Append values using geometric progression until x exceeds xMax
        dx = dx0
        while x <= self.border.xMax:
            if x not in seen:
                xs.append(x)
                seen.add(x)
            dx *= lambda_
            x += dx
        # Append values greater than xMax using the same geometric progression
        while x <= self.border.xMax:
            if x not in seen:
                xs.append(x)
                seen.add(x)
            dx *= lambda_
            x += dx
        # make first at 0 and last at xMax
        # NOTE: this is a manipulation of the original code
        xs[0] = self.border.xMin
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
        seen = set()
        y_len = self.border.get_distance(x)
        y_star = self.border.get_upper(x) - y_len / 2
        y = self.border.get_lower(x)
        dy0 = 10**-4
        lambda_initial = 1.1  # Initial guess for lambda
        # Calculate lambda using the private method
        lambda_ = self._calculate_lambda(x, dy0, lambda_initial, 181)
        if y not in seen:
            ys.append(y)
            seen.add(y)
        dy = dy0
        while y < y_star:
            y += dy
            if y not in seen:
                ys.append(y)
                seen.add(y)
            dy *= lambda_
        # Delete the last value and replace it with y_star
        ys = ys[:-1]
        ys.append(y_star)
        # Calculate the symmetric points against y_star
        ys = ys[:-1]  # remove the last value
        upper_part = [
            2 * y_star - y for y in reversed(ys) if (2 * y_star - y) not in seen
        ]
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

    def plotGrid(self):
        """
        Plot the vector of y-values for all x-values in the same plot,
        connecting them to appear as a grid.
        """
        xs = self.calcX()
        plt.figure(figsize=(10, 6))
        # Plot points and lines to form a grid
        for x in xs:
            try:
                ys = self.calcY(x)
                x_values = [x] * len(ys)  # x-values set to the current x for each y
                # Plot individual nodes
                # plt.scatter(x_values, ys, color="black", s=1)
                # Connect nodes in a row along the y-axis
                plt.plot([x] * len(ys), ys, color="black", linewidth=0.1)
            except ValueError as e:
                print(f"Skipping x={x}: {e}")
        # Connect nodes in columns along the x-axis
        for i in range(len(xs) - 1):
            try:
                ys_current = self.calcY(xs[i])
                ys_next = self.calcY(xs[i + 1])
                # Only connect nodes if both y-lists are of the same length
                if len(ys_current) == len(ys_next):
                    for y1, y2 in zip(ys_current, ys_next):
                        plt.plot(
                            [xs[i], xs[i + 1]], [y1, y2], color="black", linewidth=0.1
                        )
            except ValueError as e:
                print(f"Skipping connection between x={xs[i]} and x={xs[i + 1]}: {e}")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Connected Nodes of the Grid")
        # plt.grid(True)
        plt.savefig("./out/connected_grid.svg")
        plt.show()

    def generate_grid_ascii(self, imax, jmax, xgrid, ygrid):
        # Open the file 'grid.ascii' for writing
        file = open("grid.ascii", "w")

        # -----
        # Nodes
        # -----
        nindx = 2
        nd = 2  # Number of dimensions
        file.write("{0:1d} {1:1d}\n".format(nindx, nd))

        nindx = 10  # Index, =10 for nodes, =12 for cells, =13 for faces
        file.write("{0:2d}\n".format(nindx))

        icol = 6  # Number of columns
        izoneid = 0  # Number of zone, defines if the face is interior or boundary
        NNT_start = 1  # Starting node
        NNT_end = imax * jmax  # Ending node
        itype = 0
        node_nd = 2
        file.write(
            "{0:1d} {1:2d} {2:1d} {3:1d} {4:8d} {5:1d} {6:1d}\n".format(
                icol, nindx, izoneid, NNT_start, NNT_end, itype, node_nd
            )
        )

        izoneid = 8
        itype = 1
        file.write("{0:2d}\n".format(nindx))
        file.write(
            "{0:1d} {1:2d} {2:1d} {3:1d} {4:8d} {5:1d} {6:1d}\n".format(
                icol, nindx, izoneid, NNT_start, NNT_end, itype, node_nd
            )
        )

        for j in range(1, jmax + 1):
            for i in range(1, imax + 1):
                file.write(
                    "{0:23.16f} {1:23.16f}\n".format(
                        xgrid[i - 1][j - 1], ygrid[i - 1][j - 1]
                    )
                )

        # -----
        # Cells
        # -----
        nindx = 12
        izoneid = 0
        itype = 0
        nte_type = 0
        NTE_start = 1
        NTE_end = (imax - 1) * (jmax - 1)
        file.write("{0:2d}\n".format(nindx))
        file.write(
            "{0:1d} {1:2d} {2:1d} {3:1d} {4:8d} {5:1d} {6:1d}\n".format(
                icol, nindx, izoneid, NTE_start, NTE_end, itype, nte_type
            )
        )

        izoneid = 9
        itype = 1
        nte_type = 3
        file.write("{0:2d}\n".format(nindx))
        file.write(
            "{0:1d} {1:2d} {2:1d} {3:1d} {4:8d} {5:1d} {6:1d}\n".format(
                icol, nindx, izoneid, NTE_start, NTE_end, itype, nte_type
            )
        )

        # -----
        # Faces
        # -----
        nindx = 13
        izoneid = 0
        itype = 0
        NBC_type = 0
        NFace_type = 0
        nte_type = 0

        # -- Calculate number of faces
        # -- Interior faces
        nfaces_y = (imax - 2) * (jmax - 1)  # y normal faces
        nfaces_x = (imax - 1) * (jmax - 2)  # x normal faces
        nfaces_int = nfaces_x + nfaces_y

        # -- Boundary faces
        nfaces_ilo = jmax - 1
        nfaces_ihi = jmax - 1
        nfaces_jlo = imax - 1
        nfaces_jhi = imax - 1  # Number of faces on jhi

        nfaces_tot = nfaces_int + nfaces_ilo + nfaces_ihi + nfaces_jlo + nfaces_jhi

        NFace_start = 1
        NFace_end = nfaces_tot
        icol = 6  # Number of columns
        file.write("{0:2d}\n".format(nindx))
        file.write(
            "{0:1d} {1:2d} {2:1d} {3:1d} {4:8d} {5:1d} {6:1d}\n".format(
                icol, nindx, izoneid, NFace_start, NFace_end, NBC_type, NFace_type
            )
        )

        # ----
        # Interior faces
        # ----
        izoneid = 10
        NBC_type = 2
        NFace_type = 2
        NFace_start = 1
        NFace_end = nfaces_int
        icol = 6  # Number of columns
        file.write("{0:2d}\n".format(nindx))
        file.write(
            "{0:1d} {1:2d} {2:2d} {3:1d} {4:8d} {5:1d} {6:1d}\n".format(
                icol, nindx, izoneid, NFace_start, NFace_end, NBC_type, NFace_type
            )
        )

        icol = 4  # Number of columns
        for j in range(1, jmax - 1):
            for i in range(2, imax):
                # vertical face
                nod1 = (j - 1) * imax + i
                nod2 = nod1 + imax
                nel1 = (j - 1) * (imax - 1) + i - 1
                nel2 = nel1 + 1
                file.write(
                    "{0:1d} {1:7d} {2:7d} {3:7d} {4:7d}\n".format(
                        icol, nod1, nod2, nel1, nel2
                    )
                )
                # horizontal face
                nod1 = nod2
                nod2 = nod1 - 1
                nel1 = nel1
                nel2 = nel1 + (imax - 1)
                file.write(
                    "{0:1d} {1:7d} {2:7d} {3:7d} {4:7d}\n".format(
                        icol, nod1, nod2, nel1, nel2
                    )
                )
            i = imax
            # horizontal face
            nod1 = j * imax + i
            nod2 = nod1 - 1
            nel1 = (j - 1) * (imax - 1) + i - 1
            nel2 = nel1 + (imax - 1)
            file.write(
                "{0:1d} {1:7d} {2:7d} {3:7d} {4:7d}\n".format(
                    icol, nod1, nod2, nel1, nel2
                )
            )

        j = jmax - 1
        for i in range(2, imax):
            # vertical face
            nod1 = (j - 1) * imax + i
            nod2 = nod1 + imax
            nel1 = (j - 1) * (imax - 1) + i - 1
            nel2 = nel1 + 1
            file.write(
                "{0:1d} {1:7d} {2:7d} {3:7d} {4:7d}\n".format(
                    icol, nod1, nod2, nel1, nel2
                )
            )

        # ----
        # j=jmax boundary - wall
        # ----
        izoneid = 4
        NBC_type = 3
        NFace_type = 2
        NFace_start = NFace_end + 1
        NFace_end = NFace_start + nfaces_jhi - 1
        icol = 6  # Number of columns
        file.write("{0:2d}\n".format(nindx))
        file.write(
            "{0:1d} {1:2d} {2:2d} {3:8d} {4:8d} {5:2d} {6:1d}\n".format(
                icol, nindx, izoneid, NFace_start, NFace_end, NBC_type, NFace_type
            )
        )

        icol = 4  # Number of columns
        for i in range(1, imax):
            nod1 = (jmax - 1) * imax + i + 1
            nod2 = nod1 - 1
            nel1 = (jmax - 2) * (imax - 1) + i
            nel2 = 0
            file.write(
                "{0:1d} {1:7d} {2:7d} {3:7d} {4:7d}\n".format(
                    icol, nod1, nod2, nel1, nel2
                )
            )

        # ----
        # j=1 boundary - wall
        # ----
        izoneid = 4
        NBC_type = 3
        NFace_type = 2
        NFace_start = NFace_end + 1
        NFace_end = NFace_start + nfaces_jlo - 1
        icol = 6  # Number of columns
        file.write("{0:2d}\n".format(nindx))
        file.write(
            "{0:1d} {1:2d} {2:2d} {3:8d} {4:8d} {5:2d} {6:1d}\n".format(
                icol, nindx, izoneid, NFace_start, NFace_end, NBC_type, NFace_type
            )
        )

        icol = 4  # Number of columns
        for i in range(1, imax):
            nod1 = i
            nod2 = i + 1
            nel1 = i
            nel2 = 0
            file.write(
                "{0:1d} {1:7d} {2:7d} {3:7d} {4:7d}\n".format(
                    icol, nod1, nod2, nel1, nel2
                )
            )

        # ----
        # i=imax boundary - Outflow
        # ----
        izoneid = 5
        NBC_type = 0
        NFace_type = 2
        NFace_start = NFace_end + 1
        NFace_end = NFace_start + nfaces_ihi - 1
        icol = 6  # Number of columns
        file.write("{0:2d}\n".format(nindx))
        file.write(
            "{0:1d} {1:2d} {2:2d} {3:8d} {4:8d} {5:2d} {6:1d}\n".format(
                icol, nindx, izoneid, NFace_start, NFace_end, NBC_type, NFace_type
            )
        )

        icol = 4  # Number of columns
        for j in range(1, jmax):
            nod1 = j * imax
            nod2 = (j + 1) * imax
            nel1 = j * (imax - 1)
            nel2 = 0
            file.write(
                "{0:1d} {1:7d} {2:7d} {3:7d} {4:7d}\n".format(
                    icol, nod1, nod2, nel1, nel2
                )
            )

        # ----
        # i=1 boundary - Inflow
        # ----
        izoneid = 4
        NBC_type = 1
        NFace_type = 2
        NFace_start = NFace_end + 1
        NFace_end = NFace_start + nfaces_ilo - 1
        icol = 6  # Number of columns
        file.write("{0:2d}\n".format(nindx))
        file.write(
            "{0:1d} {1:2d} {2:2d} {3:8d} {4:8d} {5:2d} {6:1d}\n".format(
                icol, nindx, izoneid, NFace_start, NFace_end, NBC_type, NFace_type
            )
        )

        icol = 4  # Number of columns
        for j in range(1, jmax):
            nod1 = j * imax + 1
            nod2 = (j - 1) * imax + 1
            nel1 = (j - 1) * (imax - 1) + 1
            nel2 = 0
            file.write(
                "{0:1d} {1:7d} {2:7d} {3:7d} {4:7d}\n".format(
                    icol, nod1, nod2, nel1, nel2
                )
            )

        file.close()
