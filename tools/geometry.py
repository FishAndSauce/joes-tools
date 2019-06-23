from math import acos, degrees, sqrt
from matplotlib import pyplot as plt
import pandas as pd


def CosineRule(a, b, c):
    cos_A = (b ** 2 + c ** 2 - a ** 2) / (2 * b * c)
    return degrees(acos(cos_A))


class Coord(object):
    """ Coordinate representing a position/point in 2D space

    Attributes:
        x (float): Value of x at position
        y (float): Value of y at position
    """

    def __init__(self, x, y, name=None):
        self.x = x
        self.y = y
        self.name = name

    def distance_to(self, other_point):
        rise = other_point.y - self.y
        run = other_point.x - self.x
        return sqrt(rise ** 2 + run ** 2)

    def find_gradient(self, other_coord):
        rise = other_coord.y - self.y
        run = other_coord.x - self.x
        if run == 0.0:
            return None
        else:
            return rise / run

    def find_lowest_alpha(self, B, other_points):
        alphas = []
        for C in other_points:
            c = self.distance_to(C)
            b = self.distance_to(B)
            a = B.distance_to(C)
            beta = CosineRule(a, b, c)
            alphas.append(180 - beta)
        index = alphas.index(min(alphas))
        return other_points[index]

    def convex_quadrant(self, other_point):

        def gradient_direction(value):
            if value > 0.0:
                result = 'increasing'
            elif value == 0.0:
                result = 'flat'
            else:
                result = 'decreasing'
            return result

        Q_matrix = pd.DataFrame(
            [
                [['Q4'], ['Q1'], ['Q4', 'Q1']],
                [['Q3'], ['Q2'], ['Q3', 'Q2']],
                [['Q3', 'Q4'], ['Q1', 'Q2'], []]
            ],
            columns=['increasing', 'decreasing', 'flat'],
            index=['increasing', 'decreasing', 'flat'])

        # Evaluate truth of x and y as increasing from self to other point
        y = gradient_direction(other_point.y - self.y)
        x = gradient_direction(other_point.x - self.x)

        quadrants = Q_matrix.loc[x, y]

        return quadrants

    def plot(self):
        plt.scatter(self.x, self.y)
        # plt.show()


class CoordCollection(object):
    '''

    '''

    def __init__(self, coords):
        self.coords = coords
        self.x_list = [p.x for p in self.coords]
        self.y_list = [p.y for p in self.coords]
        self.name_list = [p.name for p in self.coords]
        self.coords_list = [(p.x, p.y) for p in self.coords]
        self.coords_dict = {coord.name: coord for coord in self.coords}
        self.x_max = max(self.x_list)
        self.envelope = []

    def max_x_point(self):
        index = self.x_list.index(max(self.x_list))
        return (self.coords[index])

    def max_y_point(self):
        index = self.y_list.index(max(self.y_list))
        return (self.coords[index])

    def min_x_point(self):
        index = self.x_list.index(min(self.x_list))
        return (self.coords[index])

    def min_y_point(self):
        index = self.y_list.index(min(self.y_list))
        return (self.coords[index])

    def plot(self):
        plt.scatter(self.x_list, self.y_list)

    def find_envelope_points(self, start_point, last_point):
        other_points = self.coords.copy()
        other_points.remove(start_point)
        next_point = start_point.find_lowest_alpha(
            last_point,
            other_points
        )
        if id(next_point) in [id(x) for x in self.envelope]:
            self.envelope.append(start_point)
            self.envelope.append(next_point)

            return self.envelope
        else:
            self.envelope.append(start_point)
            self.find_envelope_points(next_point, start_point)

    def get_envelope(self):

        start_point = self.max_x_point()
        dummy_last_point = Coord(start_point.x - 1.0, start_point.y)
        self.find_envelope_points(start_point, dummy_last_point)
        self.envelope = CoordCollection(self.envelope)

    def plot_envelope(self, quadrant=None, **kwargs):
        if not self.envelope:
            self.get_envelope()

        if quadrant:
            plot_list = []
            coords_zip = list(zip(self.envelope.coords[:-1], self.envelope.coords[1:]))
            for points in coords_zip:
                quadrants = points[0].convex_quadrant(points[1])
                if quadrant in quadrants:
                    plot_list.append(points[0])
                    plot_list.append(points[1])

            plot_coords = CoordCollection(plot_list)
        else:
            plot_coords = self.envelope

        for i, x in enumerate(plot_coords.x_list):
            plt.scatter(x, plot_coords.y_list[i], **kwargs)
            plt.text(x, plot_coords.y_list[i], plot_coords.name_list[i])

        plt.plot(plot_coords.x_list, plot_coords.y_list, **kwargs)


class Line(object):
    """Line in the form y = mx + b:

    Attributes:
        gradient (float): Gradient of line
        y_intercept (float): Value at which line crosses the y-axis
    """

    def __init__(self, gradient, y_intercept):
        self.gradient = gradient
        self.y_intercept = y_intercept

    def find_y_at_x(self, x):
        """ Finds the value of y for a given value of x

        Args:
            x (float): Nominal x value
        Returns:
            float: Value of y at nominal x value
        """
        # y = mx + b
        return x * self.gradient + self.y_intercept

    def find_x_at_y(self, y):
        """ Finds the value of x for a given value of y

        Args:
            y (float): Nominal y value
        Returns:
            float: Value of x at nominal x value
        """
        # rearange y = mx + b to find x
        return (y - self.y_intercept) / self.gradient

    def find_intercept_on_line(self, other_line):
        ''' Finds the coordinates of intecept of self and another Line
            (if one exists)

        Args:
            other_line (Line): Other Line object which may intersect with self

        Returns:
            Coord: Coord object representing intercept x and y values
        '''

        m_self = self.gradient
        m_other = other_line.gradient

        # Test for parallel lines
        if m_self != m_other:
            # Find intercept two lines in the form y = mx + b
            # Just high school algebra...
            b_self = self.y_intercept
            b_other = other_line.y_intercept
            x = (b_self - b_other) / (m_other - m_self)
            intercept = Coord(x, self.find_y_at_x(x))
        else:
            intercept = Coord(None, None)  # no intercept for parallel lines
        return intercept

    def set_plot(self, x_bounds, **kwargs):

        x_vals = x_bounds

        y_val_1 = self.find_y_at_x(x_bounds[0])
        y_val_2 = self.find_y_at_x(x_bounds[1])
        y_vals = [y_val_1, y_val_2]

        return XYline(x_vals, y_vals, label=self.equation_string, **kwargs)

