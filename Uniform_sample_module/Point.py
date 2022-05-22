from math import sqrt

class Point:
    def __init__(self, x, y):
        """
        Constructor for point with x coord. and y coord.
        :param x: x coordinate of point
        :param y: y coordinate of point
        """
        assert isinstance(x, (int, float)) and isinstance(y, (int, float)), "Inputs x and y must be integers or floats."
        self.x = x
        self.y = y

    def copy(self):
        return Point(self.x, self.y)

    def to_list(self):
        return [self.x, self.y, 0]

def distance(p1: Point, p2: Point):
    """
    Distance between two points.
    :param p1: point 1
    :param p2: point 2
    """
    return sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2)

def slope(p1: Point, p2: Point):
    """
    Slope between two points.
    :param p1: point 1
    :param p2: point 2
    """
    return (p1.y - p2.y) / (p1.x - p2.x)

def y_intercept(p1: Point, p2: Point):
    """
    y intercept of the line between two points.
    :param p1: point 1
    :param p2: point 2
    """
    return  p1.y - p1.x * slope(p1, p2)

def line(p1: Point, p2: Point):
    """
    Function for line between two points.
    :param p1: point 1
    :param p2: point 2
    """
    return lambda x: slope(p1, p2) * x + y_intercept(p1, p2)

def midpoint(p1: Point, p2: Point):
    """
    Midpoint between two points
    :param p1: point 1
    :param p2: point 2
    """
    return Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2)
