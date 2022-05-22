import Triangle
import Point
import math
import random

class CDF_Inverse:
    def __init__(self, tri: Triangle):
        # m1 is the slope of the line formed by points a and c.
        # m1 is also equal to tan(theta) where theta is angle cab.
        self.m1 = Point.slope(tri.a, tri.c)
        # m2 is the slope of the line formed by points b and c.
        # m2 is also equal to tan(phi) where phi is the angle formed by line ab and the x-axis.
        self.m2 = Point.slope(tri.b, tri.c)
        # b is the y intercept of the line formed by points b and c.
        self.b = Point.y_intercept(tri.b, tri.c)
        # uniformly distributed points in triangle generated so far
        self.points = []

    def generate_samples(self, n: int = 1, seed: int = 11):
        # using starting seed for random, so we can have the same dataset if desired.
        random.seed(seed)
        for i in range(n):
            monte_theta = random.random()
            monte_r = random.random()
            theta = math.atan(monte_theta * self.m1 * self.m2 / (self.m2 + self.m1 * (monte_theta - 1)))
            r = math.sqrt(monte_r) * self.b / (math.sin(theta) - self.m2 * math.cos(theta))
            self.points.append([r, theta])
