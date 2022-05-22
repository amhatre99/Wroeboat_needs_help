import Triangle
import Point
import math
import random

class Rectangle_sample:
    def __init__(self, tri: Triangle):
        # m1 is the slope of the line formed by points a and c.
        # m1 is also equal to tan(theta) where theta is angle cab.
        self.m1 = Point.slope(tri.a, tri.c)
        # m2 is the slope of the line formed by points b and c.
        # m2 is also equal to tan(phi) where phi is the angle formed by line ab and the x-axis.
        self.m2 = Point.slope(tri.b, tri.c)
        # b is the y intercept of the line formed by points b and c.
        self.x1 = tri.b.x
        self.x2 = tri.c.x
        self.y2 = tri.c.y
        # uniformly distributed points in triangle generated so far
        self.points = []

    def generate_samples(self, n: int = 1, seed: int = 11):
        # using starting seed for random, so we can have the same dataset if desired.
        random.seed(seed)
        for i in range(n):
            # sample random points x and y
            randx = random.random()
            randy = random.random()
            # reflect points so that all samples lie on triangle
            randx, randy = max(randx, randy), min(randx, randy)

            # xtransform = self.x1 * (2 / (self.m1 ** 2 + 1) - 1)
            # ytransform = self.x1 * 2 * self.m1 / (self.m1 ** 2 + 1)

            self.points.append([randx * self.x1 + randy * (self.x2 - self.x1), randy * self.y2])
