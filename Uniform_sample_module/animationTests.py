from manim import *
from Point import *
from Triangle import *
from CDF_inverse_method import *
from Rectangle_sampling_method import *
from math import *

class testTriangle(Scene):
    def construct(self):
        grid = NumberPlane(x_range=(-10, 10, 1), y_range=(-6.0, 6.0, 1))
        self.add(grid)
        p1 = Point.Point(2, 3)
        p2 = Point.Point(4, -1)
        p3 = Point.Point(1, 1)
        t = Triangle.Triangle(p1, p2, p3)

        a = Line(p1.to_list(), p2.to_list())
        a.append_points(Line(p2.to_list(), p3.to_list()).points)
        a.append_points(Line(p3.to_list(), p1.to_list()).points)
        b = Line(t.a.to_list(), t.b.to_list())
        b.append_points(Line(t.b.to_list(), t.c.to_list()).points)
        b.append_points(Line(t.c.to_list(), t.a.to_list()).points)
        d = t.denorm(t.a).to_list()
        e = t.denorm(t.b).to_list()
        f = t.denorm(t.c).to_list()
        c = Line(d, e)
        c.append_points(Line(e, f).points)
        c.append_points(Line(f, d).points)
        self.play(Create(grid, run_time=2, lag_ratio=0.1))
        self.play(Create(a))
        self.wait()
        self.play(Transform(a, b))
        self.wait()
        self.play(Transform(a, c))

class testCDFinverse(Scene):
    def construct(self):
        grid = NumberPlane(x_range=(-.1, 1.9, .1), y_range=(-0.1, 1.1, .1), x_length = 20, y_length = 12)
        p1 = Point.Point(2, 3)
        p2 = Point.Point(4, -1)
        p3 = Point.Point(1, 1)
        t = Triangle.Triangle(p1, p2, p3)

        a = Line(p1.to_list(), p2.to_list())
        a.append_points(Line(p2.to_list(), p3.to_list()).points)
        a.append_points(Line(p3.to_list(), p1.to_list()).points)
        b = Line(t.a.to_list(), t.b.to_list())
        b.append_points(Line(t.b.to_list(), t.c.to_list()).points)
        b.append_points(Line(t.c.to_list(), t.a.to_list()).points)
        # self.add(b)
        self.play(Create(grid, run_time=2, lag_ratio=0.1))
        self.play(Create(a))
        self.wait()
        self.play(Transform(a, b))
        self.wait()

        CDF_generator = CDF_Inverse(t)
        CDF_generator.generate_samples(10000)
        x_y_points = [[r * cos(theta), r * sin(theta), 0] for (r, theta) in CDF_generator.points]
        centroid_sum = [0, 0]
        for p in x_y_points:
            d = Dot(p, radius=0.005)
            d.scale(5)
            #self.add(d)
            self.play(Create(d), run_time=4/len(x_y_points))
            centroid_sum[0] += p[0]
            centroid_sum[1] += p[1]

        estimated_centroid = Point.Point(centroid_sum[0] / len(x_y_points), centroid_sum[1] / len(x_y_points))
        print("estimated center", t.denorm(estimated_centroid).to_list())
        print("actual center", t.centroid().to_list())
        print("distance from estimated center to actual center", Point.distance(t.denorm(estimated_centroid), t.centroid()))


class testRectangle_sampling(Scene):
    def construct(self):
        grid = NumberPlane(x_range=(-.1, 1.9, .1), y_range=(-0.1, 1.1, .1), x_length = 20, y_length = 12)
        p1 = Point.Point(2, 3)
        p2 = Point.Point(4, -1)
        p3 = Point.Point(1, 1)
        t = Triangle.Triangle(p1, p2, p3)

        a = Line(p1.to_list(), p2.to_list())
        a.append_points(Line(p2.to_list(), p3.to_list()).points)
        a.append_points(Line(p3.to_list(), p1.to_list()).points)
        b = Line(t.a.to_list(), t.b.to_list())
        b.append_points(Line(t.b.to_list(), t.c.to_list()).points)
        b.append_points(Line(t.c.to_list(), t.a.to_list()).points)
        #self.add(b)
        self.play(Create(grid, run_time=2, lag_ratio=0.1))
        self.play(Create(a))
        self.wait()
        self.play(Transform(a, b))
        self.wait()

        rectangle_samp_gen = Rectangle_sample(t)
        rectangle_samp_gen.generate_samples(1000)
        x_y_points = rectangle_samp_gen.points
        centroid_sum = [0, 0]
        for p in x_y_points:
            d = Dot(p + [0], radius=0.005)
            d.scale(5)
            #self.add(d)
            self.play(Create(d), run_time=4/len(x_y_points))
            centroid_sum[0] += p[0]
            centroid_sum[1] += p[1]

        estimated_centroid = Point.Point(centroid_sum[0] / len(x_y_points), centroid_sum[1] / len(x_y_points))
        print("estimated center", t.denorm(estimated_centroid).to_list())
        print("actual center", t.centroid().to_list())
        print("distance from estimated center to actual center", Point.distance(t.denorm(estimated_centroid), t.centroid()))
