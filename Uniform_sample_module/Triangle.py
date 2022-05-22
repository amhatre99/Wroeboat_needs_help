import Point

class Triangle:
    def __init__(self, p1: Point, p2: Point, p3: Point):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

        # Creating transformation function to normalize triangle
        # A, B, C are side lengths of triangle
        # Find point to be mapped to origin and point to check reflection condition
        A = Point.distance(p2, p3)
        B = Point.distance(p1, p3)
        C = Point.distance(p1, p2)
        origin_point = p2
        reflect_point = p1

        # Check reflection condition and which point to be mapped to origin
        if (A - B) * (C ** 2 - (A + B) ** 2) > 0:
            origin_point = p1
            reflect_point = p2

        # Create rotation tuple: (sin(theta), cos(theta))
        rotate_by = ((p3.y - origin_point.y) / Point.distance(p3, origin_point),
                     (p3.x - origin_point.x) / Point.distance(p3, origin_point))
        # Create rotation function
        rotate_clockwise = lambda x, y: (rotate_by[1] * x + rotate_by[0] * y, -rotate_by[0] * x + rotate_by[1] * y)
        rotate_counter_clockwise = lambda x, y: (rotate_by[1] * x - rotate_by[0] * y, rotate_by[0] * x + rotate_by[1] * y)

        # Set reflection constant
        reflect = 1
        if rotate_clockwise(reflect_point.x - origin_point.x, reflect_point.y - origin_point.y)[1] < 0:
            reflect = -1

        # Set scaling constant
        scale = rotate_clockwise(p3.x - origin_point.x, p3.y - origin_point.y)[0]

        # Define normalization transformation function and set it as an instance variable
        # Order of transformation: translate, rotate, scale & reflect if needed
        def norm(p: Point):
            newx = p.x
            newy = p.y
            newx -= origin_point.x
            newy -= origin_point.y
            newx, newy = rotate_clockwise(newx, newy)
            newx /= (scale / 5)
            newy /= (scale / 5) * reflect
            return Point.Point(newx, newy)

        self.norm = norm

        # Define denormaliztion transformation function and set it as an instance variable
        # Takes any point in the normalized triangle and maps it back to the corresponding point
        # Order of transformation: scale & reflect, rotate, translate  if needed
        def denorm(p: Point):
            newx = p.x
            newy = p.y
            newx *= scale / 5
            newy *= scale / 5 * reflect
            newx, newy = rotate_counter_clockwise(newx, newy)
            newx += origin_point.x
            newy += origin_point.y
            return Point.Point(newx, newy)

        self.denorm = denorm

        self.a = self.norm(origin_point)
        self.b = self.norm(self.p3)
        self.c = self.norm(reflect_point)


    def centroid(self):
        """
        Centroid of triangle object.
        """
        return Point.Point((self.p1.x + self.p2.x + self.p3.x) / 3, (self.p1.y + self.p2.y + self.p3.y) / 3)


