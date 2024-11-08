#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <array>


class Line;
class Polygon;
class Ellipse;

namespace {
    const double eps = 1e-9;
    const double TO_RADIANS = M_PI / 180;

    bool equal(double x, double y) {
        return std::abs(x - y) < eps;
    }

    int sign(double x) {
        return std::abs(x) < eps ? 0 : 2 * (x > 0) - 1;
    }

    template<typename T>
    using iter = typename std::vector<T>::const_iterator;

    template<typename T>
    void cycle_next(iter<T> &it, const std::vector<T> &a) {
        ++it;
        if (it == a.end()) {
            it = a.begin();
        }
    }

    template<typename T>
    bool cycle_equal(const std::vector<T> &a, const std::vector<T> &b) {
        for (auto it = a.begin(); it != a.end(); ++it) {
            if (*it == b[0]) {
                auto t1 = it;
                cycle_next(t1, a);
                for (auto t2 = b.end() - 1; t2 != b.begin() && *t1 == *t2; --t2) {
                    cycle_next(t1, a);
                }
                if (t1 == it) {
                    return true;
                }

                t1 = it;
                cycle_next(t1, a);
                for (auto t2 = b.begin() + 1; t2 != b.end() && *t1 == *t2; ++t2) {
                    cycle_next(t1, a);
                }
                if (t1 == it) {
                    return true;
                }
            }
        }
        return false;
    }
}

struct Point {
    double x;
    double y;

    Point(double x, double y) : x(x), y(y) {}

    void rotate(const Point& center, double a, double b) {
        Point t{x - center.x, y - center.y};
        x = a * t.x - b * t.y + center.x;
        y = a * t.y + b * t.x + center.y;
    }
    void reflect(const Line& axis);
    void scale(const Point& center, double coefficient) {
        x = coefficient * (x - center.x) + center.x;
        y = coefficient * (y - center.y) + center.y;
    }
};

Point middle(const Point &A, const Point &B) {
    return Point{(A.x + B.x) / 2, (A.y + B.y) / 2};
}

bool operator==(const Point &A, const Point &B) {
    return equal(A.x, B.x) && equal(A.y, B.y);
}

bool operator!=(const Point &A, const Point &B) {
    return !(A == B);
}

double scalar_product(const Point &a, const Point &b) {
    return a.x * b.x + a.y * b.y;
}

double scalar_product(const Point &A, const Point &B, const Point &C, const Point &D) {
    return (B.x - A.x) * (D.x - C.x) + (B.y - A.y) * (D.y - C.y);
}

double vector_product(const Point &a, const Point &b) {
    return a.x * b.y - b.x * a.y;
}

double vector_product(const Point &A, const Point &B, const Point &C, const Point &D) {
    return (B.x - A.x) * (D.y - C.y) - (D.x - C.x) * (B.y - A.y);
}

double dist2(const Point &A, const Point &B) {
    return (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y);
}

double dist(const Point &A, const Point &B) {
    return sqrt(dist2(A, B));
}
class Line {
public:
    double a;
    double b;
    double c;

    Line(const Point& A, const Point& B): a(B.y - A.y), b(A.x - B.x), c(-(a * A.x + b * A.y)) {}
    Line(double k, double b): a(k), b(-1), c(b) {}
    Line(const Point& A, double k): a(k), b(-1), c(-(k * A.x + b * A.y)) {}

    double operator()(const Point& A) const {
        return a * A.x + b * A.y + c;
    }
    double operator()(double x, double y) const {
        return a * x + b * y + c;
    }
};

void Point::reflect(const Line &axis) {
    double k = 2 * axis(*this) / (axis.a * axis.a + axis.b * axis.b);
    x -= k * axis.a;
    y -= k * axis.b;
}

namespace {
    bool isIntersect(const Point &point, const Point &A, const Point &B) {
        auto [y1, y2] = std::minmax(A.y, B.y);
        return y1 - eps < point.y && point.y < y2 - 2 * eps &&
               point.x < A.x + (B.x - A.x) / (B.y - A.y) * (point.y - A.y) + eps;
    }
}

bool segmentContainsPoint(const Point &A, const Point &B, const Point &point) {
    return equal(dist(point, A) + dist(point, B), dist(A, B));
}

bool operator==(const Line& l, const Line& k) {
    return equal(l.a * k.b, l.b * k.a) && equal(l.b * k.c, l.c * k.b) && equal(l.a * k.c, l.c * k.a);
}
bool operator!=(const Line& l, const Line& k) {
    return !(l == k);
}

class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;

    virtual bool isEqualTo(const Shape& another) const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;

    virtual bool containsPoint(const Point& point) const = 0;

    virtual void rotate(const Point& center, double angle) = 0;
    virtual void reflect(const Point& center) = 0;
    virtual void reflect(const Line& axis) = 0;
    virtual void scale(const Point& center, double coefficient) = 0;

    virtual ~Shape() {};
};

class Ellipse : public Shape {
protected:
    Point F1;
    Point F2;
    double half_big_axis;
public:
    Ellipse(const Point& F1, const Point& F2, double a_2) : F1(F1), F2(F2), half_big_axis(a_2 / 2) {}
    double get_half_big_axis() const {
        return half_big_axis;
    }
    double get_half_small_axis() const {
        //NOLINTNEXTLINE(readability-magic-numbers)
        return sqrt(half_big_axis * half_big_axis - dist2(F1, F2) / 4);
    }
    std::pair<Point,Point> focuses() const {
        return {F1, F2};
    }
    Point center() const {
        return middle(F1, F2);
    }
    std::pair<Line, Line> directrices() const {
        Point O = center();
        double dx = (F2.x - F1.x) / dist(F1, F2) * (half_big_axis / eccentricity());
        double dy = (F2.y - F1.y) / dist(F1, F2) * (half_big_axis / eccentricity());

        Point A{O.x + dx, O.y + dy};
        Line D1(A, Point{A.x + (F2.y - F1.y), A.y + (F1.x - F2.x)});

        Point B{O.x - dx, O.y - dy};
        Line D2(B, Point{B.x + (F2.y - F1.y), B.y + (F1.x - F2.x)});
        return {D1, D2};
    }
    double eccentricity() const {
        return dist(F1, F2) / (2 * half_big_axis);
    }
    double perimeter() const override {
        double b = get_half_small_axis();
        //NOLINTNEXTLINE(readability-magic-numbers)
        return 4 * half_big_axis * std::comp_ellint_2(sqrt(1 - b * b / half_big_axis / half_big_axis));
    }
    double area() const override {
        return M_PI * half_big_axis * get_half_small_axis();
    }

    bool isEqualTo(const Shape& another) const override {
        const Ellipse* ptr = dynamic_cast<const Ellipse*>(&another);
        return ptr && equal(half_big_axis, ptr->half_big_axis) && ((F1 == ptr->F1 && F2 == ptr->F2) || (F1 == ptr->F2 && F2 == ptr->F1));
    }

    bool isCongruentTo(const Shape& another) const override {
        const Ellipse* ptr = dynamic_cast<const Ellipse*>(&another);
        return ptr && equal(half_big_axis, ptr->half_big_axis) && equal(dist2(F1, F2), dist2(ptr->F1, ptr->F2));
    }

    bool isSimilarTo(const Shape& another) const override {
        const Ellipse* ptr = dynamic_cast<const Ellipse*>(&another);
        return ptr && equal(eccentricity(), ptr->eccentricity());
    }

    bool containsPoint(const Point& point) const override {
        return dist(F1, point) + dist(F2, point) < 2 * half_big_axis + eps;
    }


    void rotate(const Point& center, double angle) override {
        angle *= TO_RADIANS;
        double x = cos(angle);
        double y = sin(angle);
        F1.rotate(center, x, y);
        F2.rotate(center, x, y);
    }
    void reflect(const Point& center) override {
        F1.rotate(center, -1, 0);
        F2.rotate(center, -1, 0);
    }
    void reflect(const Line& axis) override {
        F1.reflect(axis);
        F2.reflect(axis);
    }
    void scale(const Point& center, double coefficient) override {
        F1.scale(center, coefficient);
        F2.scale(center, coefficient);
        half_big_axis *= std::abs(coefficient);
    }
};

class Circle : public Ellipse {
public:
    Circle(const Point& center, double r) : Ellipse(center, center, 2 * r) {}

    double radius() const {
        return half_big_axis;
    }
    double perimeter() const override {
        return 2 * M_PI * half_big_axis;
    }
    double area() const override {
        return M_PI * half_big_axis * half_big_axis;
    }
    bool containsPoint(const Point& point) const override {
        return dist2(F1, point) < half_big_axis * half_big_axis + eps;
    }
};

class Polygon : public Shape {
private:

    struct Congruent {
        double scal = 0;
        double vect = 0;
        double ratio = 0;

        Congruent() {}
        Congruent(Point A, Point B, Point C) :
                scal(scalar_product(A, B, B, C)),
                vect(vector_product(A, B, B, C)),
                ratio(dist2(B, C) / dist2(A, B)) {}

        void reflect() {
            vect *= -1;
        }
        bool operator==(const Congruent& x) const {
            return equal(scal, x.scal) && equal(vect, x.vect) && equal(ratio, x.ratio);
        }
    };

    struct Similar {
        double cosinus = 0;
        double sinus = 0;
        double ratio = 0;

        Similar() {}
        Similar(Point A, Point B, Point C) :
                cosinus(scalar_product(A, B, B, C) / (dist(A, B) * dist(B, C))),
                sinus(vector_product(A, B, B, C) / (dist(A, B) * dist(B, C))),
                ratio(dist(B, C) / dist(A, B)) {}

        void reflect() {
            sinus *= -1;
        }
        bool operator==(const Similar& x) const {
            return equal(cosinus, x.cosinus) && equal(sinus, x.sinus) && equal(ratio, x.ratio);
        }
    };

    template <typename Angle>
    bool checkEquality(const Polygon& another) const {
        if (is_convex != another.is_convex || vertices.size() != another.vertices.size()) {
            return false;
        }

        std::vector<Angle> a(vertices.size());
        std::vector<Angle> b(vertices.size());
        a[0] = Angle(vertices.back(), vertices[0], vertices[1]);
        b[0] = Angle(another.vertices.back(), another.vertices[0], another.vertices[1]);
        a.back() = Angle(vertices[vertices.size() - 2], vertices.back(), vertices[0]);
        b.back() = Angle(another.vertices[vertices.size() - 2], another.vertices.back(), another.vertices[0]);

        for (size_t i = 1; i + 1 < a.size(); ++i) {
            a[i] = Angle(vertices[i - 1], vertices[i], vertices[i + 1]);
            b[i] = Angle(another.vertices[i - 1], another.vertices[i], another.vertices[i + 1]);
        }
        if (cycle_equal(a, b)) {
            return true;
        }
        for (size_t i = 0; i < a.size(); ++i) {
            a[i].reflect();
        }
        if (cycle_equal(a, b)) {
            return true;
        }
        for (size_t i = 0; i < a.size(); ++i) {
            a[i].ratio = 1 / a[i].ratio;
        }
        if (cycle_equal(a, b)) {
            return true;
        }
        for (size_t i = 0; i < a.size(); ++i) {
            a[i].reflect();
        }
        return cycle_equal(a, b);
    }

protected:
    std::vector<Point> vertices;
    bool is_convex = true;

public:
    Polygon() {};

    Polygon(const std::vector<Point>& points) : vertices(points) {
        auto it = vertices.begin();
        const int s = sign(vector_product(*it, vertices.back(), *it, *(it + 1)));
        ++it;
        while (it + 1 < vertices.end() && sign(vector_product(*it, *(it - 1), *it, *(it + 1))) == s) {
            ++it;
        }
        is_convex = it + 1 == vertices.end() && sign(vector_product(*it, *(it - 1), *it, vertices[0])) == s;
    }

    template <typename... Ts>
    Polygon(const Ts&... args) : Polygon(std::vector<Point>({args...})) {}

    Polygon(std::initializer_list<Point> args) : Polygon(std::vector<Point>(args)) {}

    size_t verticesCount() const {
        return vertices.size();
    }
    const std::vector<Point>& getVertices() const {
        return vertices;
    }
    bool isConvex() const {
        return is_convex;
    }

    double perimeter() const override {
        double p = dist(vertices[0], vertices.back());
        for (size_t i = 1; i < vertices.size(); ++i) {
            p += dist(vertices[i - 1], vertices[i]);
        }
        return p;
    }
    double area() const override {
        double ans = 0;
        for (size_t i = 1; i + 1 < vertices.size(); ++i) {
            ans += vector_product(vertices[0], vertices[i], vertices[0], vertices[i + 1]);
        }
        return std::abs(ans) / 2;
    }

    bool isEqualTo(const Shape& another) const override {
        const Polygon* ptr = dynamic_cast<const Polygon*>(&another);
        return (ptr && is_convex == ptr->is_convex && vertices.size() == ptr->vertices.size() &&
                cycle_equal(vertices, ptr->vertices));
    }

    bool isCongruentTo(const Shape& another) const override {
        const Polygon* ptr = dynamic_cast<const Polygon*>(&another);
        return ptr && checkEquality<Congruent>(*ptr);
    }

    bool isSimilarTo(const Shape& another) const override {
        const Polygon* ptr = dynamic_cast<const Polygon*>(&another);
        return ptr && checkEquality<Similar>(*ptr);
    }

    bool containsPoint(const Point& point) const override {
        int ans = isIntersect(point, vertices[0], vertices.back());
        if (segmentContainsPoint(vertices[0], vertices.back(), point)) {
            return true;
        }
        for (size_t i = 1; i < vertices.size(); ++i) {
            if (segmentContainsPoint(vertices[i - 1], vertices[i], point)) {
                return true;
            }
            ans += isIntersect(point, vertices[i - 1], vertices[i]);
        }
        return ans % 2;
    }

    void rotate(const Point& center, double angle) override {
        angle *= TO_RADIANS;
        double x = cos(angle);
        double y = sin(angle);
        for (auto& point : vertices) {
            point.rotate(center, x, y);
        }
    }
    void reflect(const Point& center) override {
        for (auto& point : vertices) {
            point.rotate(center, -1, 0);
        }
    }
    void reflect(const Line& axis) override {
        for (auto& point : vertices) {
            point.reflect(axis);
        }
    }
    void scale(const Point& center, double coefficient) override {
        for (auto& point : vertices) {
            point.scale(center, coefficient);
        }
    }
};

class Rectangle : public Polygon {
public:
    Rectangle(const Point& A, const Point& C, double k) {
        Point B = C;
        if (k < 1) {
            k = 1 / k;
        }
        double angle = atan(k);
        k = sqrt(1 + k * k);
        B.rotate(A, cos(angle) / k, sin(angle) / k);
        Point D = B;
        D.rotate(middle(A, C), -1, 0);

        vertices = {A, B, C, D};
    }
    Point center() const {
        return middle(vertices[0], vertices[1]);
    }
    std::pair<Line, Line> diagonals() const {
        //NOLINTNEXTLINE(readability-magic-numbers)
        return {Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3])};
    }
};

class Square : public Rectangle {
public:
    Square(const Point& A, const Point& B) : Rectangle(A, B, 1) {}
    Circle circumscribedCircle() const {
        return Circle(center(), dist(vertices[0], vertices[2]) / 2);
    }
    Circle inscribedCircle() const {
        return Circle(center(), dist(vertices[0], vertices[1]) / 2);
    }
};

namespace {
    double orthocenter_formula(double x1, double y1, double x2, double y2, double x3, double y3) {
        return ((y1 * y1 + x2 * x3) * (y3 - y2) + (y2 * y2 + x1 * x3) * (y1 - y3) + (y3 * y3 + x1 * x2) * (y2 - y1)) /
               (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
    }
}

class Triangle : public Polygon {
private:
    //NOLINTNEXTLINE(readability-magic-numbers)
    std::array<double, 6> coordinates() const {
        return {vertices[0].x, vertices[0].y, vertices[1].x, vertices[1].y, vertices[2].x, vertices[2].y};
    }
public:
    Triangle(const Point& A, const Point& B, const Point& C) : Polygon({A, B, C}) {}
    Point circumscribedCenter() const {
        auto [x1, y1, x2, y2, x3, y3] = coordinates();

        double z1 = x1 * x1 + y1 * y1;
        double z2 = x2 * x2 + y2 * y2;
        double z3 = x3 * x3 + y3 * y3;

        double zx = (y2 - y3) * z1 + (y3 - y1) * z2 + (y1 - y2) * z3;
        double zy = (x2 - x3) * z1 + (x3 - x1) * z2 + (x1 - x2) * z3;
        double z = 2 * ((x1 - x2) * (y3 - y1) - (y1 - y2) * (x3 - x1));
        return Point{-zx / z, zy / z};
    }
    Circle circumscribedCircle() const {
        Point O = circumscribedCenter();
        return Circle(O, dist(O, vertices[0]));
    }
    Circle inscribedCircle() const {
        auto [x1, y1, x2, y2, x3, y3] = coordinates();

        double a = dist(vertices[1], vertices[2]);
        double b = dist(vertices[2], vertices[0]);
        double c = dist(vertices[0], vertices[1]);
        double p = a + b + c;

        Point I{(a * x1 + b * x2 + c * x3) / p, (a * y1 + b * y2 + c * y3) / p};
        p /= 2;
        return Circle(I, sqrt((p - a) * (p - b) * (p - c) / p));

    }
    Point centroid() const {
        auto [x1, y1, x2, y2, x3, y3] = coordinates();
        //NOLINTNEXTLINE(readability-magic-numbers)
        return Point{(x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3};
    }
    Point orthocenter() const {
        auto [x1, y1, x2, y2, x3, y3] = coordinates();
        return Point{orthocenter_formula(x1, y1, x2, y2, x3, y3), orthocenter_formula(y1, x1, y2, x2, y3, x3)};
    }
    Line EulerLine() const {
        return Line(circumscribedCenter(), centroid());
    }
    Circle ninePointsCircle() const {
        Point H = orthocenter();
        Point O = circumscribedCenter();
        return Circle(middle(O, H), dist(O, vertices[0]) / 2);
    }
};

bool operator==(const Shape& a, const Shape& b) {
    return a.isEqualTo(b);
}
