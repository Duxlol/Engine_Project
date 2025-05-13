#include <algorithm>
#include "easy_image.h"
#include "ini_configuration.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <stack>
#include "vector3d.h"
#include "l_parser.h"

// 2D structures
struct Point2D {
    double x, y;
};

struct Line2D {
    Point2D start, end;
    std::vector<double> color;
    double z1 = 0.0; // Z-coordinate of start point (for Z-buffering)
    double z2 = 0.0; // Z-coordinate of end point (for Z-buffering)
};

typedef std::vector<Line2D> Lines2D;

// 3D figure structures
struct Face {
    std::vector<int> point_indices;
};

struct Figure {
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    std::vector<std::pair<int, int>> lines;
    std::vector<double> color;
};

typedef std::vector<Figure> Figures3D;

// Class to represent a simple color
struct Color {
    double red, green, blue;
};

// ZBuffer class for Z-buffering
class ZBuffer {
private:
    std::vector<std::vector<double>> buffer;
    unsigned int width, height;

public:
    ZBuffer(unsigned int width, unsigned int height)
        : width(width), height(height), buffer(height, std::vector<double>(width, std::numeric_limits<double>::infinity())) {}

    void set(unsigned int x, unsigned int y, double z) {
        if (x < width && y < height && z < buffer[y][x]) {
            buffer[y][x] = z;
        }
    }

    double get(unsigned int x, unsigned int y) const {
        if (x < width && y < height) {
            return buffer[y][x];
        }
        return std::numeric_limits<double>::infinity();
    }
};

// Transformation matrices
Matrix scalingMatrix(double scale) {
    Matrix S;
    S(1, 1) = scale;
    S(2, 2) = scale;
    S(3, 3) = scale;
    return S;
}

Matrix rotationX(double angle) {
    // Convert angle from degrees to radians
    double radians = angle * M_PI / 180.0;

    Matrix R;
    R(1,1) = 1;
    R(2, 2) = cos(radians);
    R(2, 3) = sin(radians);
    R(3, 2) = -sin(radians);
    R(3, 3) = cos(radians);

    return R;
}

Matrix rotationY(double angle) {
    // Convert angle from degrees to radians
    double radians = angle * M_PI / 180.0;

    Matrix R;
    R(1, 1) = cos(radians);
    R(1, 3) = -sin(radians);
    R(3, 1) = sin(radians);
    R(3, 3) = cos(radians);

    return R;
}

Matrix rotationZ(double angle) {
    // Convert angle from degrees to radians
    double radians = angle * M_PI / 180.0;

    Matrix R;
    R(1, 1) = cos(radians);
    R(1, 2) = sin(radians);
    R(2, 1) = -sin(radians);
    R(2, 2) = cos(radians);

    return R;
}

Matrix translationMatrix(const Vector3D& translation) {
    Matrix T;
    T(4, 1) = translation.x;
    T(4, 2) = translation.y;
    T(4, 3) = translation.z;

    return T;
}

Matrix eyePointTrans(const Vector3D& eye) {
    double r = sqrt(eye.x * eye.x + eye.y * eye.y + eye.z * eye.z);
    double theta = atan2(eye.y, eye.x);
    double phi = acos(eye.z / r);

    Matrix V;
    V(1, 1) = -sin(theta);
    V(1, 2) = -cos(theta) * cos(phi);
    V(1, 3) = cos(theta) * sin(phi);
    V(2, 1) = cos(theta);
    V(2, 2) = -sin(theta) * cos(phi);
    V(2, 3) = sin(theta) * sin(phi);
    V(3, 2) = sin(phi);
    V(3, 3) = cos(phi);
    V(4, 3) = -r;

    return V;
}

Point2D doProjection(const Vector3D& point, double d = 1.0) {
    double x = -(d * point.x) / point.z;
    double y = -(d * point.y) / point.z;
    return {x, y};
}

Lines2D doProjection(const Figures3D& figures, double d = 1.0) {
    Lines2D lines;

    for (const auto& figure : figures) {
        // project lines from figure.lines
        for (const auto& line : figure.lines) {
            const Vector3D& p1_3d = figure.points[line.first];
            const Vector3D& p2_3d = figure.points[line.second];

            // skip lines behind eye
            if (p1_3d.z >= 0 || p2_3d.z >= 0) {
                continue;
            }

            // project points to 2D n preserve Z values for z-buffering
            Point2D p1 = doProjection(p1_3d, d);
            Point2D p2 = doProjection(p2_3d, d);

            Line2D line2d;
            line2d.start = p1;
            line2d.end = p2;
            line2d.color = figure.color;
            line2d.z1 = p1_3d.z;
            line2d.z2 = p2_3d.z;

            lines.push_back(line2d);
        }

        // project faces for filled polygons
        for (const auto& face : figure.faces) {
            if (face.point_indices.size() >= 2) {
                bool allPointsVisible = true;

                // check if any point in the face is behind the eye
                for (size_t i = 0; i < face.point_indices.size(); i++) {
                    if (figure.points[face.point_indices[i]].z >= 0) {
                        allPointsVisible = false;
                        break;
                    }
                }

                if (!allPointsVisible) {
                    continue;
                }

                // convert face to lines
                for (size_t i = 0; i < face.point_indices.size(); i++) {
                    size_t j = (i + 1) % face.point_indices.size();
                    const Vector3D& p1_3d = figure.points[face.point_indices[i]];
                    const Vector3D& p2_3d = figure.points[face.point_indices[j]];

                    Point2D p1 = doProjection(p1_3d, d);
                    Point2D p2 = doProjection(p2_3d, d);

                    Line2D line;
                    line.start = p1;
                    line.end = p2;
                    line.color = figure.color;
                    line.z1 = p1_3d.z;
                    line.z2 = p2_3d.z;

                    lines.push_back(line);
                }
            }
        }
    }

    return lines;
}

std::vector<int> scaleColor(const std::vector<double> &color) {
    std::vector<int> scaled(3);
    for (int i = 0; i < 3; ++i) {
        scaled[i] = static_cast<int>(std::round(color[i] * 255));
        scaled[i] = std::clamp(scaled[i], 0, 255);
    }
    return scaled;
}

// 3D shapes
Figure createCube(const std::vector<double>& color) {
    Figure cube;
    cube.color = color;

    // points of the cube
    cube.points = {
        Vector3D::point(1, -1, -1),  // 0
        Vector3D::point(-1, 1, -1),  // 1
        Vector3D::point(1, 1, 1),    // 2
        Vector3D::point(-1, -1, 1),  // 3
        Vector3D::point(1, 1, -1),   // 4
        Vector3D::point(-1, -1, -1), // 5
        Vector3D::point(1, -1, 1),   // 6
        Vector3D::point(-1, 1, 1)    // 7
    };

    // define faces with point indices
    Face f1, f2, f3, f4, f5, f6;
    f1.point_indices = {0, 4, 2, 6};
    f2.point_indices = {4, 1, 7, 2};
    f3.point_indices = {1, 5, 3, 7};
    f4.point_indices = {5, 0, 6, 3};
    f5.point_indices = {6, 2, 7, 3};
    f6.point_indices = {0, 5, 1, 4};

    cube.faces = {f1, f2, f3, f4, f5, f6};

    // define lines for wireframe drawing
    cube.lines = {
        {0, 4}, {4, 1}, {1, 5}, {5, 0},
        {2, 6}, {6, 3}, {3, 7}, {7, 2},
        {0, 6}, {4, 2}, {1, 7}, {5, 3}
    };

    return cube;
}

Figure createTetra(const std::vector<double>& color) {
    Figure tetra;
    tetra.color = color;

    // points of the tetrahedron
    tetra.points = {
        Vector3D::point(1, -1, -1),  // 0
        Vector3D::point(-1, 1, -1),  // 1
        Vector3D::point(1, 1, 1),    // 2
        Vector3D::point(-1, -1, 1)   // 3
    };

    // define faces with point indices
    Face f1, f2, f3, f4;
    f1.point_indices = {0, 1, 2};
    f2.point_indices = {1, 3, 2};
    f3.point_indices = {0, 3, 1};
    f4.point_indices = {0, 2, 3};

    tetra.faces = {f1, f2, f3, f4};

    // define lines for wireframe drawing
    tetra.lines = {
        {0, 1}, {1, 2}, {2, 0},
        {0, 3}, {1, 3}, {2, 3}
    };

    return tetra;
}

Figure createOcta(const std::vector<double>& color) {
    Figure octa;
    octa.color = color;

    // points of the octahedron
    octa.points = {
        Vector3D::point(1, 0, 0),   // 0
        Vector3D::point(0, 1, 0),   // 1
        Vector3D::point(-1, 0, 0),  // 2
        Vector3D::point(0, -1, 0),  // 3
        Vector3D::point(0, 0, -1),  // 4
        Vector3D::point(0, 0, 1)    // 5
    };

    // define faces using point indices
    Face f1, f2, f3, f4, f5, f6, f7, f8;
    f1.point_indices = {0, 1, 5};
    f2.point_indices = {1, 2, 5};
    f3.point_indices = {2, 3, 5};
    f4.point_indices = {3, 0, 5};
    f5.point_indices = {1, 0, 4};
    f6.point_indices = {2, 1, 4};
    f7.point_indices = {3, 2, 4};
    f8.point_indices = {0, 3, 4};

    octa.faces = {f1, f2, f3, f4, f5, f6, f7, f8};

    // define lines for wireframe drawing
    octa.lines = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {0, 5}, {1, 5}, {2, 5}, {3, 5},
        {0, 4}, {1, 4}, {2, 4}, {3, 4}
    };

    return octa;
}

Figure createIco(const std::vector<double>& color) {
    Figure ico;
    ico.color = color;

    const double sqrt5_2 = sqrt(5.0) / 2.0;

    ico.points = {
        Vector3D::point(0, 0, sqrt5_2), // Punt 0 (top)
        Vector3D::point(cos(0*2*M_PI/5), sin(0*2*M_PI/5), 0.5), // Punt 1
        Vector3D::point(cos(1*2*M_PI/5), sin(1*2*M_PI/5), 0.5), // Punt 2
        Vector3D::point(cos(2*2*M_PI/5), sin(2*2*M_PI/5), 0.5), // Punt 3
        Vector3D::point(cos(3*2*M_PI/5), sin(3*2*M_PI/5), 0.5), // Punt 4
        Vector3D::point(cos(4*2*M_PI/5), sin(4*2*M_PI/5), 0.5), // Punt 5
        Vector3D::point(cos(M_PI/5 + 0*2*M_PI/5), sin(M_PI/5 + 0*2*M_PI/5), -0.5), // Punt 6
        Vector3D::point(cos(M_PI/5 + 1*2*M_PI/5), sin(M_PI/5 + 1*2*M_PI/5), -0.5), // Punt 7
        Vector3D::point(cos(M_PI/5 + 2*2*M_PI/5), sin(M_PI/5 + 2*2*M_PI/5), -0.5), // Punt 8
        Vector3D::point(cos(M_PI/5 + 3*2*M_PI/5), sin(M_PI/5 + 3*2*M_PI/5), -0.5), // Punt 9
        Vector3D::point(cos(M_PI/5 + 4*2*M_PI/5), sin(M_PI/5 + 4*2*M_PI/5), -0.5), // Punt 10
        Vector3D::point(0, 0, -sqrt5_2)
    };

    ico.faces = {
        Face{std::vector<int>{0, 1, 2}},
        Face{std::vector<int>{0, 2, 3}},
        Face{std::vector<int>{0, 3, 4}},
        Face{std::vector<int>{0, 4, 5}},
        Face{std::vector<int>{0, 5, 1}},
        Face{std::vector<int>{1, 6, 2}},
        Face{std::vector<int>{2, 6, 7}},
        Face{std::vector<int>{2, 7, 3}},
        Face{std::vector<int>{3, 7, 8}},
        Face{std::vector<int>{3, 8, 4}},
        Face{std::vector<int>{4, 8, 9}},
        Face{std::vector<int>{4, 9, 5}},
        Face{std::vector<int>{5, 9, 10}},
        Face{std::vector<int>{5, 10, 1}},
        Face{std::vector<int>{1, 10, 6}},
        Face{std::vector<int>{11, 7, 6}},
        Face{std::vector<int>{11, 8, 7}},
        Face{std::vector<int>{11, 9, 8}},
        Face{std::vector<int>{11, 10, 9}},
        Face{std::vector<int>{11, 6, 10}}
    };

    return ico;
}

Figure createDodeca(const std::vector<double>& color) {
    Figure dode;
    dode.color = color;

    // gebruik punten van icosaëder
    Figure ico = createIco(color);

    // bereken centroïden voor elk vlak
    for (const Face& face : ico.faces) {
        Vector3D centroid;
        for (int idx : face.point_indices) {
            centroid = centroid + ico.points[idx];
        }
        centroid = centroid * (1.0 / 3.0);
        dode.points.push_back(centroid);
    }

    // definieer vlakken (12 vijfhoeken)
    dode.faces = {
        Face{{0, 1, 2, 3, 4}},
        Face{{0, 5, 6, 7, 1}},
        Face{{1, 7, 8, 9, 2}},
        Face{{2, 9, 10, 11, 3}},
        Face{{3, 11, 12, 13, 4}},
        Face{{4, 13, 14, 5, 0}},
        Face{{19, 18, 17, 16, 15}},
        Face{{19, 14, 13, 12, 18}},
        Face{{18, 12, 11, 10, 17}},
        Face{{17, 10, 9, 8, 16}},
        Face{{16, 8, 7, 6, 15}},
        Face{{15, 6, 5, 14, 19}}
    };

    return dode;
}

Figure createCylinder(int n, double h, const std::vector<double>& color) {
    Figure cyl;
    cyl.color = color;

    // punten voor onder en bovenvlak
    for(int i = 0; i < n; i++) {
        double theta = 2 * M_PI * i / n;
        double x = cos(theta);
        double y = sin(theta);
        cyl.points.push_back(Vector3D::point(x, y, 0)); // onder
        cyl.points.push_back(Vector3D::point(x, y, h)); // boven
    }

    // zijvlakken
    for(int i = 0; i < n; i++) {
        int curr = 2 * i;
        int next = (2 * (i + 1)) % (2 * n);
        cyl.faces.push_back({{curr, next, next + 1, curr + 1}});
    }

    // onder en bovenvlak
    Face bottom, top;
    for(int i = 0; i < n; i++) {
        bottom.point_indices.push_back(2 * i);
        top.point_indices.push_back(2 * i + 1);
    }
    reverse(top.point_indices.begin(), top.point_indices.end());
    cyl.faces.push_back(bottom);
    cyl.faces.push_back(top);

    return cyl;
}

Figure createCone(const int n, const double h, const std::vector<double>& color) {
    Figure cone;
    cone.color = color;

    for (int i = 0; i < n; i++) {
        double angle = 2 * M_PI * i / n;
        double x = cos(angle);
        double y = sin(angle);


        cone.points.push_back(Vector3D::point(x, y, 0));
    }


    int top_index = cone.points.size();
    cone.points.push_back(Vector3D::point(0, 0, h));

    for (int i = 0; i < n; i++) {
        int current = i;
        int next = (i + 1) % n;

        cone.faces.push_back(Face{{current, next, top_index}});
    }

    std::vector<int> bottom_face;
    for (int i = n-1; i >= 0; i--) {
        bottom_face.push_back(i);
    }
    cone.faces.push_back(Face{bottom_face});

    return cone;
}

Figure createTorus(const double r, const double R, const int n, const int m, const std::vector<double>& color) {
    Figure torus;
    torus.color = color;

    // genereer punten
    for (int i = 0; i < n; i++) {
        double u = 2 * M_PI * i / n;
        for (int j = 0; j < m; j++) {
            double v = 2 * M_PI * j / m;

            // bereken punten volgens parametervergelijking
            double x = (R + r * cos(v)) * cos(u);
            double y = (R + r * cos(v)) * sin(u);
            double z = r * sin(v);

            torus.points.push_back(Vector3D::point(x, y, z));
        }
    }

    // genereer faces (vierhoeken)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            int current = i * m + j;
            int next_i = (i + 1) % n;
            int next_j = (j + 1) % m;

            int p1 = current;
            int p2 = next_i * m + j;
            int p3 = next_i * m + next_j;
            int p4 = i * m + next_j;

            // voeg vierhoek toe in tegenwijzerzin
            torus.faces.push_back(Face{{p1, p2, p3, p4}});
        }
    }

    return torus;
}

Figure createSphere(int iterations, const std::vector<double>& color) {
    Figure sphere = createIco(color); // Start met een icosaëder

    for (int it = 0; it < iterations; it++) {
        std::vector<Face> newFaces;
        std::map<std::pair<int, int>, int> edgeMidpointCache;

        // Hulpfunctie om middenpunten te vinden of te genereren
        auto getMidpoint = [&](int p1Idx, int p2Idx) {
            if (p1Idx > p2Idx) std::swap(p1Idx, p2Idx);
            auto key = std::make_pair(p1Idx, p2Idx);

            if (edgeMidpointCache.count(key)) {
                return edgeMidpointCache[key];
            }

            Vector3D p1 = sphere.points[p1Idx];
            Vector3D p2 = sphere.points[p2Idx];
            Vector3D midpoint = (p1 + p2) * 0.5;
            midpoint.normalise(); // Zorg dat het punt op de bol ligt

            sphere.points.push_back(midpoint);
            int newIdx = sphere.points.size() - 1;
            edgeMidpointCache[key] = newIdx;
            return newIdx;
        };

        // Subdivide elk driehoekig vlak
        for (const Face& face : sphere.faces) {
            int a = face.point_indices[0];
            int b = face.point_indices[1];
            int c = face.point_indices[2];

            int ab = getMidpoint(a, b);
            int bc = getMidpoint(b, c);
            int ca = getMidpoint(c, a);

            // Voeg 4 nieuwe driehoeken toe
            newFaces.push_back(Face{{a, ab, ca}});
            newFaces.push_back(Face{{b, bc, ab}});
            newFaces.push_back(Face{{c, ca, bc}});
            newFaces.push_back(Face{{ab, bc, ca}});
        }
        sphere.faces = newFaces;
    }

    // Genereer lijnen voor wireframe
    for (const Face& face : sphere.faces) {
        for (size_t i = 0; i < face.point_indices.size(); i++) {
            int j = (i + 1) % face.point_indices.size();
            sphere.lines.push_back({face.point_indices[i], face.point_indices[j]});
        }
    }

    return sphere;
}

// Figure parse3DLSystem(const std::string& filename, const std::vector<double>& color) {
//     std::ifstream in(filename);
//     LParser::LSystem3D lSystem;
//     in >> lSystem;
//     in.close();
//
//     Figure figure;
//     figure.color = color;
//
//     std::string current = lSystem.get_initiator();
//     for (unsigned int i = 0; i < lSystem.get_nr_iterations(); i++) {
//         std::string next;
//         for (char c : current) {
//             if (lSystem.get_alphabet().count(c)) {
//                 next += lSystem.get_replacement(c);
//             } else {
//                 next += c;
//             }
//         }
//         current = next;
//     }
//
//     Vector3D position(0, 0, 0);
//     Vector3D H(1, 0, 0); // Heading
//     Vector3D L(0, 1, 0); // Left
//     Vector3D U(0, 0, 1); // Up
//     std::stack<std::tuple<Vector3D, Vector3D, Vector3D, Vector3D>> stateStack;
//     double angle = lSystem.get_angle() * (M_PI / 180.0);
//
//     std::vector<Vector3D> points;
//     std::vector<std::pair<int, int>> lines;
//     int pointIdx = 0;
//
//     for (char c : current) {
//         switch (c) {
//             case '+': // Draai linksom om H-as
//                 L = L * cos(angle) + U * sin(angle);
//                 U = L.cross(H);
//                 break;
//             case '-': // Draai rechtsom om H-as
//                 L = L * cos(-angle) + U * sin(-angle);
//                 U = L.cross(H);
//                 break;
//             case '^': // Draai om L-as
//                 H = H * cos(angle) + U * sin(angle);
//                 U = H.cross(L);
//                 break;
//             case '&': // Draai om L-as (negatief)
//                 H = H * cos(-angle) + U * sin(-angle);
//                 U = H.cross(L);
//                 break;
//             case '\\': // Roll linksom om U-as
//                 H = H * cos(angle) + L * sin(angle);
//                 L = H.cross(U);
//                 break;
//             case '/': // Roll rechtsom om U-as
//                 H = H * cos(-angle) + L * sin(-angle);
//                 L = H.cross(U);
//                 break;
//             case '|': // Keer richting om
//                 H = H * -1;
//                 L = L * -1;
//                 break;
//             case '(': // Sla staat op
//                 stateStack.push({position, H, L, U});
//                 break;
//             case ')': // Herstel staat
//                 if (!stateStack.empty()) {
//                     auto [savedPos, savedH, savedL, savedU] = stateStack.top();
//                     stateStack.pop();
//                     position = savedPos;
//                     H = savedH;
//                     L = savedL;
//                     U = savedU;
//                 }
//                 break;
//             default:
//                 if (lSystem.draw(c)) {
//                     Vector3D newPos = position + H;
//                     points.push_back(position);
//                     points.push_back(newPos);
//                     lines.push_back({pointIdx, pointIdx + 1});
//                     pointIdx += 2;
//                     position = newPos;
//                 } else {
//                     position = position + H;
//                 }
//                 break;
//         }
//     }
//
//     figure.points = points;
//     figure.lines = lines;
//     return figure;
// }

std::vector<Line2D> generateLinesFromLSystem(const LParser::LSystem2D &lsystem, const std::vector<double> &color) {
    std::vector<Line2D> lines;
    std::stack<std::tuple<double, double, double> > stateStack; // For branching
    std::string current = lsystem.get_initiator();

    for (unsigned int i = 0; i < lsystem.get_nr_iterations(); ++i) {
        std::string next;
        for (char c: current) {
            if (lsystem.get_alphabet().count(c)) {
                next += lsystem.get_replacement(c);
            } else {
                next += c;
            }
        }
        current = next;
    }
    double x = 0.0, y = 0.0;
    double angle = (-lsystem.get_starting_angle()) * (M_PI / 180.0);
    const double delta = lsystem.get_angle() * (M_PI / 180.0);

    for (char c: current) {
        if (c == '(') {
            // save current state
            stateStack.push(std::make_tuple(x, y, angle));
        } else if (c == ')') {
            // restore previous state
            if (!stateStack.empty()) {
                auto [savedX, savedY, savedAngle] = stateStack.top();
                stateStack.pop();
                x = savedX;
                y = savedY;
                angle = savedAngle;
            }
        } else if (lsystem.get_alphabet().count(c)) {
            bool drawLine = lsystem.draw(c);
            double newX = x + cos(angle);
            double newY = y + sin(angle);
            if (drawLine) {
                Line2D line;
                line.start = {x, y};
                line.end = {newX, newY};
                line.color = color;
                line.z1 = 1.0;
                line.z2 = 1.0;
                lines.push_back(line);
            }
            x = newX;
            y = newY;
        } else if (c == '+') {
            angle -= delta;
        } else if (c == '-') {
            angle += delta;
        }
    }

    return lines;
}

Figure parseFigure(const ini::Configuration& configuration, int index) {
    std::string sectionName = "Figure" + std::to_string(index);
    Figure figure;

    std::string type = configuration[sectionName]["type"].as_string_or_die();
    figure.color = configuration[sectionName]["color"].as_double_tuple_or_die();

    if (type == "Cube") {
        return createCube(figure.color);
    }
    else if (type == "Tetrahedron") {
        return createTetra(figure.color);
    }
    else if (type == "Octahedron") {
        return createOcta(figure.color);
    }
    else if (type == "Icosahedron") {
        return createIco(figure.color);
    }
    else if (type == "Dodecahedron") {
        return createDodeca(figure.color);
    }
    else if (type == "Cylinder") {
        int n = configuration[sectionName]["n"].as_int_or_die();
        double h = configuration[sectionName]["height"].as_double_or_die();
        return createCylinder(n, h, figure.color);
    }
    else if (type == "Cone") {
        int n = configuration[sectionName]["n"].as_int_or_die();
        double h = configuration[sectionName]["height"].as_double_or_die();
        return createCone(n, h, figure.color);
    }
    else if (type == "Sphere") {
        int n = configuration[sectionName]["n"].as_int_or_die();
        return createSphere(n, figure.color);
    }
    else if (type == "Torus") {
        double R = configuration[sectionName]["R"].as_double_or_die();
        double r = configuration[sectionName]["r"].as_double_or_die();
        int n = configuration[sectionName]["n"].as_int_or_die();
        int m = configuration[sectionName]["m"].as_int_or_die();
        return createTorus(R, r, n, m, figure.color);
    }

    // else if (type == "3DLSystem") {
    //     std::string inputfile = configuration[sectionName]["inputfile"].as_string_or_die();
    //     return parse3DLSystem(inputfile, figure.color);
    // }

    else if (type == "LineDrawing") {
        // parse points for custom figure
        int nrPoints = configuration[sectionName]["nrPoints"].as_int_or_die();
        figure.points.resize(nrPoints);

        for (int i = 0; i < nrPoints; ++i) {
            std::string pointName = "point" + std::to_string(i);
            std::vector<double> point = configuration[sectionName][pointName].as_double_tuple_or_die();
            figure.points[i] = Vector3D::point(point[0], point[1], point[2]);
        }

        // parse lines
        int nrLines = configuration[sectionName]["nrLines"].as_int_or_die();
        figure.lines.resize(nrLines);

        for (int i = 0; i < nrLines; ++i) {
            std::string lineName = "line" + std::to_string(i);
            std::vector<int> line = configuration[sectionName][lineName].as_int_tuple_or_die();
            figure.lines[i] = std::make_pair(line[0], line[1]);
        }

        int nrFaces = configuration[sectionName]["nrFaces"].as_int_or_default(-1);
        if (nrFaces >= 0) {
            figure.faces.resize(nrFaces);

            for (int i = 0; i < nrFaces; ++i) {
                std::string faceName = "face" + std::to_string(i);
                std::vector<int> faceIndices = configuration[sectionName][faceName].as_int_tuple_or_die();
                Face face;
                face.point_indices = faceIndices;
                figure.faces[i] = face;
            }
        }
    }

    return figure;
}

void transformFigure(Figure& figure, const ini::Configuration& configuration, int index) {
    std::string sectionName = "Figure" + std::to_string(index);

    double scale = configuration[sectionName]["scale"].as_double_or_die();
    double rotateX = configuration[sectionName]["rotateX"].as_double_or_die();
    double rotateY = configuration[sectionName]["rotateY"].as_double_or_die();
    double rotateZ = configuration[sectionName]["rotateZ"].as_double_or_die();
    std::vector<double> center = configuration[sectionName]["center"].as_double_tuple_or_die();

    Matrix S = scalingMatrix(scale);
    Matrix Rx = rotationX(rotateX);
    Matrix Ry = rotationY(rotateY);
    Matrix Rz = rotationZ(rotateZ);
    Matrix T = translationMatrix(Vector3D::point(center[0], center[1], center[2]));
    Matrix transform = S * Rx * Ry * Rz * T;

    for (Vector3D& point : figure.points) {
        point = point * transform;
    }
}
img::EasyImage draw2DLines(const std::vector<Line2D>& lines, int size, const std::vector<double>& bgColor) {
    if (lines.empty()) {
        std::cerr << "Lines vector is empty" << std::endl;
        return img::EasyImage(1, 1, img::Color(0, 0, 0));
    }
    double xmin = lines[0].start.x, xmax = lines[0].start.x;
    double ymin = lines[0].start.y, ymax = lines[0].start.y;
    for (const auto& line : lines) {
        xmin = std::min({xmin, line.start.x, line.end.x});
        xmax = std::max({xmax, line.start.x, line.end.x});
        ymin = std::min({ymin, line.start.y, line.end.y});
        ymax = std::max({ymax, line.start.y, line.end.y});
    }
    double xrange = xmax - xmin;
    double yrange = ymax - ymin;
    if (xrange == 0) xrange = 1;
    if (yrange == 0) yrange = 1;

    int imageWidth, imageHeight;
    if (xrange > yrange) {
        imageWidth = size;
        imageHeight = static_cast<int>(size * (yrange / xrange));
    } else {
        imageHeight = size;
        imageWidth = static_cast<int>(size * (xrange / yrange));
    }

    double scale = 0.95 * std::min(imageWidth / xrange, imageHeight / yrange);
    double dx = (imageWidth / 2.0) - (scale * (xmin + xmax)) / 2.0;
    double dy = (imageHeight / 2.0) - (scale * (ymin + ymax)) / 2.0;

    img::EasyImage image(imageWidth, imageHeight);
    std::vector<int> bgScaled = scaleColor(bgColor);
    image.clear(img::Color(bgScaled[0], bgScaled[1], bgScaled[2]));

    for (const auto& line : lines) {
        Point2D p1 = line.start;
        Point2D p2 = line.end;
        p1.x = p1.x * scale + dx;
        p1.y = p1.y * scale + dy;
        p2.x = p2.x * scale + dx;
        p2.y = p2.y * scale + dy;
        int y1 = image.get_height() - static_cast<int>(std::round(p1.y)) - 1;
        int y2 = image.get_height() - static_cast<int>(std::round(p2.y)) - 1;
        int x1 = static_cast<int>(std::round(p1.x));
        int x2 = static_cast<int>(std::round(p2.x));

        std::vector<int> lineColor = scaleColor(line.color);
        img::Color color(lineColor[0], lineColor[1], lineColor[2]);
        image.draw_line(x1, y1, x2, y2, color);
    }

    return image;
}

img::EasyImage generate_image(const ini::Configuration &configuration) {
    std::string type = configuration["General"]["type"].as_string_or_die();
    int size = configuration["General"]["size"].as_int_or_die();
    std::vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

    if (type == "Wireframe") {
        // parse eye position
        std::vector<double> eyePos = configuration["General"]["eye"].as_double_tuple_or_die();
        Vector3D eye = Vector3D::point(eyePos[0], eyePos[1], eyePos[2]);

        // parse figures
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        std::vector<Figure> figures(nrFigures);

        for (int i = 0; i < nrFigures; ++i) {
            figures[i] = parseFigure(configuration, i);
            transformFigure(figures[i], configuration, i);
        }

        // apply eye transformation
        Matrix eyeMatrix = eyePointTrans(eye);

        // transform all points to eye coordinates
        for (auto& figure : figures) {
            for (auto& point : figure.points) {
                point = point * eyeMatrix;
            }
        }

        // generate 2D lines through projection
        std::vector<Line2D> lines2D = doProjection(figures, 1.0);

        // draw the 2D lines
        return draw2DLines(lines2D, size, backgroundColor);
    }
    else if (type == "2DLSystem") {
        std::string inputfile = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        std::vector<double> color = configuration["2DLSystem"]["color"].as_double_tuple_or_die();

        std::ifstream input_stream(inputfile);
        LParser::LSystem2D lsystem(input_stream);

        std::vector<Line2D> lines2D = generateLinesFromLSystem(lsystem, color);
        return draw2DLines(lines2D, size, backgroundColor);
    }

    return img::EasyImage();
}

int main(int argc, char const *argv[]) {
    int retVal = 0;
    try {
        std::vector<std::string> args = std::vector<std::string>(argv + 1, argv + argc);
        if (args.empty()) {
            std::ifstream fileIn("filelist");
            std::string filelistName;
            while (std::getline(fileIn, filelistName)) {
                args.push_back(filelistName);
            }
        }
        for (std::string fileName: args) {
            ini::Configuration conf;
            try {
                std::ifstream fin(fileName);
                if (fin.peek() == std::istream::traits_type::eof()) {
                    std::cout << "Ini file appears empty. Does '" << fileName << "' exist?" << std::endl;
                    continue;
                }
                fin >> conf;
                fin.close();
            } catch (ini::ParseException &ex) {
                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                retVal = 1;
                continue;
            }

            img::EasyImage image = generate_image(conf);
            if (image.get_height() > 0 && image.get_width() > 0) {
                std::string::size_type pos = fileName.rfind('.');
                if (pos == std::string::npos) {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    fileName += ".bmp";
                } else {
                    fileName = fileName.substr(0, pos) + ".bmp";
                }
                try {
                    std::ofstream f_out(fileName.c_str(), std::ios::trunc | std::ios::out | std::ios::binary);
                    f_out << image;
                } catch (std::exception &ex) {
                    std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                    retVal = 1;
                }
            } else {
                std::cout << "Could not generate image for " << fileName << std::endl;
            }
        }
    } catch (const std::bad_alloc &exception) {
        //When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
        //Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
        //(Unless of course you are already consuming the maximum allowed amount of memory)
        //If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
        //mark the test as failed while in reality it just needed a bit more memory
        std::cerr << "Error: insufficient memory" << std::endl;
        retVal = 100;
    }
    return retVal;
}