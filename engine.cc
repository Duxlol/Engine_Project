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

struct Point2D {
    double x, y;
};

struct Line2D {
    Point2D start, end;
    std::vector<double> color;
};

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
    R(4,4) = 1;

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

// Generates an eye transformation matrix
Matrix eyeTransformationMatrix(const Vector3D& eye) {
    // Calculate spherical coordinates from Cartesian
    double r = sqrt(eye.x * eye.x + eye.y * eye.y + eye.z * eye.z);
    double theta = atan2(eye.y, eye.x);
    double phi = acos(eye.z / r);

    Matrix eyepointMatrix;
    eyepointMatrix(1, 1) = -sin(theta);
    eyepointMatrix(1, 2) = -cos(theta)*cos(phi);
    eyepointMatrix(1, 3) = cos(theta)*sin(phi);
    eyepointMatrix(1,4) = 0;

    eyepointMatrix(2, 1) = cos(theta);
    eyepointMatrix(2, 2) = -sin(theta)*cos(phi);
    eyepointMatrix(2, 3) = sin(theta)*sin(phi);
    eyepointMatrix(2,4) = 0;

    eyepointMatrix(3,1) = 0;
    eyepointMatrix(3, 2) = sin(phi);
    eyepointMatrix(3, 3) = cos(phi);
    eyepointMatrix(3,4) = 0;

    eyepointMatrix(4,1) = 0;
    eyepointMatrix(4,2) = 0;
    eyepointMatrix(4,3) = -r;
    eyepointMatrix(4,4) = 1;

    return eyepointMatrix;
}

// Perspective projection (assumes point is already in eye coordinates)
Point2D projectPoint(const Vector3D& point) {
    double d = 1.0; // Fixed distance value as specified in the assignment

    // Perspective projection formulas
    double x = (d * point.x) / -point.z;
    double y = (d * -point.y) / -point.z;

    return {x, y};
}

std::vector<int> scaleColor(const std::vector<double> &color) {
    std::vector<int> scaled(3);
    for (int i = 0; i < 3; ++i) {
        scaled[i] = static_cast<int>(std::round(color[i] * 255));
        scaled[i] = std::clamp(scaled[i], 0, 255);
    }
    return scaled;
}

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
            // Save current state
            stateStack.push(std::make_tuple(x, y, angle));
        } else if (c == ')') {
            // Restore previous state
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
                lines.push_back({{x, y}, {newX, newY}, color});
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

    // Parse color
    figure.color = configuration[sectionName]["color"].as_double_tuple_or_die();

    // Parse points
    int nrPoints = configuration[sectionName]["nrPoints"].as_int_or_die();
    figure.points.resize(nrPoints);

    for (int i = 0; i < nrPoints; ++i) {
        std::string pointName = "point" + std::to_string(i);
        std::vector<double> point = configuration[sectionName][pointName].as_double_tuple_or_die();
        figure.points[i] = Vector3D::point(point[0], point[1], point[2]);
    }

    // Parse lines
    int nrLines = configuration[sectionName]["nrLines"].as_int_or_die();
    figure.lines.resize(nrLines);

    for (int i = 0; i < nrLines; ++i) {
        std::string lineName = "line" + std::to_string(i);
        std::vector<int> line = configuration[sectionName][lineName].as_int_tuple_or_die();
        figure.lines[i] = std::make_pair(line[0], line[1]);
    }

    return figure;
}

// Apply all transformations to a figure
void transformFigure(Figure& figure, const ini::Configuration& configuration, int index) {
    std::string sectionName = "Figure" + std::to_string(index);

    double scale = configuration[sectionName]["scale"].as_double_or_die();
    double rotateX = configuration[sectionName]["rotateX"].as_double_or_die();
    double rotateY = configuration[sectionName]["rotateY"].as_double_or_die();
    double rotateZ = configuration[sectionName]["rotateZ"].as_double_or_die();
    std::vector<double> center = configuration[sectionName]["center"].as_double_tuple_or_die();

    // Create transformation matrices
    Matrix S = scalingMatrix(scale);
    Matrix Rx = rotationX(rotateX);
    Matrix Ry = rotationY(rotateY);
    Matrix Rz = rotationZ(rotateZ);
    Matrix T = translationMatrix(Vector3D::point(center[0], center[1], center[2]));

    // Combined transformation matrix
    Matrix transform = S * Rx * Ry * Rz * T;

    // Apply transformation to all points
    for (Vector3D& point : figure.points) {
        point = point * transform;
    }
}

std::vector<Line2D> applyEyeAndProjection(const std::vector<Figure>& figures, const Vector3D& eye) {
    // Create eye transformation matrix
    Matrix eyeMatrix = eyeTransformationMatrix(eye);

    std::vector<Line2D> lines2D;

    for (const Figure& figure : figures) {
        // Transform all points to eye coordinates
        std::vector<Vector3D> transformedPoints;
        transformedPoints.reserve(figure.points.size());

        for (const Vector3D& point : figure.points) {
            transformedPoints.push_back(point * eyeMatrix);
        }

        // Project all points and create 2D lines
        for (const std::pair<int, int>& line : figure.lines) {
            int start = line.first;
            int end = line.second;

            // Skip lines with negative z (behind the eye)
            if (transformedPoints[start].z >= 0 || transformedPoints[end].z >= 0) continue;

            Point2D startPoint = projectPoint(transformedPoints[start]);
            Point2D endPoint = projectPoint(transformedPoints[end]);

            lines2D.push_back({startPoint, endPoint, figure.color});
        }
    }

    return lines2D;
}

img::EasyImage draw2DLines(const std::vector<Line2D>& lines, int size, const std::vector<double>& bgColor) {
    if (lines.empty()) {
        std::cerr << "lines is empty" << std::endl;
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
        // Parse eye position
        std::vector<double> eyePos = configuration["General"]["eye"].as_double_tuple_or_die();
        Vector3D eye = Vector3D::point(eyePos[0], eyePos[1], eyePos[2]);

        // Parse figures
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        std::vector<Figure> figures(nrFigures);

        for (int i = 0; i < nrFigures; ++i) {
            figures[i] = parseFigure(configuration, i);
            transformFigure(figures[i], configuration, i);
        }

        // Apply eye transformation and projection
        std::vector<Line2D> lines2D = applyEyeAndProjection(figures, eye);


        // Draw the 2D lines
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
                    std::cerr << "Failed to write  image to file: " << ex.what() << std::endl;
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