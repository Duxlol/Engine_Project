//
// Created by kali on 5/14/25.
//

#ifndef ZBUFFER_H
#define ZBUFFER_H
#include <limits>
#include <vector>

class ZBuffer {
private:
    std::vector<std::vector<double>> buffer;
    unsigned int width, height;

public:
    ZBuffer(unsigned int width, unsigned int height);

    void set(unsigned int x, unsigned int y, double z);

    double get(unsigned int x, unsigned int y) const;
};



#endif //ZBUFFER_H
