//
// Created by kali on 5/14/25.
//

#include "zbuffer.h"

ZBuffer::ZBuffer(unsigned int width, unsigned int height)
        : width(width), height(height), buffer(height, std::vector<double>(width, std::numeric_limits<double>::infinity())) {}

void ZBuffer::set(unsigned int x, unsigned int y, double z) {
    if (x < width && y < height && z < buffer[y][x]) {
        buffer[y][x] = z;
    }
}

double ZBuffer::get(unsigned int x, unsigned int y) const {
    if (x < width && y < height) {
        return buffer[y][x];
    }
    return std::numeric_limits<double>::infinity();
}