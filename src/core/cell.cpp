#include "cell.h"
#include <cmath>
#include <cstring>

Cell::Cell() : type(CellType::orthorhombic) {
    std::memset(matrix, 0, sizeof(matrix));
    matrix[0] = 1.0; matrix[4] = 1.0; matrix[8] = 1.0;
}

Cell::Cell(double a, double b, double c) : type(CellType::orthorhombic) {
    std::memset(matrix, 0, sizeof(matrix));
    matrix[0] = a; matrix[4] = b; matrix[8] = c;
}

Cell::Cell(const double m[9], CellType t) : type(t) {
    std::memcpy(matrix, m, sizeof(matrix));
}

double Cell::volume() const noexcept {
    double a1 = matrix[0], a2 = matrix[1], a3 = matrix[2];
    double b1 = matrix[3], b2 = matrix[4], b3 = matrix[5];
    double c1 = matrix[6], c2 = matrix[7], c3 = matrix[8];
    return a1*(b2*c3 - b3*c2) - a2*(b1*c3 - b3*c1) + a3*(b1*c2 - b2*c1);
}

void Cell::minimum_image(double dx, double dy, double dz,
                          double& ix, double& iy, double& iz) const noexcept {
    if (type == CellType::orthorhombic) {
        double inv_a = 1.0 / matrix[0];
        double inv_b = 1.0 / matrix[4];
        double inv_c = 1.0 / matrix[8];
        ix = dx - std::round(dx * inv_a) * matrix[0];
        iy = dy - std::round(dy * inv_b) * matrix[4];
        iz = dz - std::round(dz * inv_c) * matrix[8];
    } else {
        ix = dx; iy = dy; iz = dz; // stub for non-orthogonal types
    }
}

void Cell::wrap(double& x, double& y, double& z) const noexcept {
    if (type == CellType::orthorhombic) {
        x -= std::floor(x / matrix[0]) * matrix[0];
        y -= std::floor(y / matrix[4]) * matrix[4];
        z -= std::floor(z / matrix[8]) * matrix[8];
    }
}

void Cell::scale(double sx, double sy, double sz) noexcept {
    matrix[0] *= sx; matrix[4] *= sy; matrix[8] *= sz;
}
