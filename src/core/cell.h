#pragma once

enum class CellType {
    orthorhombic,
    triclinic,
    hexagonal,
    rhombic_dodecahedron
};

struct Cell {
    double matrix[9]; // 3x3, row-major
    CellType type;

    Cell();
    Cell(double a, double b, double c);
    Cell(const double m[9], CellType t);

    double volume() const noexcept;
    void minimum_image(double dx, double dy, double dz,
                       double& ix, double& iy, double& iz) const noexcept;
    void wrap(double& x, double& y, double& z) const noexcept;
    void scale(double sx, double sy, double sz) noexcept;
};
