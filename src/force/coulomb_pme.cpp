#include "coulomb_pme.h"
#include "core/cell.h"
#include <cmath>
#include <algorithm>
#include <cstring>

CoulombPME::CoulombPME(PMEParams params)
    : params_(std::move(params))
{
    int nx = params_.nx, ny = params_.ny, nz = params_.nz;
    int nx_half = nx / 2 + 1;
    size_t grid_size = static_cast<size_t>(nz) * ny * (2 * nx_half);
    grid_.resize(grid_size, 0.0);

    plan_fwd_ = fftw_plan_dft_r2c_3d(
        nz, ny, nx,
        grid_.data(),
        reinterpret_cast<fftw_complex*>(grid_.data()),
        FFTW_ESTIMATE);

    plan_bwd_ = fftw_plan_dft_c2r_3d(
        nz, ny, nx,
        reinterpret_cast<fftw_complex*>(grid_.data()),
        grid_.data(),
        FFTW_ESTIMATE);
}

CoulombPME::~CoulombPME() {
    fftw_destroy_plan(plan_fwd_);
    fftw_destroy_plan(plan_bwd_);
}

std::string CoulombPME::name() const { return "CoulombPME"; }

void CoulombPME::set_charges(std::span<const double> charges) {
    charges_.assign(charges.begin(), charges.end());
    self_energy_ = calc_self_energy(params_.ewald_coefficient, charges_, k_e_);
}

// ---- B-spline M_4 (cubic cardinal B-spline) and its derivative ----
// Support: |u| < 2, non-zero on 4 consecutive grid points

double CoulombPME::theta4(double u) {
    double a = std::abs(u);
    if (a >= 2.0) return 0.0;
    if (a >= 1.0) {
        double t = 2.0 - a;
        return t * t * t / 6.0;
    }
    return (4.0 - 6.0 * a * a + 3.0 * a * a * a) / 6.0;
}

double CoulombPME::dtheta4(double u) {
    double a = std::abs(u);
    if (a >= 2.0) return 0.0;
    double s = (u >= 0.0) ? 1.0 : -1.0;
    if (a >= 1.0) {
        double t = 2.0 - a;
        return -s * t * t / 2.0;
    }
    return s * (9.0 * a * a - 12.0 * a) / 6.0;
}

double CoulombPME::calc_self_energy(double alpha, std::span<const double> charges, double k_e) {
    double sum_q2 = 0.0;
    for (double q : charges) sum_q2 += q * q;
    return -k_e * alpha / std::sqrt(M_PI) * sum_q2;
}

// ---- Spread point charges onto grid using B-spline interpolation ----

void CoulombPME::spread_charges(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    const Cell& cell)
{
    size_t n = pos_x.size();
    double Lx = cell.matrix[0];
    double Ly = cell.matrix[4];
    double Lz = cell.matrix[8];
    int nx = params_.nx, ny = params_.ny, nz = params_.nz;

    std::fill(grid_.begin(), grid_.end(), 0.0);

    for (size_t i = 0; i < n; ++i) {
        double q = charges_[i];
        if (q == 0.0) continue;

        double frac_x = pos_x[i] / Lx;
        double frac_y = pos_y[i] / Ly;
        double frac_z = pos_z[i] / Lz;

        frac_x -= std::floor(frac_x);
        frac_y -= std::floor(frac_y);
        frac_z -= std::floor(frac_z);

        double ux = frac_x * nx;
        double uy = frac_y * ny;
        double uz = frac_z * nz;

        int m0x = static_cast<int>(std::floor(ux)) - 1;
        int m0y = static_cast<int>(std::floor(uy)) - 1;
        int m0z = static_cast<int>(std::floor(uz)) - 1;

        double wx[4], wy[4], wz[4];
        for (int d = 0; d < 4; ++d) {
            wx[d] = theta4(ux - (m0x + d));
            wy[d] = theta4(uy - (m0y + d));
            wz[d] = theta4(uz - (m0z + d));
        }

        for (int iz = 0; iz < 4; ++iz) {
            int gz = (m0z + iz) % nz;
            if (gz < 0) gz += nz;
            double w_z = wz[iz];
            for (int iy = 0; iy < 4; ++iy) {
                int gy = (m0y + iy) % ny;
                if (gy < 0) gy += ny;
                double w_yz = wy[iy] * w_z;
                for (int ix = 0; ix < 4; ++ix) {
                    int gx = (m0x + ix) % nx;
                    if (gx < 0) gx += nx;
                    grid_[gz * ny * nx + gy * nx + gx] += q * wx[ix] * w_yz;
                }
            }
        }
    }
}

// ---- Modify Q̃(k) by Green's function in Fourier space, compute E_rec ----
// Uses FFTW r2c: only kx = 0 … nx/2 stored (last dim truncated)
// Hermitian pairing: multiply stored |Q̃|² by w = 2 for 0 < kx < nx/2,
// w = 1 for kx = 0 or kx = nx/2 (Nyquist, self-paired).

double CoulombPME::compute_reciprocal_energy_and_potential(const Cell& cell) {
    int nx = params_.nx, ny = params_.ny, nz = params_.nz;
    int nz_half = nz / 2 + 1;
    int nx_half = nx / 2 + 1;
    double alpha = params_.ewald_coefficient;
    double V = cell.volume();
    double Lx = cell.matrix[0];
    double Ly = cell.matrix[4];
    double Lz = cell.matrix[8];
    int spline_order = params_.spline_order;

    fftw_execute(plan_fwd_);
    fftw_complex* grid_f = reinterpret_cast<fftw_complex*>(grid_.data());

    double energy = 0.0;

    auto b_factor = [&](int m, int K) -> double {
        if (m == 0) return 1.0;
        double t = M_PI * std::abs(m) / K;
        return std::pow(std::sin(t) / t, 2 * spline_order);
    };

    for (int kz = 0; kz < nz; ++kz) {
        int mz = (kz <= nz / 2) ? kz : kz - nz;
        for (int ky = 0; ky < ny; ++ky) {
            int my = (ky <= ny / 2) ? ky : ky - ny;
            size_t row_base = static_cast<size_t>(kz) * ny * nx_half
                            + static_cast<size_t>(ky) * nx_half;
            for (int kx = 0; kx <= nx / 2; ++kx) {
                int mx = kx;

                if (mx == 0 && my == 0 && mz == 0) continue;

                double w = (kx == 0 || (nx % 2 == 0 && kx == nx / 2)) ? 1.0 : 2.0;

                double msq = (static_cast<double>(mx) * mx) / (Lx * Lx)
                           + (static_cast<double>(my) * my) / (Ly * Ly)
                           + (static_cast<double>(mz) * mz) / (Lz * Lz);
                double k2 = 4.0 * M_PI * M_PI * msq;

                double B = b_factor(mx, nx) * b_factor(my, ny) * b_factor(mz, nz);

                double green = 4.0 * M_PI * k_e_ / V * std::exp(-k2 / (4.0 * alpha * alpha)) * B / k2;

                size_t idx = row_base + static_cast<size_t>(kx);
                double re = grid_f[idx][0];
                double im = grid_f[idx][1];
                energy += 0.5 * w * green * (re * re + im * im);

                grid_f[idx][0] = re * green;
                grid_f[idx][1] = im * green;
            }
        }
    }

    fftw_execute(plan_bwd_);

    return energy;
}

// ---- Interpolate forces from reciprocal potential grid ----

void CoulombPME::interpolate_forces(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    const Cell& cell,
    std::span<double> forces_x,
    std::span<double> forces_y,
    std::span<double> forces_z,
    double fx_factor,
    double fy_factor,
    double fz_factor)
{
    size_t n = pos_x.size();
    double Lx = cell.matrix[0];
    double Ly = cell.matrix[4];
    double Lz = cell.matrix[8];
    int nx = params_.nx, ny = params_.ny, nz = params_.nz;

    for (size_t i = 0; i < n; ++i) {
        double q = charges_[i];
        if (q == 0.0) continue;

        double frac_x = pos_x[i] / Lx;
        double frac_y = pos_y[i] / Ly;
        double frac_z = pos_z[i] / Lz;

        frac_x -= std::floor(frac_x);
        frac_y -= std::floor(frac_y);
        frac_z -= std::floor(frac_z);

        double ux = frac_x * nx;
        double uy = frac_y * ny;
        double uz = frac_z * nz;

        int m0x = static_cast<int>(std::floor(ux)) - 1;
        int m0y = static_cast<int>(std::floor(uy)) - 1;
        int m0z = static_cast<int>(std::floor(uz)) - 1;

        double wx[4], wy[4], wz[4];
        double dwx[4], dwy[4], dwz[4];
        for (int d = 0; d < 4; ++d) {
            wx[d] = theta4(ux - (m0x + d));
            wy[d] = theta4(uy - (m0y + d));
            wz[d] = theta4(uz - (m0z + d));
            dwx[d] = dtheta4(ux - (m0x + d));
            dwy[d] = dtheta4(uy - (m0y + d));
            dwz[d] = dtheta4(uz - (m0z + d));
        }

        double fx = 0.0, fy = 0.0, fz = 0.0;

        for (int iz = 0; iz < 4; ++iz) {
            int gz = (m0z + iz) % nz;
            if (gz < 0) gz += nz;
            double w_z = wz[iz];
            double dw_z = dwz[iz];
            for (int iy = 0; iy < 4; ++iy) {
                int gy = (m0y + iy) % ny;
                if (gy < 0) gy += ny;
                double w_y = wy[iy];
                double dw_y = dwy[iy];
                size_t row_base = static_cast<size_t>(gz) * ny * nx
                                + static_cast<size_t>(gy) * nx;
                for (int ix = 0; ix < 4; ++ix) {
                    int gx = (m0x + ix) % nx;
                    if (gx < 0) gx += nx;
                    double phi = grid_[row_base + static_cast<size_t>(gx)];
                    fx += phi * dwx[ix] * w_y * w_z;
                    fy += phi * wx[ix] * dw_y * w_z;
                    fz += phi * wx[ix] * w_y * dw_z;
                }
            }
        }

        forces_x[i] -= q * fx_factor * fx;
        forces_y[i] -= q * fy_factor * fy;
        forces_z[i] -= q * fz_factor * fz;
    }
}

// ---- Main compute ----

void CoulombPME::compute(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    const Cell& cell,
    std::span<double> forces_x,
    std::span<double> forces_y,
    std::span<double> forces_z,
    double& energy_out)
{
    int nx = params_.nx, ny = params_.ny, nz = params_.nz;

    spread_charges(pos_x, pos_y, pos_z, cell);

    double e_recip = compute_reciprocal_energy_and_potential(cell);

    double norm = 1.0 / (static_cast<double>(nx) * ny * nz);
    for (size_t i = 0; i < static_cast<size_t>(nx) * ny * nz; ++i) {
        grid_[i] *= norm;
    }

    double fx_factor = static_cast<double>(nx) / cell.matrix[0];
    double fy_factor = static_cast<double>(ny) / cell.matrix[4];
    double fz_factor = static_cast<double>(nz) / cell.matrix[8];

    interpolate_forces(pos_x, pos_y, pos_z, cell,
                       forces_x, forces_y, forces_z,
                       fx_factor, fy_factor, fz_factor);

    double energy = e_recip + self_energy_;
    energy_out += energy;
}
