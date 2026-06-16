#pragma once

#include "force/force_component.h"
#include <vector>
#include <fftw3.h>

struct PMEParams {
    double cutoff{1.2};
    double ewald_coefficient{3.5};
    int nx{32}, ny{32}, nz{32};
    int spline_order{4};
};

class CoulombPME : public ForceComponent {
public:
    explicit CoulombPME(PMEParams params);
    ~CoulombPME() override;

    CoulombPME(const CoulombPME&) = delete;
    CoulombPME& operator=(const CoulombPME&) = delete;
    CoulombPME(CoulombPME&&) = default;
    CoulombPME& operator=(CoulombPME&&) = default;

    std::string name() const override;

    void set_charges(std::span<const double> charges);

    void compute(
        std::span<const double> pos_x,
        std::span<const double> pos_y,
        std::span<const double> pos_z,
        const Cell& cell,
        std::span<double> forces_x,
        std::span<double> forces_y,
        std::span<double> forces_z,
        double& energy_out
    ) override;

    static double theta4(double u);
    static double dtheta4(double u);
    static double calc_self_energy(double alpha, std::span<const double> charges, double k_e);

private:
    PMEParams params_;
    double k_e_{138.935456};
    std::vector<double> charges_;
    std::vector<double> grid_;
    fftw_plan plan_fwd_{nullptr};
    fftw_plan plan_bwd_{nullptr};

    double self_energy_;

    void spread_charges(
        std::span<const double> pos_x,
        std::span<const double> pos_y,
        std::span<const double> pos_z,
        const Cell& cell);
    double compute_reciprocal_energy_and_potential(const Cell& cell);
    void interpolate_forces(
        std::span<const double> pos_x,
        std::span<const double> pos_y,
        std::span<const double> pos_z,
        const Cell& cell,
        std::span<double> forces_x,
        std::span<double> forces_y,
        std::span<double> forces_z,
        double fx_factor,
        double fy_factor,
        double fz_factor);
};
