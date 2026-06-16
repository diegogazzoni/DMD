#include "coulomb_pme.h"
#include "core/cell.h"
#include <cmath>

CoulombPME::CoulombPME(PMEGridParams params) : params_(params) {}

std::string CoulombPME::name() const { return "CoulombPME"; }

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
    // stub — FFTW integration pending
    double energy = 0.0;
    energy_out += energy;
}
