#include "integrate/integrator.h"
#include "core/cell.h"
#include "core/system_data.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(IntegratorTest, OneDimensionalHarmonic) {
    SystemData sys(1);
    sys.masses[0] = 1.0;

    double k = 100.0;
    double dt = 0.001;
    double x0 = 1.0;

    sys.pos_x[0] = x0;
    sys.vel_x[0] = 0.0;
    sys.forces_x[0] = -k * x0;

    Integrator integrator;
    Cell cell(10.0, 10.0, 10.0);

    double omega = std::sqrt(k);
    double e0 = 0.5 * k * x0 * x0;
    double t = 0.0;

    for (int step = 0; step < 1000; ++step) {
        integrator.step(sys, dt);

        sys.forces_x[0] = -k * sys.pos_x[0];
        t += dt;

        double ke = 0.5 * sys.masses[0] * sys.vel_x[0] * sys.vel_x[0];
        double pe = 0.5 * k * sys.pos_x[0] * sys.pos_x[0];

        if (step % 200 == 0) {
            EXPECT_NEAR(ke + pe, e0, 1e-6);
        }
    }

    double t_exact = omega * t;
    double x_exact = std::cos(t_exact);
    EXPECT_NEAR(sys.pos_x[0], x_exact, 1e-3);
}
