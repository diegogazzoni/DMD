#include "integrate/integrator.h"
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

    double e0 = 0.5 * k * x0 * x0;

    for (int step = 0; step < 1000; ++step) {
        integrator.half_kick(sys, dt);
        integrator.advance(sys, dt);
        sys.forces_x[0] = -k * sys.pos_x[0];
        integrator.half_kick(sys, dt);

        double pe = 0.5 * k * sys.pos_x[0] * sys.pos_x[0];
        double ke = 0.5 * sys.vel_x[0] * sys.vel_x[0];

        if (step % 200 == 0) {
            EXPECT_NEAR(ke + pe, e0, 2e-3);
        }
    }

    double omega = std::sqrt(k);
    double t = 1000 * dt;
    double x_exact = x0 * std::cos(omega * t);
    EXPECT_NEAR(sys.pos_x[0], x_exact, 1e-2);
}
