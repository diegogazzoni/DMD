import time
from _dmd_core import Engine


class SimulationStatus:
    """Tracks and prints simulation progress (GROMACS-style)."""

    def __init__(self, engine: Engine, dt: float = 0.002, interval: int = 100):
        self.engine = engine
        self.dt = dt
        self.interval = interval
        self.t0 = time.time()
        self.step0 = 0

    def start(self):
        self.t0 = time.time()
        self.step0 = self.engine.current_step
        self._print_header()

    def step(self):
        step = self.engine.current_step
        if step % self.interval == 0:
            self._print(step)

    def finish(self):
        elapsed = time.time() - self.t0
        steps = self.engine.current_step - self.step0
        ns_day = _ns_per_day(steps, elapsed, self.dt)
        print(f"\r  [{steps} steps in {elapsed:.1f}s]  ns/day={ns_day:.1f}")
        print()

    def _print_header(self):
        cols = f"{'Step':>10} {'Time (ps)':>10} {'PE (kJ/mol)':>14} {'Box (nm)':>10} {'ns/day':>8}"
        print(cols)
        print("-" * len(cols))

    def _print(self, step: int):
        elapsed = time.time() - self.t0
        ns_day = _ns_per_day(step - self.step0, elapsed, self.dt)
        print(
            f"{step:>10} {self.engine.current_time:>10.4f} "
            f"{self.engine.potential_energy:>14.4f} "
            f"{self.engine.box_size:>10.4f} "
            f"{ns_day:>8.1f}"
        )


def _ns_per_day(steps: int, elapsed_s: float, dt_ps: float) -> float:
    if elapsed_s <= 0 or steps <= 0:
        return 0.0
    return dt_ps * steps / elapsed_s * 86.4
