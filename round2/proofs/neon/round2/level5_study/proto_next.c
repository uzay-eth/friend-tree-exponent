/* Exact L4 -> L5 clipped-rate identity, exercised in double precision.
 * This fixture enumerates physical 5th neighbors, including D5-1 -> D5.
 * It is not a coupled L5 table solver. See REPORT.md for the row equations. */
#include <assert.h>
#include <math.h>
#include <stdio.h>

static double excess(int d, int d5) {
    return d < d5 ? 1.0/d - 1.0/d5 : 0.0;
}

static double kept(int d, int d5, int eligible) {
    return 1.0/d - (eligible ? excess(d, d5) : 0.0);
}

int main(void) {
    unsigned cases = 0;
    double max_error = 0.0;
    /* All other L4 slots qualify unless the chosen maximum degree crosses D5. */
    for (int d5 = 6; d5 <= 8; ++d5)
      for (int max_d4 = 2; max_d4 <= 20; ++max_d4)
        for (int da = 2; da <= 120; ++da)
          for (int db = 2; db <= 120; ++db)
            for (int delta_code = -3; delta_code <= 3; ++delta_code) {
                double delta = delta_code * 0.01823;
                int eligible = max_d4 < d5;
                double ka = kept(da, d5, eligible), kb = kept(db, d5, eligible);
                if (eligible) {
                    assert(ka <= 1.0/d5 && kb <= 1.0/d5);
                    assert(ka >= 0.0 && kb >= 0.0);
                }
                double direct = (1.0/da + 1.0/db) * delta;
                double split = (ka + kb) * delta +
                    (eligible ? excess(da, d5)*delta + excess(db, d5)*delta : 0.0);
                double err = fabs(direct - split);
                assert(err < 1e-15);
                max_error = fmax(max_error, err);
                ++cases;
            }
    printf("clipped-rate identity: %u cases; max_abs_error=%.17g\n", cases, max_error);
    return 0;
}
