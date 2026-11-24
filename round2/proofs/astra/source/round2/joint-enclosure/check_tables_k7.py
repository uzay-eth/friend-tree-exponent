"""Run the independent Python table checks only for the unsigned K7/144 inputs."""

from pathlib import Path
import sys


if __name__ == "__main__":
    tree, side = Path(sys.argv[1]), sys.argv[2]
    sys.path.insert(0, str(tree/"round2/moments"))
    import check_tables as c

    c.resources()
    import numpy

    c.np = numpy
    t = c.Tables(side)
    print(f"{side}: independently checked {len(t.finite)} mbar entries and analytic tail; "
          f"profile={dict(t.inputs(7, 144))}, overlap={dict(t.overlap(7))}, "
          f"conditioned={t.conditioned(7)}, tail-lines={t.tail_lines(7, 144)}", flush=True)
