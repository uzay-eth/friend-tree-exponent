"""Run the occupation proposals and replays on one core within the memory budget."""

import os
import resource


def resources():
    cpus=os.sched_getaffinity(0)
    cpu=int(os.environ.get('MOMENTS_CPU',13 if 13 in cpus else min(cpus)))
    os.sched_setaffinity(0,{cpu})
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[key]='1'
    soft,hard=resource.getrlimit(resource.RLIMIT_AS)
    limits=[6<<30]+[v for v in (soft,hard) if v!=resource.RLIM_INFINITY]
    resource.setrlimit(resource.RLIMIT_AS,(min(limits),hard))
