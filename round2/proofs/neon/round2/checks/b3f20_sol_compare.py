"""Parameterize the earlier independent row-major comparison sweep for D0=20."""
from pathlib import Path
import sys

s = Path(__file__).with_name('b3f18_sol_comparison.c').read_text()
changes = {
    '#include "../cert/cert_b3f_format.h"': '#include "round2/cert/cert_b3f_format.h"',
    'N=210, NG=435, P=397, HU=187, CU=352, P4=35, YT=20, BU=60, BW=35, S3=24':
        'N=253, NG=523, P=481, HU=228, CU=432, P4=54, YT=35, BU=96, BW=54, S3=45',
    'c.D0!=18||c.KT!=17||c.LF!=14000||c.KF!=180':
        'c.D0!=20||c.KT!=19||c.LF!=16000||c.KF!=200',
    'Ext tail[19],ta[2][18];': 'Ext tail[21],ta[2][20];',
    'k<=18;k++': 'k<=c.KT+1;k++',
    'k<=17': 'k<=c.KT',
    'k<=181;k++': 'k<=c.KF+1;k++',
    'k==181?0:14000': 'k>c.KF?0:c.LF',
    'j=D<18?': 'j=D<c.D0?',
    'tail[k<18?k:18]': 'tail[k<=c.KT?k:c.KT+1]',
    'k<=6&&D<9': 'k<=c.KYT&&D<c.DYT',
    'k<=6&&D<8': 'k<=c.KS3&&D<c.DS3',
    'v+(p*YT+i)*np': 'v+((size_t)p*YT+i)*np',
}
for a, b in changes.items():
    assert a in s
    s = s.replace(a, b)
sys.stdout.write(s)
