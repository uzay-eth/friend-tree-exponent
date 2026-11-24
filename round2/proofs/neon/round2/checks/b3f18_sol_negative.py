"""Additional fail-closed checks of the actual KT offset and arithmetic guardrails."""
from pathlib import Path
import struct
import subprocess

base = Path('./smoke_lo.bin').read_bytes()
for name, offset, fmt, value, marker in (
    ('true-kt',60,'i',8,'invalid header (err=-7)'),
    ('huge-sr',40,'q',2**63-1,'1024:77: runtime error: signed integer overflow'),
    ('huge-kt',60,'i',2**31-1,'invalid header (err=-7)'),
    ('huge-lf',48,'i',2**31-1,'cert_b3f_format.h:175:'),
):
    path=Path(f'./sol-{name}.bin')
    data=bytearray(base)
    struct.pack_into('<'+fmt,data,offset,value)
    path.write_bytes(data)
    r=subprocess.run(['./sol-b3f-ubsan',str(path),'--threads','2','--tree-check','0'],
                     text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(f'{name}: rc={r.returncode}\n{r.stdout}',end='')
    assert r.returncode==1 and marker in r.stdout
