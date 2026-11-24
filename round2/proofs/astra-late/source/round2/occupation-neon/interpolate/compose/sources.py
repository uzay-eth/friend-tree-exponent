"""Frozen, hash-checked rational operands for the repaired corrector."""
import ast
from fractions import Fraction as Q
import hashlib
from pathlib import Path

OCC=Path(__file__).resolve().parents[2]
PINS={
    'renewal/check.py':'423cfdddf34559c705b55b7f50f8c4a57f655915226d23cf7170028a0ff81388',
    'rows/charged-f.json':'415a712878e03b71c431a36bf4cf8eb11a047c2a47cafb5d9298fe23bc1561da',
}


def authenticated(name):
    data=(OCC/name).read_bytes()
    if hashlib.sha256(data).hexdigest()!=PINS[name]:raise ValueError(f'wrong frozen source: {name}')
    return data


def bands():
    src=ast.parse(authenticated('renewal/check.py'))
    a,=(node.value for node in src.body if isinstance(node,ast.Assign) and
        any(isinstance(t,ast.Name) and t.id=='BANDS' for t in node.targets))
    out={s:tuple(map(Q,v)) for s,v in ast.literal_eval(a).items()}
    assert out.keys()=={'lower','upper'} and all(len(v)==5 and min(v)>0 for v in out.values())
    return out
