"""Compatibility with the pinned inner/tail row extractors and endpoint sum."""
from argparse import ArgumentParser
from pathlib import Path
from fractions import Fraction as Q
from rows import Rows, Cert, SPEC, embedding


def check(blob):
    c=Cert(blob)
    e=Rows(c)
    cases=[('hV',((4,3,1),(1,2,1),(1,1,1))),
           ('Yt',((3,3,1),(1,2,1),(0,1,1),(1,1,1))),
           ('P4',((3,3,1),(0,1,1),(0,1,1),(1,1,1)))]
    count=0
    for family,x in cases:
        a,ch,p=e.affine(family,x)
        old,free=e.row(family,x)
        assert a==old/c.S
        for j,parent,_,coeff in free:
            assert (p if parent else ch[j])==Q(coeff,c.S)
        for d in (18,20,32,64,65,90):
            z=((d-6,5,1),)+x[1:]
            a,g,p=e.controls(family,z)
            assert a+sum(min(v for _,v in q) for q in g)+min(0,p/2)>=0
            assert next(ids for f,ids in embedding(SPEC[family]) if f==family)==tuple(range(len(x)))
            count+=1
    print(blob.name,'inner identity and tail signed corner controls',count)


if __name__=='__main__':
    p=ArgumentParser();p.add_argument('blobs',nargs='+',type=Path)
    for b in p.parse_args().blobs:check(b)
    paid=(Q('3.06204334292295e-9')+Q('1.49696294463021e-12')+Q('3.16565034087608e-11'),
          Q('4.00077770896578e-9')+Q('2.90279385539250e-12')+Q('2.18012978802676e-11'))
    old=(Q('.565676634491'),Q('.565771874464'))
    den=(Q('1.356833021311'),Q('1.367169995210'))
    vals=[]
    for side,b,d,r in zip((1,-1),old,den,paid):
        n=side*((side*(b+side*r/d)*10**12)//1)
        vals.append(f'{n//10**12}.{n%10**12:012d}')
    assert vals==['0.565676636772','0.565771871520']
    print('combined outward',*vals)
