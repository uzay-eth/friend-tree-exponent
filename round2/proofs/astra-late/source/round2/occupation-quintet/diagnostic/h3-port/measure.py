from pathlib import Path

print('side cohort mother-band cores H3-rows ; direct sibling-box increment one-child-minorant Jensen-majorant (per core)')
for side in ('lower','upper'):
    rows=[x.split('\t') for x in Path(__file__).with_name(f'sample-{side}.tsv').read_text().splitlines()]
    assert len(rows)==18
    for cohort in ('pooled','heldout'):
        rr=[r for r in rows if r[0].endswith('/numeric4-root-tree.tsv')==(cohort=='heldout')]
        assert len(rr)==(16 if cohort=='pooled' else 2)
        for band,key in (('19–65',None),('19–37','0'),('38–65','1')):
            data=[r for r in rr if key is None or r[1]==key]
            cores=sum({r[0]:int(r[2]) for r in data}.values())
            n=sum(int(r[3]) for r in data)
            v=[sum(float(r[i]) for r in data)/cores for i in range(4,9)]
            assert abs(v[0]-v[1]-v[2])<1e-17 and 0<=v[3]<=v[2]<=v[4]
            print(side,cohort,band,cores,n,' '.join(f'{x:.12g}' for x in v))
