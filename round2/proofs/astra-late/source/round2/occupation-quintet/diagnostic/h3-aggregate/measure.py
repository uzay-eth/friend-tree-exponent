from pathlib import Path

print('side cohort cut band count: physical base joint sibling upstream box_gap compressed; ppm/core')
for s in ('lower','upper'):
    rows=[r.split('\t') for r in Path(__file__).with_name(f'sample-{s}.tsv').read_text().splitlines()]
    for name in ('pooled','heldout'):
        v=[r for r in rows if r[0].endswith('numeric4-root-tree.tsv')==(name=='heldout')]
        cores=sum({r[0]:int(r[3]) for r in v}.values())
        assert cores==({'pooled':896887,'heldout':235709}[name])
        for k in range(3):
            for b in (37,65,100,179,360,1998,'all','66+'):
                group=[r for r in v if int(r[1])==k and (b=='all' or (int(r[2])>=100 if b=='66+' else int(r[2])==b))]
                n=sum(int(r[4]) for r in group)
                t=[sum(float(r[i]) for r in group)*1e6/cores for i in range(5,12)]
                assert abs(t[0]-t[1]-t[2])<1e-10 and abs(t[2]-t[3]-t[4])<1e-10
                assert -1e-10<=t[6]<=t[2]+1e-10<=t[5]+2e-10
                print(s,name,k,b,n,*(f'{x:.8f}' for x in t))
