from fractions import Fraction as F

# Analytic witness calculations discussed in REPORT.md, and a compact exact summary.
def read():
    with open('exact.log') as f:
        for line in f:
            a = line.split()
            if a[0]=='m':
                n=int(a[1]); dp=int(a[7]); dm=int(a[9]); print('size',n,'histories',a[3])
            elif a[0]=='rooted': print(*a)
            elif '-' in a[0]:
                u,v=map(int,a[0].split('-'))
                if u+1==v: print(a[0],*[f'{float(F(int(x),d)):.10f}' for x,d in zip(a[1:],[dp,dp,dm,dm])])
            else:
                print(a[0],*[f'{float(F(int(x),d)):.10f}' for x,d in zip(a[1:],[dp,dp,dm,dm])])

if __name__=='__main__': read()
