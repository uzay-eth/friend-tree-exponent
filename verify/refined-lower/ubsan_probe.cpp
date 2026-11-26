#include <cstdint>
#include <cstdio>
#include <cstdlib>
int main(int c,char**v){ int64_t a=INT64_MAX - (c-1); __int128 b=((__int128)1<<126); b = b + b*(c); a = a + atoi(v[1]); printf("%lld %d\n",(long long)a,(int)(b>0)); }
