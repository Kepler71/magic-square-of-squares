from pathlib import Path
import json

ns = {'__name__': 'claude_symbol_functions'}
exec(Path('s2_fp_independence.py').read_text().split('out = {}')[0], ns)
p = 10007
basis = []
certificate = []
attempts = 0
while len(basis) < 34 and attempts < 1000:
    attempts += 1
    a, b, c, r = ns['sample_point'](p)
    rr = ns['row_bits'](p, a, b, c, r)
    if rr is None:
        continue
    mask, ncols = rr
    x = mask
    for v in basis:
        x = min(x, x ^ v)
    if not x:
        continue
    basis.append(x)
    assert all((r[i,j]**2-a-i*b-j*c) % p == 0 for i,j in r)
    certificate.append({'a': a, 'b': b, 'c': c,
                        'roots_lex': [r[i,j] for i in (-1,0,1) for j in (-1,0,1)],
                        'mask': mask})
assert len(basis) == 34
rows = [r['mask'] for r in certificate]
# Independent elementary matrix elimination in explicit 0/1 arrays.
A = [[(m >> j) & 1 for j in range(34)] for m in rows]
rank = 0
for j in range(34):
    k = next((i for i in range(rank,34) if A[i][j]), None)
    if k is None:
        continue
    A[rank], A[k] = A[k], A[rank]
    for i in range(34):
        if i != rank and A[i][j]:
            A[i] = [x ^ y for x,y in zip(A[i],A[rank])]
    rank += 1
assert rank == 34
out = {'prime': p, 'columns': 34, 'rank': rank, 'attempts': attempts,
       'column_order': '16 half coordinates x,x-t in grid.lines_of order; b,c,b+c,b-c,b+2c,b-2c,2b+c,2b-c; 9 roots lex; constant bit 1',
       'witnesses': certificate}
Path('independence_certificate.json').write_text(json.dumps(out, indent=2))
print(json.dumps({k:v for k,v in out.items() if k != 'witnesses'}, indent=2))
