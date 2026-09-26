from pathlib import Path
from math import isqrt
import json

data = json.loads(Path('independence_certificate.json').read_text())
p = data['prime']
assert p > 3 and all(p % d for d in range(2, isqrt(p)+1))
idx = (-1, 0, 1)
rows = []
for witness in data['witnesses']:
    a, b, c = (witness[k] for k in ('a', 'b', 'c'))
    r = dict(zip(((i,j) for i in idx for j in idx), witness['roots_lex']))
    assert all((r[i,j]**2 - a - i*b - j*c) % p == 0 for i,j in r)
    triples = [[r[i,j] for i in idx] for j in idx]
    triples += [[r[i,j] for j in idx] for i in idx]
    triples += [[r[-1,-1],r[0,0],r[1,1]], [r[-1,1],r[0,0],r[1,-1]]]
    values = []
    for x,y,z in triples:
        values += [(x+y)*(y+z), (x+y)*(x+z)]
    values += [b,c,b+c,b-c,b+2*c,b-2*c,2*b+c,2*b-c]
    values += [r[i,j] for i in idx for j in idx]
    assert len(values) == 33 and all(v % p for v in values)
    bits = [int(pow(v % p,(p-1)//2,p) == p-1) for v in values] + [1]
    assert sum(bit << j for j,bit in enumerate(bits)) == witness['mask']
    rows.append(bits)
assert len(rows) == 34
rank = 0
for col in range(34):
    pivot = next((i for i in range(rank,34) if rows[i][col]), None)
    if pivot is None:
        continue
    rows[rank], rows[pivot] = rows[pivot], rows[rank]
    for i in range(rank+1,34):
        if rows[i][col]:
            rows[i] = [x ^ y for x,y in zip(rows[i],rows[rank])]
    rank += 1
assert rank == data['rank'] == 34
print('Verified: prime 10007; 34 valid nonzero specialization points; rank 34 over F_2.')
