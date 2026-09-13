from math import gcd
from odd_criterion import odd_bound
import json
def closes(m,n):
    a,b,_=odd_bound(min(m*m,n*n),max(m*m,n*n)); return a*b<=4
known=[(m,n) for m in range(2,21) for n in range(1,m) if gcd(m,n)==1]
proved20=[p for p in known if closes(*p)]
print(f"G1, m,n ≤ 20: пар {len(known)}; нечётный критерий доказывает rank B = 0 у {len(proved20)} (PARI: 62, критерий Fable с 2-адикой: 44)")
for N in (50,100,200):
    pairs=[(m,n) for m in range(2,N+1) for n in range(1,m) if gcd(m,n)==1]
    c=sum(closes(*p) for p in pairs)
    print(f"G1, m,n ≤ {N}: пар {len(pairs)}; закрыто доказанно {c} ({100*c/len(pairs):.1f}%)")
json.dump({'m_le_20':proved20},open('g1_odd.json','w'))
