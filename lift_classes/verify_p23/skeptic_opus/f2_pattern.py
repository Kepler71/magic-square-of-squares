# Независимая проверка шага 1 (линейная алгебра над F2) — скептик Claude (Opus), 2026-09-26.
# Клетки (i,j), i,j in {-1,0,1}; значение 1+(i r + j s) z. Центр (0,0) = 1, класс 0.
# Восемь линий АРИФМЕТИЧЕСКОЙ сетки: 3 строки, 3 столбца, 2 диагонали.
# Перебираем ВСЕ 2^8 векторов классов (по одному биту — достаточно, т.к. пространство
# Q*/Q*^2 — F2-векторное и уравнения линейны покоординатно).
from itertools import product
cells=[(i,j) for i in (-1,0,1) for j in (-1,0,1) if (i,j)!=(0,0)]
lines=[[(i,-1),(i,0),(i,1)] for i in (-1,0,1)]+[[(-1,j),(0,j),(1,j)] for j in (-1,0,1)]
lines+=[[(-1,-1),(0,0),(1,1)],[(-1,1),(0,0),(1,-1)]]
assert len(lines)==8
T=[(-1,0),(1,1),(0,-1)]   # (1-rz)(1+(r+s)z)(1-sz)
L=[(-1,0),(1,-1),(0,1)]   # (1-rz)(1+(r-s)z)(1+sz)
sols=[]
for bits in product((0,1),repeat=8):
    x=dict(zip(cells,bits)); x[(0,0)]=0
    if all(sum(x[c] for c in ln)%2==0 for ln in lines):
        sols.append(x)
print('число решений над F2:',len(sols),'(ожидается 4 = размерность 2)')
weights=sorted(sum(x[c] for c in cells) for x in sols)
print('веса векторов:',weights,'(ожидается [0,4,6,6])')
ok=True
for x in sols:
    S=x[(1,1)]; U=x[(0,1)]
    pat={(-1,-1):S,(-1,0):U,(-1,1):S^U,(0,-1):U,(0,0):0,(0,1):U,(1,-1):S^U,(1,0):U,(1,1):S}
    ok&= all(x[c]==pat[c] for c in pat)
    ok&= all(x[(i,j)]==x[(-i,-j)] for (i,j) in pat)          # противоположные клетки равны
    ok&= (sum(x[c] for c in T)%2==S) and (sum(x[c] for c in L)%2==(S^U))
    # ни один ненулевой вектор не сидит на <=3 клетках
    w=sum(x[c] for c in cells); ok&= (w==0 or w>=4)
print('матрица S,U,S+U; [T]=S; [L]=S+U; противоположные равны; вес 0 или >=4:',ok)
# Какие подмножества линий уже дают те же 4 решения (для справки: насколько жёстко используются 8 линий)
from itertools import combinations
minimal=[]
for k in range(1,9):
    for sub in combinations(range(8),k):
        cnt=0
        for bits in product((0,1),repeat=8):
            x=dict(zip(cells,bits)); x[(0,0)]=0
            if all(sum(x[c] for c in lines[t])%2==0 for t in sub): cnt+=1
        if cnt==4: minimal.append(sub)
    if minimal: break
print('минимальное число линий, дающих ту же картину:',len(minimal[0]),'; примеров:',len(minimal))
print('(восемь уравнений имеют ранг 6 — два уравнения избыточны)')
