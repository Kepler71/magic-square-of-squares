"""Независимая проверка четырёх 2-адических лемм Codex (без его check_mod128.py).
Перебираем ОСТАТКИ (m, n, d) mod 128, удовлетворяющие условиям леммы (условия d | mn и взаимной простоты не накладываем —
это лишь расширяет перебор и для доказательства отсутствия точек безопасно). Для каждого набора коэффициентов
квартики W² = c4 U⁴ + c2 U²V² + c0 V⁴ (mod 128) проверяем все примитивные (U,V): (U,1) и (1,V) с чётным V.
Если ни при одной паре W² ≡ RHS (mod 128) не решается — точек над Q_2 нет (редукция 2-адической точки)."""
M=128
SQ=set((w*w)%M for w in range(M))
def v2(x):
    if x%M==0: return 99
    e=0
    while x%2==0: x//=2; e+=1
    return e
def inv(x): return pow(x,-1,M)
def empty(c4,c2,c0):
    for U in range(M):
        if (c4*pow(U,4,M)+c2*U*U+c0)%M in SQ: return False
    for V in range(0,M,2):
        if (c4+c2*V*V+c0*pow(V,4,M))%M in SQ: return False
    return True
def run(name,triples,control=None):
    tr=set(triples); bad=[t for t in tr if not empty(*t)]
    print(f"{name}: наборов коэффициентов {len(tr)}; с возможной 2-адической точкой: {len(bad)} {sorted(bad)[:4]}")
    if control:
        ctr=set(control); ok=sum(1 for t in ctr if not empty(*t))
        print(f"   положительный контроль (условие леммы нарушено): наборов {len(ctr)}, с возможной точкой {ok} (должно быть > 0)")
L1=[];C1=[]
for m in range(1,M,2):
    for n in range(1,M,2):
        H=((m-n)%8==0 or (m+n)%8==0)
        A=(m*m+n*n)%M; B=(m*m*n*n)%M
        for d in range(1,M,2):
            t=((d)%M,A,(B*inv(d))%M)
            if H and d%8 in (3,5): L1.append(t)
            elif H and d%8 in (1,7): C1.append(t)
run("Лемма 1 (C_d, H, d≡3,5 mod 8)",L1,C1)
L2=[];C2=[]
for m in range(4,M,8):                      # v2(m)=2
    for n in range(1,M,2):
        A=(m*m+n*n)%M
        mn2=(m*m*n*n)%(2*M)//2               # m²n²/2 mod 128 (m²n² mod 256 определено остатками mod 128)
        for dl in range(1,M,2):              # d = 2δ
            L2.append(((2*dl)%M,A,(mn2*inv(dl))%M))
        for d in range(1,M,2):
            C2.append((d,A,(m*m*n*n*inv(d))%M))
run("Лемма 2 (C_d, v2(mn)=2, d чётно)",L2,C2)
L3=[];C3=[]
for m in range(0,M,2):
    e=v2(m) if m else 99
    for n in range(1,M,2):
        A=(m*m+n*n)%M; D=(n*n-m*m)%(2*M); D2=(D*D)%M
        for d in range(1,M,2):
            t=(d,(-2*A)%M,(D2*inv(d))%M)
            if e!=2 and d%8==5: L3.append(t)
            elif e==2 and d%8==5: C3.append(t)
run("Лемма 3 (C'_d, один чётный, v2≠2, d≡5 mod 8)",L3,C3)
L4=[];C4=[]
for m in range(1,M,2):
    for n in range(1,M,2):
        H=((m-n)%8==0 or (m+n)%8==0)
        A=(m*m+n*n)%M; D=(n*n-m*m)%(4*M)
        D2half=((D*D)%(2*M))//2               # D²/2 mod 128
        for dl in range(1,M,2):
            t=((2*dl)%M,(-2*A)%M,(D2half*inv(dl))%M)
            (C4 if H else L4).append(t)
run("Лемма 4 (C'_d, m,n нечётны, не H, d чётно)",L4,C4)
