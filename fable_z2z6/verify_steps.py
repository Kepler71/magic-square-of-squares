# Численная проверка арифметических шагов ручных доказательств (все взаимно простые 1<=a<b<=200)
from math import gcd
def v(n,p):
    n=abs(n); e=0
    while n%p==0: n//=p; e+=1
    return e
bad = 0; n = 0
for b in range(2, 201):
    for a in range(1, b):
        if gcd(a,b) != 1: continue
        n += 1
        M = b**3*(2*a+b); N = a**3*(a+2*b)
        alpha = a*a+a*b+b*b; beta = a*a*b*b*(a+b)**2
        # (i) p>=5 делит не более одной формы
        forms = [a, b, a-b, a+b, a+2*b, 2*a+b]
        for p in range(5, 700):
            if all(p % q for q in range(2, int(p**0.5)+1)):
                if sum(1 for F in forms if F % p == 0) > 1: bad += 1; print("two forms", a, b, p)
        # (ii) p=3: ровно одна из групп ab(a+b) / (a-b)(a+2b)(2a+b) делится на 3; в типе 2: v3(alpha)=1, beta≡4 mod 9, 4ab(a+b)≡(2a)^3 mod 9,
        t1 = (a*b*(a+b)) % 3 == 0; t2 = ((a-b)*(a+2*b)*(2*a+b)) % 3 == 0
        if t1 == t2: bad += 1; print("3-type", a, b)
        if t2:
            if v(alpha,3) != 1 or beta % 9 != 4 or (4*a*b*(a+b) - (2*a)**3) % 9 != 0: bad += 1; print("type2 mod9", a, b)
            f = a**6 + 3*a**5*b - 5*a**3*b**3 + 3*a*b**5 + b**6
            if v(f,3) != 1: bad += 1; print("v3(f)", a, b, v(f,3))
            # все три точки порядка 2 сравнимы с -1 mod 3 (сдвинутая модель): -a^2b^2, -b^2c^2, -a^2c^2
            c = a+b
            if any((e+1) % 3 for e in (-a*a*b*b, -b*b*c*c, -a*a*c*c)): bad += 1; print("2-tors cusp", a, b)
        else:
            # тип 1 при 3: M или -M — квадрат единицы mod 3 (расщеплённый узел) — проверяем, что узловая касательная рациональна
            if a % 3 == 0 or b % 3 == 0:
                if (M % 3 not in (1,) and a % 3 == 0) or (N % 3 not in (1,) and b % 3 == 0): bad += 1; print("split3", a, b)
        # (iii) p=2: модель Фрея. Случай 2|a: A=-M, B=N; 2|b: A=-N, B=M; a,b нечётные: A=M, B=N-M
        if a % 2 == 0: A, B = -M, N; va = v(a,2)
        elif b % 2 == 0: A, B = -N, M; va = v(b,2)
        else: A, B = M, N-M; va = v(a+b,2)
        if (B - A - 1) % 4 != 0 or B % 16 != 0: bad += 1; print("frey", a, b)
        cpar = ((B - A - 1)//4) % 2      # 0 => расщеплённая, 1 => нерасщеплённая
        type1 = (va >= 2)
        if type1 != (cpar == 0): bad += 1; print("split2", a, b, va, cpar)
        # T в новой модели: X = x1/4, x1 = a^2b^2 (+M если сдвиг для нечётных; для 2|b сдвиг на N)
        if a % 2 == 0: x1 = a*a*b*b
        elif b % 2 == 0: x1 = a*a*b*b
        else: x1 = a*a*b*b + M
        X = x1 // 4 if x1 % 4 == 0 else None
        y1 = a*a*b*b*(a+b)**2
        if type1:
            if X is None or X % 2 != 0 or (y1 - 4*X) % 16 != 0: bad += 1; print("T node", a, b)
        else:
            if x1 % 4 != 0 or (x1//4) % 2 != 1: bad += 1; print("T smooth", a, b)
print("pairs", n, "bad", bad)
