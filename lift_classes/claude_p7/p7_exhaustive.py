# Claude, 26.09: независимая проверка леммы «при p = 7 восемь АП-произведений квадраты в Q_7 ⇒ все девять клеток квадраты в Q_7».
# Полный перебор A = rz, B = sz по модулю 7^M: целый случай (z ∈ Z_7) и полюс v(z) = −n (A = a/7^n, B = b/7^n, не оба a, b ≡ 0 mod 7).
# Класс клетки в Q_7: (v mod 2, χ(unit mod 7)); квадрат ⇔ v чётно и unit — вычет. Если клетка ≡ 0 mod 7^M — «не определено» (считаем отдельно).
import sys
p = int(sys.argv[1]) if len(sys.argv) > 1 else 7
M = int(sys.argv[2]) if len(sys.argv) > 2 else 3
Q = p**M
QR = {x * x % p for x in range(1, p)}
def cls(num, shift):          # значение num / p^shift, num по модулю p^M
    num %= Q
    if num == 0: return None
    v = 0
    while num % p == 0: num //= p; v += 1
    return ((v - shift) % 2, 0 if (num % p) in QR else 1)   # (чётность оценки, невычет?)
cells = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
# 8 линий АРИФМЕТИЧЕСКОЙ сетки: строки i, столбцы j, две диагонали
lines = [[(i, j) for j in (-1, 0, 1)] for i in (-1, 0, 1)] + [[(i, j) for i in (-1, 0, 1)] for j in (-1, 0, 1)] \
      + [[(-1, -1), (0, 0), (1, 1)], [(-1, 1), (0, 0), (1, -1)]]
stats = dict(total=0, undetermined=0, eight_ok=0, all9_squares=0, counterexamples=0)
examples = []
def check(num_of, shift):
    c = {}
    for (i, j) in cells:
        k = cls(num_of(i, j), shift)
        if k is None: stats['undetermined'] += 1; return
        c[(i, j)] = k
    stats['total'] += 1
    for L in lines:
        pv = sum(c[x][0] for x in L) % 2; pc = sum(c[x][1] for x in L) % 2
        if pv or pc: return
    stats['eight_ok'] += 1
    if all(c[x] == (0, 0) for x in cells): stats['all9_squares'] += 1
    else:
        stats['counterexamples'] += 1
        if len(examples) < 5: examples.append((shift, c))
# целый случай: z ∈ Z_p  (A, B любые по модулю p^M)
for A in range(Q):
    for B in range(Q):
        check(lambda i, j: 1 + i * A + j * B, 0)
print('целый случай:', dict(stats))
# полюс v(z) = −n: клетка = (p^n + i a + j b) / p^n, a, b не оба делятся на p
for n in (1, 2):
    for a in range(Q):
        for b in range(Q):
            if a % p == 0 and b % p == 0: continue
            check(lambda i, j: p**n + i * a + j * b, n)
    print(f'после полюса n={n}:', dict(stats))
print('примеры контрпримеров:', examples)
