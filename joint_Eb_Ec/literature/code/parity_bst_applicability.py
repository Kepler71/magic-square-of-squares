# Контроль: применима ли теорема 1.1 Bremner–Silverman–Tzanakis (n бесквадратное, ЦЕЛЫЕ точки)
# к примитивному целому магическому квадрату из квадратов {a+ib+jc}.
# Утверждение: в примитивном квадрате все клетки нечётны, а b, c, b+c, b-c делятся на 8,
# поэтому ни один шаг не бесквадратный; а при переходе к бесквадратной модели x/t^2 (t чётно, a нечётно)
# целочисленность пропадает.
# Проверка: полный перебор a,b,c по модулю 2^k: все девять a+ib+jc — квадраты по модулю 2^k.
import itertools
for k in (4, 5, 6, 7):
    M = 2**k
    sq = {(x*x) % M for x in range(M)}
    mixed = 0; odd_bad_step = 0; odd_total = 0; even_total = 0
    for a, b, c in itertools.product(range(M), repeat=3):
        cells = [(a + i*b + j*c) % M for i in (-1, 0, 1) for j in (-1, 0, 1)]
        if not all(x in sq for x in cells):
            continue
        par = {x % 2 for x in cells}
        if len(par) == 2:
            mixed += 1
        elif par == {1}:
            odd_total += 1
            if any(s % 8 for s in (b, c, (b + c) % M, (b - c) % M)):
                odd_bad_step += 1
        else:
            even_total += 1
            # все клетки чётны => каждая ≡ 0 mod 4 (чётный квадрат), квадрат не примитивен
            assert all(x % 4 == 0 for x in cells) or M < 4
    print(f"mod 2^{k}: наборов (a,b,c) со всеми 9 квадратами: нечётных={odd_total}, чётных={even_total}, "
          f"смешанной чётности={mixed}, нечётных с шагом не ≡0 mod 8={odd_bad_step}")
print("Вывод, если mixed=0 и odd_bad_step=0 при всех k: примитивный квадрат имеет все клетки нечётными и 8 | b, c, b±c.")
