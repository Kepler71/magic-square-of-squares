load('/home/kep/magicKube/descent/n61/independent/lib.sage')
import random
random.seed(int(2026))
k = kk; r = rr
fails = []
def ck(name, cond, info=''):
    if not cond:
        fails.append((name, info)); print("FAIL", name, info)

print("=== A. is_square_P2: сверка с p-адическим полем Sage (независимый путь) ===")
Rx = PolynomialRing(Qp(2, 120), 'X'); X = Rx.gen()
K2 = Qp(2, 120).extension(X**2 - X - 41, names='w'); w2 = K2.gen()
def to_K2(x):
    c0, c1 = coords_om(x)
    return K2(QQ(c0)) + K2(QQ(c1)) * w2
bad = 0
for _ in range(300):
    a0 = QQ(random.randint(-60, 60)); a1 = QQ(random.randint(-60, 60))
    x = a0 + a1 * r
    if x == 0: continue
    mine = is_square_P2(x)
    theirs = to_K2(x).is_square()
    if mine != bool(theirs): bad += 1; print("  mismatch", x, mine, theirs)
ck("is_square_P2 vs Qp", bad == 0, f"bad={bad}")
# статистика: доля квадратов среди единиц ~ 1/8
cnt = 0; tot = 0
for _ in range(400):
    x = QQ(random.randint(-99, 99)) + QQ(random.randint(-99, 99)) * r
    if x == 0 or val2(x) != 0: continue
    tot += 1; cnt += 1 if is_square_P2(x) else 0
print(f"  доля квадратов среди 2-единиц: {cnt}/{tot} (ожидание ~1/8)")
ck("square density P2", tot > 0 and abs(cnt/tot - 1/8) < 0.06, f"{cnt}/{tot}")
# квадраты обязаны распознаваться
ck("squares are squares P2", all(is_square_P2((QQ(random.randint(1,40)) + QQ(random.randint(-40,40))*r)**2) for _ in range(50)))

print("=== B. is_square_odd: квадраты и плотность ===")
oddP = [P for q in [3, 5, 7, 11, 13, 31, 61, 337] for P in k.primes_above(q)]
bad = 0
for P in oddP:
    for _ in range(30):
        t = QQ(random.randint(-50, 50)) + QQ(random.randint(-50, 50)) * r
        if t == 0 or t.valuation(P) == Infinity: continue
        if not is_square_odd(P, t**2): bad += 1
ck("squares are squares odd", bad == 0, f"bad={bad}")

print("=== C. символ Гильберта: тождества ===")
places_all = [('inf', 1), ('inf', -1)] + [P for q in primes(2, 70) for P in k.primes_above(q)]
def rnd():
    while True:
        x = QQ(random.randint(-40, 40)) + QQ(random.randint(-40, 40)) * r
        if x != 0: return x
bad1 = bad2 = bad3 = bad4 = 0
for _ in range(25):
    a = rnd(); b = rnd(); c = rnd()
    for pl in places_all[:14]:
        try:
            if hsym(pl, a, -a) != 1: bad1 += 1; print("  (a,-a)", pl, a)
            if a != 1 and hsym(pl, a, 1 - a) != 1: bad2 += 1; print("  (a,1-a)", pl, a)
            if hsym(pl, a, b) != hsym(pl, b, a): bad3 += 1; print("  sym", pl, a, b)
            if hsym(pl, a * b, c) != hsym(pl, a, c) * hsym(pl, b, c): bad4 += 1; print("  bilin", pl, a, b, c)
        except Exception as ex:
            print("  EXC", pl, ex); bad1 += 1
ck("(a,-a)=1", bad1 == 0); ck("(a,1-a)=1", bad2 == 0)
ck("симметрия", bad3 == 0); ck("билинейность", bad4 == 0)

print("=== D. формула произведения (глобальный контроль) ===")
badp = 0
for _ in range(20):
    a = rnd(); b = rnd()
    S = set()
    for x in (a, b, k(2)):
        for P, e_ in k.ideal(x).factor(): S.add(P)
    prod_ = 1
    for P in S: prod_ *= hsym(P, a, b)
    for sr in (1, -1): prod_ *= hilbert_real(a, b, sr)
    if prod_ != 1: badp += 1; print("  product != 1:", a, b)
ck("формула произведения", badp == 0, f"bad={badp}")

print("=== E. невырожденность и кэш при P|2 ===")
basis, G = p2_basis()
print("  базис k_P2^*/(k_P2^*)^2:", basis)
print("  Грам:", G)
MG = matrix(GF(2), [[0 if G[i][j] == 1 else 1 for j in range(4)] for i in range(4)])
ck("Gram невырождена", MG.rank() == 4, str(MG.rank()))
# быстрый hilbert_2 == прямой поиск
bad = 0
for _ in range(30):
    a = rnd(); b = rnd()
    if hilbert_2(a, b) != hilbert_2_search(a, b): bad += 1; print("  fast!=search", a, b)
ck("hilbert_2 == прямой поиск", bad == 0)
# нечётные единицы: символ = 1
bad = 0
for P in oddP:
    for _ in range(20):
        a = rnd(); b = rnd()
        if a.valuation(P) == 0 and b.valuation(P) == 0 and hilbert_odd(P, a, b) != 1: bad += 1
ck("нечётные единицы -> 1", bad == 0)

print("\nИТОГ тестов:", "ВСЕ ОК" if not fails else f"ОШИБКИ: {fails}")
