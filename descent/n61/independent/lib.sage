# Независимая реализация локальных квадратов и символов Гильберта над k = Q(sqrt(165)).
# НЕ используется k.hilbert_symbol и никакой код из ctp.sage / ek_descent.sage.
#
# k = Q(r), r^2 = 165, O_k = Z[om], om = (1+r)/2, om^2 = om + 41.
# 2 инертно: P2 = (2), e = 1, f = 2, F_4 = O/2.

import itertools

kk = QuadraticField(165, 'r')
rr = kk.gen()

# ---------- координаты в базисе (1, om) ----------
def coords_om(x):
    """x = c0 + c1*om, c0,c1 in QQ (om = (1+r)/2, r = 2om-1)."""
    a0, a1 = QQ(list(kk(x))[0]), QQ(list(kk(x))[1])   # x = a0 + a1*r
    return (a0 - a1, 2*a1)

def from_om(c0, c1):
    return kk(c0) + kk(c1) * (1 + rr)/2

# ---------- арифметика в O/2^m = (Z/2^m)[om]/(om^2-om-41) ----------
def red2(x, m):
    """x в k, 2-целый -> (c0,c1) mod 2^m."""
    c0, c1 = coords_om(x)
    M = 2**m
    out = []
    for c in (c0, c1):
        num, den = ZZ(c.numerator()), ZZ(c.denominator())
        if den % 2 == 0:
            raise ValueError("not 2-integral")
        out.append(int((num * inverse_mod(den, M)) % M))
    return (out[0], out[1])

def mul2(a, b, m):
    M = 2**m
    a0, a1 = a; b0, b1 = b
    # (a0+a1 w)(b0+b1 w) = a0b0 + (a0b1+a1b0) w + a1b1 w^2,  w^2 = w + 41
    c0 = (a0*b0 + 41*a1*b1) % M
    c1 = (a0*b1 + a1*b0 + a1*b1) % M
    return (c0, c1)

def val2(x):
    """v_{P2}(x) для x в k*, P2 = (2) инертно => v = 2-адическая оценка координат."""
    c0, c1 = coords_om(x)
    if c0 == 0 and c1 == 0:
        return Infinity
    vs = []
    for c in (c0, c1):
        if c == 0: continue
        vs.append(ZZ(c.numerator()).valuation(2) - ZZ(c.denominator()).valuation(2))
    return min(vs)

# множество квадратов единиц в O/2^3 (2e+1 = 3)
_SQU3 = None
def squares_mod8():
    global _SQU3
    if _SQU3 is None:
        S = set()
        for c0 in range(8):
            for c1 in range(8):
                # единица <=> ненулевой образ в F_4 = O/2 <=> (c0,c1) != (0,0) mod 2
                if c0 % 2 == 0 and c1 % 2 == 0: continue
                t = (c0, c1)
                S.add(mul2(t, t, 3))
        _SQU3 = S
    return _SQU3

def is_square_P2(x):
    """x в k*: является ли x квадратом в k_{P2} = Q_2(sqrt165)?"""
    x = kk(x)
    if x == 0: raise ValueError("zero")
    v = val2(x)
    if v % 2 != 0: return False
    w = x / kk(2)**v
    return red2(w, 3) in squares_mod8()

# ---------- локальные квадраты в нечётных P ----------
_ODD_CACHE = {}
def odd_data(P):
    key = str(P)
    if key not in _ODD_CACHE:
        F = P.residue_field()
        q = F.cardinality()
        pi = kk.uniformizer(P, others='positive')
        _ODD_CACHE[key] = (F, q, pi)
    return _ODD_CACHE[key]

def is_square_odd(P, x):
    x = kk(x)
    if x == 0: raise ValueError("zero")
    F, q, pi = odd_data(P)
    v = x.valuation(P)
    if v % 2 != 0: return False
    u = x / pi**v
    return F(u)**((q - 1) // 2) == F(1)

def chi(P, u):
    """квадратичный вычет u (P-единица) в F_q: +-1"""
    F, q, pi = odd_data(P)
    t = F(u)**((q - 1) // 2)
    return 1 if t == F(1) else -1

# ---------- символ Гильберта: нечётное P (ручная ручная формула) ----------
def hilbert_odd(P, a, b):
    """(a,b)_P = (-1)^{al*be*(q-1)/2} * chi(u)^be * chi(v)^al,  a = pi^al u, b = pi^be v."""
    a = kk(a); b = kk(b)
    if a == 0 or b == 0: raise ValueError("zero")
    F, q, pi = odd_data(P)
    al = ZZ(a.valuation(P)); be = ZZ(b.valuation(P))
    u = a / pi**al; v = b / pi**be
    s = 1
    if (al * be) % 2 == 1 and (q - 1) // 2 % 2 == 1:
        s = -1
    s *= chi(P, u)**int(be % 2)
    s *= chi(P, v)**int(al % 2)
    return s

# ---------- символ Гильберта: P | 2, прямой поиск нуля a X^2 + b Y^2 = Z^2 ----------
_M_SEARCH = 6      # y,z пробегают O/2^6 (нужно >= 4, см. отчёт)
_REPS2 = None
def reps_mod(m):
    """представители O/2^m"""
    M = 2**m
    return [from_om(c0, c1) for c0 in range(M) for c1 in range(M)]

def _norm_sq_class(x):
    """x * (квадрат) c v(x) в {0,1}"""
    v = val2(x)
    return x / kk(2)**(2 * (v // 2))

def hilbert_2_search(a, b):
    """(a,b)_{P2} = 1 <=> a X^2 + b Y^2 = Z^2 нетривиально разрешимо в k_{P2}."""
    a = kk(a); b = kk(b)
    if a == 0 or b == 0: raise ValueError("zero")
    if is_square_P2(a) or is_square_P2(b): return 1
    an = _norm_sq_class(a); bn = _norm_sq_class(b)
    m = _M_SEARCH
    M = 2**m
    # карта A: y = 1, z пробегает O/2^m
    for c0 in range(M):
        for c1 in range(M):
            z = from_om(c0, c1)
            t = z*z - bn
            if t == 0: continue
            if is_square_P2(t / an): return 1
    # карта B: z = 1, y в P2 (v>=1) mod 2^m
    for c0 in range(0, M, 2):
        for c1 in range(0, M, 2):
            y = from_om(c0, c1)
            if y == 0: continue
            t = 1 - bn * y * y
            if t == 0: continue
            if is_square_P2(t / an): return 1
    return -1

# ---- ускорение: билинейная форма на k_{P2}^*/(k_{P2}^*)^2 (dim 4) ----
_P2_BASIS = None
_P2_GRAM = None
def p2_basis():
    global _P2_BASIS, _P2_GRAM
    if _P2_BASIS is not None: return _P2_BASIS, _P2_GRAM
    cands = [kk(2), kk(-1), kk(5), kk(3), kk(7), rr, 1 + rr, 2 + rr, kk(-3), kk(11), kk(13), 3 + rr]
    basis = []
    for c in cands:
        if len(basis) == 4: break
        indep = True
        for mask in range(1, 2**len(basis)):
            p = kk(1)
            for i in range(len(basis)):
                if mask >> i & 1: p *= basis[i]
            if is_square_P2(c * p): indep = False; break
        if indep and not is_square_P2(c):
            basis.append(c)
    assert len(basis) == 4, ("basis too small", basis)
    G = [[hilbert_2_search(basis[i], basis[j]) for j in range(4)] for i in range(4)]
    _P2_BASIS, _P2_GRAM = basis, G
    return basis, G

def p2_coords(x):
    """координаты класса x в базисе p2_basis (вектор из 4 бит)"""
    basis, G = p2_basis()
    for mask in range(16):
        p = kk(1)
        for i in range(4):
            if mask >> i & 1: p *= basis[i]
        if is_square_P2(x * p):
            return [(mask >> i) & 1 for i in range(4)]
    raise RuntimeError("class not found")

_P2_COORD_CACHE = {}
def hilbert_2(a, b):
    basis, G = p2_basis()
    out = 1
    ca = _P2_COORD_CACHE.get(str(a))
    if ca is None:
        ca = p2_coords(a); _P2_COORD_CACHE[str(a)] = ca
    cb = _P2_COORD_CACHE.get(str(b))
    if cb is None:
        cb = p2_coords(b); _P2_COORD_CACHE[str(b)] = cb
    for i in range(4):
        for j in range(4):
            if ca[i] and cb[j] and G[i][j] == -1: out = -out
    return out

# ---------- вещественные места ----------
def exact_sign_emb(x, sr):
    """знак x = a0 + a1*r при вложении r -> sr*sqrt(165), sr in {+1,-1}; точная арифметика."""
    x = kk(x)
    a0, a1 = QQ(list(x)[0]), QQ(list(x)[1])
    b = a1 * sr          # x -> a0 + b*sqrt(165)
    if b == 0: return int(sign(a0))
    if a0 == 0: return int(sign(b))
    if a0 > 0 and b > 0: return 1
    if a0 < 0 and b < 0: return -1
    # разные знаки: сравниваем a0^2 и 165 b^2
    d = a0**2 - 165 * b**2
    if d == 0: return 0
    return int(sign(a0)) if d > 0 else int(sign(b))

def hilbert_real(a, b, sr):
    return -1 if (exact_sign_emb(a, sr) < 0 and exact_sign_emb(b, sr) < 0) else 1

# ---------- единый диспетчер ----------
def hsym(place, a, b):
    """place: ('inf', sr) или идеал P"""
    if isinstance(place, tuple) and place[0] == 'inf':
        return hilbert_real(a, b, place[1])
    P = place
    if P.smallest_integer() == 2:
        return hilbert_2(a, b)
    return hilbert_odd(P, a, b)

def is_local_square(place, x):
    if isinstance(place, tuple) and place[0] == 'inf':
        return exact_sign_emb(x, place[1]) > 0
    P = place
    if P.smallest_integer() == 2: return is_square_P2(x)
    return is_square_odd(P, x)
