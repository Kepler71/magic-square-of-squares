# Fable, 14.09.2026. Геометрия границы «клетка = 0» на поверхности девяти квадратов
# и проверка условий Корвахи–Цаннье / Левина. Запуск: python3 geometry.py  (Sage через sage.all)
from sage.all import *
import itertools, json

out = {}
R = PolynomialRing(QQ, 'a,b,c'); a, b, c = R.gens()
# девять клеток как линейные формы на P^2 = пространство магических квадратов [a:b:c], центр a
cells = {'M': a, 'a+b': a+b, 'a-b': a-b, 'a+c': a+c, 'a-c': a-c,
         'a+b+c': a+b+c, 'a-b-c': a-b-c, 'a+b-c': a+b-c, 'a-b+c': a-b+c}
names = list(cells)
lines = [cells[n] for n in names]

# ---------- 1. инцидентность девяти прямых ----------
def meet(l1, l2):
    M = matrix(QQ, [[l.monomial_coefficient(v) for v in (a, b, c)] for l in (l1, l2)])
    K = M.right_kernel().basis()
    assert len(K) == 1
    v = K[0]; v = v / v[[i for i in range(3) if v[i] != 0][0]]
    return tuple(v)
pts = {}
for i, j in itertools.combinations(range(9), 2):
    pts.setdefault(meet(lines[i], lines[j]), set()).update([i, j])
triple = {p: sorted(names[i] for i in s) for p, s in pts.items() if len(s) >= 3}
out['intersection_points'] = len(pts)
out['triple_points'] = {str(p): s for p, s in triple.items()}
out['max_multiplicity'] = max(len(s) for s in pts.values())
# ОДТ на V: над тройной точкой 2^(9-3)/2 = 32 точки
out['ODP_count_predicted'] = 32 * len(triple)
# 5-подмножества в общем положении
gp5 = [T for T in itertools.combinations(range(9), 5)
       if not any(len(set(T) & s) >= 3 for s in pts.values())]
out['five_subsets_general_position'] = len(gp5)
out['example_gp5'] = [names[i] for i in gp5[0]]
gp4_corners = not any(len(set([1,2,3,4]) & s) >= 3 for s in pts.values())
out['four_corner_lines_general_position'] = gp4_corners

# ---------- 2. компоненты и роды прообразов прямых на V ----------
# На прямой ℓ_i = 0 берём параметр t; остальные 8 клеток — линейные формы от t (или константы).
# Классы по модулю квадратов над Qbar: вектор чётностей нулей; над Q добавляем константы (-1, 2, 3, ...).
Rt = PolynomialRing(QQ, 't'); t = Rt.gen()
def restrict(name):
    l = cells[name]
    # параметризация прямой: две точки
    M = matrix(QQ, [[l.monomial_coefficient(v) for v in (a, b, c)]])
    K = M.right_kernel().basis()
    P0, P1 = K[0], K[1]
    sub = {a: P0[0] + t*P1[0], b: P0[1] + t*P1[1], c: P0[2] + t*P1[2]}
    return {n: cells[n].subs(sub) for n in names if n != name}, (P0, P1)

comp = {}
for name in names:
    rest, par = restrict(name)
    polys = list(rest.values())
    # точки ветвления кандидаты: нули линейных форм и бесконечность
    roots = sorted(set(-f[0]/f[1] for f in polys if f.degree() == 1))
    # вектор чётностей над Qbar: по каждому корню + бесконечность (степень mod 2)
    def vec(f):
        v = [0]*(len(roots)+1)
        if f.degree() == 1:
            v[roots.index(-f[0]/f[1])] = 1
            v[-1] = 1
        return vector(GF(2), v)
    # классы относительно первой клетки (проективное отождествление: отношения ℓ_j/ℓ_1)
    f1 = polys[0]
    vecs = [vec(f) + vec(f1) for f in polys[1:]]
    Mv = matrix(GF(2), vecs)
    k = Mv.rank()
    # точки ветвления = координаты, где какой-то элемент образа имеет 1
    span = Mv.row_space()
    branch = [i for i in range(len(roots)+1) if any(v[i] == 1 for v in span.basis())]
    B = len(branch)
    g = 1 - 2**k + B * 2**(k-2)
    ncomp_geom = 2**7 // 2**k          # 128 листов над прямой (x_i = 0), деление на степень компоненты
    # над Q: добавить классы констант (ведущие коэффициенты и значения констант)
    consts = set()
    for f in polys:
        lc = f.leading_coefficient() if f.degree() == 1 else f[0]
        consts.add(QQ(lc / (f1.leading_coefficient() if f1.degree()==1 else f1[0])))
    sqfree = set()
    for q in consts:
        q = QQ(q); s = sign(q) * prod(p for p, e in factor(abs(q.numerator()*q.denominator())) if e % 2)
        if s != 1: sqfree.add(s)
    comp[name] = dict(rank_Qbar=k, degree_component=2**k, components_Qbar=ncomp_geom,
                      branch_points=B, genus=g, roots=[str(r) for r in roots],
                      constant_classes_over_Q=sorted(str(s) for s in sqfree))
out['boundary_curves'] = comp

# прямой контроль рода через поля функций (центр и ребро; угол — степень 32)
def genus_functionfield(name, maxdeg=16):
    rest, _ = restrict(name)
    polys = list(rest.values())
    f1 = polys[0]
    F = FunctionField(QQ, 't'); tt = F.gen()
    K = F
    gens_added = 0
    for f in polys[1:]:
        h = (f / f1)(tt)
        # добавляем sqrt(h), если ещё не квадрат (проверка через норму/факторизацию — грубо: пробуем)
        Ry = PolynomialRing(K, 'y'); y = Ry.gen()
        pol = y**2 - h
        if pol.is_irreducible():
            K = K.extension(pol, 'y%d' % gens_added); gens_added += 1
            if 2**gens_added > maxdeg: return None
    return int(K.genus()), 2**gens_added
try:
    out['genus_functionfield_center'] = genus_functionfield('M')
    out['genus_functionfield_edge'] = genus_functionfield('a+b+c')
except Exception as e:
    out['genus_functionfield_error'] = str(e)

# ---------- 3. Квадрика Q: x^2+y^2 = w^2+t^2 (модель для T = {1+p, 1+q, 1+p+q}), сечения класса (1,1) ----------
def cz_check(r, Dsq=2, Dij=2, p=1):
    # D = p*sum D_i, D_i^2 = Dsq, D_i.D_j = Dij
    D2 = p**2 * (r*Dsq + r*(r-1)*Dij)
    DDi = p * (Dsq + (r-1)*Dij)
    x = polygen(QQ, 'x')
    sols = [s for s in (Dsq*x**2 - 2*DDi*x + D2).roots(ring=AA, multiplicities=False) if s > 0]
    xi = min(sols)
    lhs = 2*D2*xi; rhs = DDi*xi**2 + 3*D2*p
    return dict(r=r, D2=D2, DDi=DDi, xi=str(xi), lhs=str(lhs), rhs=str(rhs), holds=bool(lhs > rhs))
out['CZ_quadric'] = [cz_check(r) for r in (2, 3, 4, 5)]
# общее положение четырёх сечений x,y,w,t на квадрике
Pq = PolynomialRing(QQ, 'x,y,w,t'); x, y, w, tq = Pq.gens()
Qq = x**2 + y**2 - w**2 - tq**2
gp = True
for trio in itertools.combinations([x, y, w, tq], 3):
    I = Pq.ideal([Qq] + list(trio))
    if I.dimension() > 0: gp = False   # проективно: размерность аффинного конуса >0 означает точку
out['quadric_four_sections_common_point_free'] = gp  # ожидаем False? см. ниже: dimension конуса
# точнее: ищем общие проективные точки
def common_points(trio):
    I = Pq.ideal([Qq] + list(trio))
    return I.dimension()  # 1 = только вершина конуса (нет проективных точек); >=2 = есть
out['quadric_trio_cone_dims'] = {str(trio): common_points(trio) for trio in itertools.combinations([x, y, w, tq], 3)}

# ---------- 4. различность корней в слоях ----------
Rs = PolynomialRing(QQ, 's'); s = Rs.gen()
# угловая клетка: p = p0 фиксировано; формы по q: 1+q,1-q,1+p0+q,1+p0-q  и  1+q,1-q,1-p0+q,1-p0-q
def bad_values(roots):
    bad = set()
    for r1, r2 in itertools.combinations(roots, 2):
        d = r1 - r2
        if d == 0: return 'always'
        bad.update(d.roots(multiplicities=False))
    return sorted(bad)
out['corner_fibre_bad_p0_set1'] = [str(v) for v in bad_values([Rs(-1), Rs(1), -1-s, 1+s])]
out['corner_fibre_bad_p0_set2'] = [str(v) for v in bad_values([Rs(-1), Rs(1), -1+s, 1-s])]
# рёберная клетка: p+q = s фиксировано; корни по p: -1, 1, 1+s, s-1, (s-1)/2, (s+1)/2
edge_roots = [Rs(-1), Rs(1), 1+s, s-1, (s-1)/2, (s+1)/2]
crit = set()
for r1, r2 in itertools.combinations(edge_roots, 2):
    crit.update((r1-r2).roots(multiplicities=False))
mind = {}
for v in sorted(crit):
    vals = set(r(v) for r in edge_roots)
    mind[str(v)] = len(vals)
out['edge_fibre_distinct_roots_at_critical_s'] = mind
out['edge_fibre_min_distinct_roots'] = min(mind.values())

# ---------- 5. численная иллюстрация леммы об одной клетке: x S-единица, 2 - x^2 = z^2 ----------
def lemma_demo(S, E):
    sols = set()
    ps = list(S)
    for exps in itertools.product(range(-E, E+1), repeat=len(ps)):
        xv = prod(QQ(p)**e for p, e in zip(ps, exps))
        for sg in (1, -1):
            xx = sg*xv
            z2 = 2 - xx**2
            if z2 >= 0 and z2.is_square():
                sols.add((xx, xx**2))
    return sorted(sols)
demo = lemma_demo([2, 3, 5, 7], 4)
out['lemma_demo_S_2357_expbound4'] = [(str(u), str(v)) for u, v in demo]
demo2 = lemma_demo([2, 3, 5, 7, 11, 13], 3)
out['lemma_demo_S_2to13_expbound3_count'] = len(demo2)
out['lemma_demo_S_2to13_cells'] = sorted(set(str(v) for u, v in demo2))

json.dump(out, default=lambda o: int(o) if hasattr(o, "__int__") else str(o), fp=open('/home/kep/magicKube/fable_integral/geometry.json', 'w'), indent=1, ensure_ascii=False)
for k_, v_ in out.items():
    print(k_, '=', v_)
