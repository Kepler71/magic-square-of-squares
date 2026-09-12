# Эллиптический Шаботи: строгий подсчёт по РАЗЛИЧНЫМ p-адическим "линиям".
#
# E(k_P) = Ehat x E(F_P) канонически (p нечётно, p | #E(F_P) нет): T = tau(Tbar) для кручения
# (Ehat без кручения), поэтому локальная точка класса (T,m0) есть
#       P(m) = tau(m0*Gbar + Tbar) + m*Ghat,   m in Z_p,
# и множество (линия) зависит ТОЛЬКО от пары редукций a = (m0*Gbar_i + Tbar_i)_{i=1,2},
# т.е. от образа точки в E(F_P1) x E(F_P2).  Классы (T,m0) с одинаковым a дают ОДНУ линию
# (это и есть "призраки при m = -1/2" в наивном перечислении: они дублируют реальные точки).
# Линий ровно |образ E(k) -> E(F_P1) x E(F_P2)|.
#
# Симметрия: -P(m) = tau(-a) + (-m)*Ghat, x(-P) = x(P).  Если -a = a (обе координаты 2-кручение),
# то theta_L(m) = theta_L(-m) чётна => нуль при m = 0 имеет ЧЁТНЫЙ порядок (>= 2).
load('/home/kep/magicKube/family/independent_chabauty/chab.sage')
import sys
b, h, n, D, u0, p = sys.argv[1:7]
prec = int(sys.argv[7]) if len(sys.argv) > 7 else 60
M = int(sys.argv[8]) if len(sys.argv) > 8 else 24
S = Section(ZZ(b), ZZ(h), ZZ(n), ZZ(D)); E = S.E; k = S.k
G = S.point_from_u(QQ(u0))
ch = Chab(S, G, ZZ(p), prec, M).setup()
pp = ZZ(p)
print(f"=== сечение ({b},{h},{n}), D={D}, k=Q(sqrt({S.dk})), p={p}, prec={prec}, M={M}")
print(f"E = {E.ainvs()};  G: u={S.u_of(G)}, h(G)={G.height()}")
print(f"|E(F_P_i)| = {[ch.Ebar(s).cardinality() for s in (1,-1)]}, ord(Gbar) = {ch.orders}, N = {ch.N}")
assert gcd(ch.N, pp) == 1, "p | N -- нужна более тонкая разбивка"
for sgn in (1, -1):
    assert ch.Ebar(sgn).cardinality() % pp != 0, "p | #E(F_P): каноническое расщепление не гарантировано"
print(f"v_p(t_R) = {[t.valuation() for t in ch.tR]}, v = {ch.v}")

tors = ch.tors
# --- известные точки с рациональным u: (ключ редукций) -> [(m, T, кратность)]
def key(Q):
    return (str(ch.redpt(Q, 1)), str(ch.redpt(Q, -1)))
def neg_key(Q):
    return (str(-ch.redpt(Q, 1)), str(-ch.redpt(Q, -1)))

known = {}
LIM = 40
for T in tors:
    for m in range(-LIM, LIM + 1):
        P = m*G + T
        uu = S.u_of(P)
        if uu is None or uu in QQ:
            kk = key(P)
            mult = 2 if (m == 0 and kk == neg_key(P)) else 1
            known.setdefault(kk, []).append((m, 'O' if T.is_zero() else str(T[0]), uu, mult))
print(f"известные точки с u in P^1(Q) (|m| <= {LIM}):")
for kk, v in known.items():
    print(f"    {kk}: {v}")

# --- группируем классы (T,m0) по ключу
groups = {}
for T in tors:
    for m0 in range(ch.N):
        groups.setdefault(key(m0*G + T), []).append((T, m0))
print(f"различных линий: {len(groups)}  (из {4*ch.N} пар (T,m0))")

total, total_known, allok = 0, 0, True
for kk, members in sorted(groups.items()):
    bounds = []
    for (T, m0) in members:
        Q = m0*G + T
        co, gamma, br, phico = ch.theta(Q)
        n0, mu, tailok, vals = strassmann(co, ch.v, ch.p, M)
        intok = all(ch.emb(c, sg).valuation() >= 0 for c in phico for sg in (1, -1))
        bndok = all(vals[j] >= j*ch.v - vpfact(j, ch.p) for j in range(M))
        bounds.append((n0, mu, tailok, intok, bndok))
        if not (tailok and intok and bndok): allok = False
    ns = set(x[0] for x in bounds)
    consistent = (len(ns) == 1)
    if not consistent: allok = False
    n0 = max(ns)
    kn = known.get(kk, [])
    mult = sum(x[3] for x in kn)
    total += n0; total_known += mult
    flag = "OK" if mult == n0 else ("!!! ЛИШНИХ " + str(n0 - mult))
    print(f"  линия {kk}: Strassmann {n0} (дубли {sorted(( 'O' if T.is_zero() else str(T[0]), m0) for T,m0 in members)}, "
          f"согласованы={consistent}), известных {mult} {flag}")
print(f"ИТОГО: нулей <= {total}; известные точки дают >= {total_known}")
print(f"технические проверки (хвост, целость коэффициентов, оценка v_p(c_n) >= n*v - v_p(n!)): {allok}")
if total == total_known and allok:
    print("ВЫВОД: полный учёт -- все нули объяснены известными точками.")
else:
    print(f"ВЫВОД: НЕ полный учёт, необъяснённых нулей {total - total_known}.")
