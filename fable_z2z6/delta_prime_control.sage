# Положительный контроль отображения delta' на E': гомоморфизм, обнуляется на phi(E(Q)), delta(T)=4ab(a+b)
load('/home/kep/magicKube/fable_z2z6/three_isog.sage')
K.<w> = QuadraticField(-3)
def dprime(Ep, tang, P):
    x0,s,m1,n1 = tang
    return P[1] - (m1*P[0]+n1)*w
def is_cube_K(z):
    if z == 0: return None
    R.<t> = K[]
    return len((t^3 - z).roots()) > 0
def delta(al, be, P):
    return P[1] - al*P[0] - be
ok = 0; tot = 0
for (a,b) in [(3,5),(7,2),(4,5),(5,7),(11,13),(3,8),(5,9)]:
    E, Ep, phi, al, be = curve_data(a,b)
    tang = tangent_Tprime(Ep)
    Em = E.minimal_model(); iso = Em.isomorphism_to(E)
    gens = [iso(P) for P in Em.gens()]
    T = E(0, be)
    # delta(T) = 4ab(a+b) mod cubes: delta(2T) = -2beta; check delta(-T) * delta(T) cube and 4ab(a+b)*(-2 beta) cube
    r1 = QQ(4*a*b*(a+b)) * QQ(-2*be)
    tot += 1; ok += QQ(r1).is_nth_power(3) if hasattr(QQ(r1),'is_nth_power') else (r1.numerator().is_nth_power(3) and r1.denominator().is_nth_power(3))
    for G in gens:
        P1 = phi(G); tot += 1; ok += bool(is_cube_K(dprime(Ep, tang, P1)))          # delta'(phi(G)) куб
        for Q in [P1]:
            Q2 = 2*Q; 
        # гомоморфизм delta' на E'(Q): возьмём образующие E'(Q)
    Epm = Ep.minimal_model(); isop = Epm.isomorphism_to(Ep)
    gp = [isop(P) for P in Epm.gens()]
    for G in gp:
        d1 = dprime(Ep, tang, G); d2 = dprime(Ep, tang, 2*G)
        tot += 1; ok += bool(is_cube_K(d1^2 / d2))                                   # delta'(2G) = delta'(G)^2 mod кубов
        # delta(phihat(G)) куб на E:
        phihat = phi.dual()
        H = phihat(G); dd = delta(al, be, H)
        if dd != 0:
            tot += 1; ok += bool(QQ(dd).numerator().is_nth_power(3) and QQ(dd).denominator().is_nth_power(3)) if False else bool(is_cube_K(K(dd)))
    print((a,b), "rank E", len(gens), "rank E'", len(gp))
print("controls passed", ok, "of", tot)
