# Claude, 26.09.2026. Контроль нормировки: строгая вилка hhat_x = lim 4^-n h(x(2^n P)) из сертификата Безу
# против Sage P.height() на 37a1 и 389a1 (таблицы Кремоны: регулятор 37a1 = 0.0511114082...).
from sage.all import *
import sys
sys.path.insert(0, '/home/kep/magicKube/rigorous_census/chunk0')
from rig_dem import dup_cert, hgt, RIF
for lab in ('37a1', '389a1', '5077a1'):
    E = EllipticCurve(lab)
    c, Fc, Gc = dup_cert(E)
    Bup = RIF(c['L']).log() / 3; Blow = RIF(c['K']).log() / 3
    Px = PolynomialRing(QQ, 'x'); F = Px(Fc); G = Px(Gc)
    for P in E.gens():
        for k in (4, 8):
            xk = P[0]
            for _ in range(k):
                xk = F(xk) / G(xk)
            hk = hgt(xk)
            lo = ((hk - Blow) / 4**k).lower(); hi = ((hk + Bup) / 4**k).upper()
            print(lab, P, 'k=%d' % k, 'вилка [%.10f, %.10f]' % (lo, hi), 'Sage %.10f' % P.height(),
                  'внутри' if lo <= P.height() <= hi else 'ВНЕ', 'B_rig=%.3f B_S=%.3f' % (max(Bup.upper(), Blow.upper()), E.silverman_height_bound()))
