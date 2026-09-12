"""Перекрёстная сверка: МОЯ тройка (g1,g2,g3) для (19,5) подаётся в ДРУГУЮ, уже имевшуюся
в проекте реализацию descent/ctp_quartic.sage (класс FisherCTP).  Значение должно совпасть.
Это сверка двух независимых кодов на одних и тех же объектах.
"""
from sage.all import *
import json
load('/home/kep/magicKube/descent/ctp_quartic.sage')

cert = json.load(open('/home/kep/magicKube/bridge/ctp_cert_19_5.json'))
I, J = QQ(cert['I']), QQ(cert['J'])
gs = [[QQ(c) for c in g] for g in cert['quartics']]
A, B = QQ(cert['model_A']), QQ(cert['model_B'])
M = EllipticCurve([0, 0, 0, A, B])
I2, J2 = IJ_of_curve(M)
print('I,J сходятся:', (I, J) == (I2, J2))
Fc = FisherCTP(QQ, I, J, verbose=True)
for g in gs:
    Fc.check_quartic(g)
print('check_quartic прошёл для всех трёх')
print('z(g1) =', Fc.z_inv(gs[0]))
for reps in (1, 2):
    v, npl = Fc.pair(gs[0], gs[1], gs[2], reps=reps)
    print('ЧУЖАЯ реализация: <g1,g2> =', v, ' (0=тривиально, 1=нетривиально), мест:', npl)
v, npl = Fc.pair(gs[1], gs[0], gs[2], reps=2)
print('ЧУЖАЯ реализация, обратный порядок: <g2,g1> =', v, ' мест:', npl)
v, npl = Fc.pair(gs[0], gs[0], Fc.quartic_from_delta(Fc.E.from_comps([QQ(1)]*3), tag='t'), reps=1)
print('ЧУЖАЯ реализация, диагональ <g1,g1> =', v)
