# Claude, 26.09.2026. Перепроверка записанных rig_*.json: только по данным JSON + пересборке кривой.
# (1) тождества Безу — точно; соответствие F/G отображению удвоения модели; (2) B = log max(K,L)/3, c1 по матрице;
# (3) матрица M: z(P) через M и x модели = z(P) через исходную модель для P0, 2P0, 3P0;
# (4) нижняя оценка ĥ(P0) заново, k = 7 удвоений (в прогоне k = 6); M0 по ней не больше записанного;
# (5) mwrank: верхняя граница ранга (rank_bound).
import sys, os, json
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from rigdem import model, verify_cert, dup_map, mob, c1_of, hgt, RIF
from sage.all import *

SLOPES = ['165/439', '422/441', '425/441', '333/442', '111/445', '329/449', '438/449', '126/451',
          '172/451', '295/452', '305/461', '236/475', '165/493', '52/499', '129/499']
report = {}
allok = True
for sl in SLOPES:
    r, s = map(int, sl.split('/'))
    D = json.load(open(os.path.join(here, 'rig_%d_%d.json' % (r, s))))
    rows = []
    for x in D['results']:
        T = x['T']; name = x['model_used']; var = x['variants'][name]; cert = var['cert']
        md = model(T); E1 = md['E1']
        assert [str(a) for a in E1.ainvs()] == x['E1']
        Emin = EllipticCurve(QQ, [QQ(a) for a in x['Emin']])
        assert Emin.is_isomorphic(E1)
        inv = Emin.isomorphism_to(E1)
        EE = Emin if name == 'min' else E1
        # (1)
        for nm in ('min', 'E1'):
            verify_cert(x['variants'][nm]['cert'])
        F, G = dup_map(EE)
        Rq = F.parent()
        Fc = Rq([QQ(c) for c in cert['F']]); Gc = Rq([QQ(c) for c in cert['G']])
        assert Fc * G == Gc * F, 'сертификат не от удвоения этой модели'
        # (2)
        B = RIF(max(ZZ(cert['K']), ZZ(cert['L']))).log() / 3
        M = matrix(ZZ, var['M'])
        c1, _ = c1_of(M)
        C = 2 * B + 2 * c1
        ok2 = abs(float(B.upper()) - x['B_rig']) < 1e-9 and abs(float(c1.upper()) - x['c1']) < 1e-9 and abs(float(C.upper()) - x['C']) < 1e-9
        # (3)
        xy = x['P0'].strip('()').split(',')
        P0 = Emin([QQ(xy[0]), QQ(xy[1])])
        P0m = P0 if name == 'min' else inv(P0)
        ok3 = True
        for k in (1, 2, 3):
            Pm = k * P0m; P1 = k * inv(P0)
            if mob(M, Pm[0]) != mob(md['mu'], P1[0]):
                ok3 = False
        # (4)
        kd = 7
        Q = P0m
        for _ in range(kd):
            Q = 2 * Q
        hlo = (hgt(Q[0]) / 4 ** kd - RIF(ZZ(cert['K'])).log() / (3 * 4 ** kd)).lower()
        M0new = int(floor(((RIF(C.upper()) / RIF(hlo) + 1) / 2).upper()))
        ok4 = M0new <= x['M0'] and hlo > 0
        # (5)
        try:
            mw = Emin.mwrank_curve(); rb = int(mw.rank_bound()); rk = int(mw.rank()); cert_mw = bool(mw.certain())
        except Exception as ex:
            rb = None; rk = None; cert_mw = None
        ok = ok2 and ok3 and ok4
        allok &= ok
        rows.append(dict(T=T, model=name, B=float(B.upper()), c1=float(c1.upper()), hlo_k7=float(hlo), hlo_k6=x['hP0_lower'],
                         M0_run=x['M0'], M0_k7=M0new, ok_numbers=ok2, ok_matrix=ok3, ok_M0=ok4,
                         mwrank_rank_bound=rb, mwrank_rank=rk, mwrank_certain=cert_mw, ellrank=x['ellrank']))
        print(sl, T, 'ok' if ok else 'FAIL', 'M0 %d/%d' % (x['M0'], M0new), 'mwrank bound', rb, 'rank', rk, cert_mw, flush=True)
    report[sl] = rows
json.dump(dict(all_ok=allok, slopes=report), open(os.path.join(here, 'verify_chunk2.json'), 'w'), ensure_ascii=False, indent=1)
print('ВСЁ СОШЛОСЬ' if allok else 'ЕСТЬ РАСХОЖДЕНИЯ')
