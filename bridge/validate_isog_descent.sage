# -*- coding: utf-8 -*-
load('/home/kep/magicKube/bridge/isog_descent_lib.sage')
from sage.schemes.elliptic_curves.descent_two_isogeny import two_descent_by_two_isogeny

tests = ['14a1','65a1','20a1','24a1','32a1','15a1','389b1','5077b1','17a1','21a1',
         '1122c1','2002f1','858k1','990h1','858f1','546f1','1050o1','3024k1',
         '256a1','256b1','256d1','1024a1','2601a1','4225a1','1155k1','1470q1']
ok_all = True; ntest = 0
for lab in tests:
    try:
        E = EllipticCurve(lab)
        if E.torsion_order() % 2 != 0: continue
        n1,n2,n1p,n2p = two_descent_by_two_isogeny(E)
    except Exception as ex:
        print(lab,"skip:",ex); continue
    T2 = [P for P in E.torsion_points() if P.order()==2]
    mine = []
    for T in T2:
        A,B = model_A_B(E,T)
        try:
            s1 = len(selmer_set(A,B,verbose=False))
            s2 = len(selmer_set(-2*A, A**2-4*B, verbose=False))   # dual side E'=E/<T>
            mine.append((s1,s2))
        except RuntimeError as ex:
            mine.append(("UND",str(ex)[:40]))
    match = any((p==(n2,n2p) or p==(n2p,n2)) for p in mine)
    ntest += 1
    if not match: ok_all = False
    print("%-8s sage (#Sel, #Sel') = (%s, %s) | mine: %s   %s"
          % (lab, n2, n2p, mine, "OK" if match else "*** MISMATCH ***"))
print("tested", ntest, "curves:", "ALL OK" if ok_all else "SOME MISMATCH")
