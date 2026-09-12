# -*- coding: utf-8 -*-
# Оценка Reg*#Sha для всех кривых класса изогении (19,16): где генератор ниже?
import time
E = load("/home/kep/magicKube/bridge/logs_19_16_claude3/E.sobj")
cls = E.isogeny_class()
print("isogeny matrix:\n", cls.matrix())
t0 = time.time()
for i, EE in enumerate(cls.curves):
    om   = EE.period_lattice().omega()          # real period (Sage: omega() = real period, incl. 2 if disc>0)
    tam  = prod([EE.tamagawa_number(p) for p in EE.conductor().prime_factors()])
    tors = EE.torsion_order()
    print("-"*60)
    print("curve %d: %s" % (i, EE.ainvs()))
    print("  omega=%s tam=%s tors=%s tam list=%s" % (om, tam, tors,
          [(p, EE.tamagawa_number(p)) for p in EE.conductor().prime_factors()]))
print("time %.1f"%(time.time()-t0))
# L'(1) via PARI on the first curve; L-function is isogeny-invariant
t0=time.time()
pari.allocatemem(4*10^9, silent=True)
Lp = pari(E).ellL1(1)     # L'(1) for rank 1
print("L'(1) =", Lp, "  (%.1f s)"%(time.time()-t0))
for i, EE in enumerate(cls.curves):
    om   = EE.period_lattice().omega()
    tam  = prod([EE.tamagawa_number(p) for p in EE.conductor().prime_factors()])
    tors = EE.torsion_order()
    RegSha = RR(Lp) * tors^2 / (RR(om) * tam)
    print("curve %d : Reg*#Sha = %s   -> if Sha=4: Reg=%s ; if Sha=16: %s" %
          (i, RegSha, RegSha/4, RegSha/16))
