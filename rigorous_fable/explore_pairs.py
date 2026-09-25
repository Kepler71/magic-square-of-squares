# разведка: для заданных наклона и класса — вырожденность пар, класс полярного дивизора u+, представление n/h
import sys, time
sys.path.insert(0, '/home/kep/magicKube/rigorous_fable')
from fmodels import *
def explore(p, q, Ts, hub, others):
    facs = [Fac(T) for T in Ts]; Et = target_model(facs[hub].E)
    print('slope %d/%d Et=%s cond=%s' % (p, q, Et.ainvs(), Et.conductor()))
    for f in facs: assert f.E.is_isomorphic(Et)
    ims = [Image(f, Et) for f in facs]
    for j in others:
        A, B, D, I = pair_AB(ims[hub], ims[j])
        chk = exact_pair_check(ims[hub], ims[j], A, B, D, QQ(1)/(10*max(abs(c) for c in cells(p,q))))
        print(' pair', Ts[hub], Ts[j], 'D=', D, 'exact check', chk)
        if len(D) == 2: print('  conic; deg A', A.numerator().degree(), A.denominator().degree()); continue
        aux = Aux(A, B, D)
        t0 = time.time(); Dp = aux.up.divisor_of_poles(); Dm = aux.um.divisor_of_poles()
        print('  E3', aux.E3.ainvs(), 'tors', aux.E3.torsion_order(), 'deg poles', Dp.degree(), Dm.degree(), 'equiv D+~D-:', len((Dp - Dm).basis_function_space()), time.time() - t0)
        for S in [None] + aux.two_torsion():
            sols = solve_linear_forms(aux.up, aux.basis_3O_plus_S(S))
            solm = solve_linear_forms(aux.um, aux.basis_3O_plus_S(S))
            print('   S=', S, 'dim sol u+:', len(sols), 'u-:', len(solm))
if __name__ == '__main__':
    explore(126, 451, [[-577,-451,-325,-126],[-325,-126,126],[-126,126,325]], 1, [2, 0])
    explore(73, 362, [[-435,-362,-289,-73],[-289,-73,73],[-73,73,289]], 1, [2, 0])
    explore(265, 298, [[-563,-298,-265,265],[-563,-33,265],[-265,33,563],[-265,265,298,563]], 0, [1, 3])
    explore(265, 298, [[-563,-298,33],[-298,-33,265,298],[-298,-265,33,298],[-33,298,563]], 2, [0, 1])
    explore(265, 298, [[-563,-298,33],[-298,-33,265,298],[-298,-265,33,298],[-33,298,563]], 1, [2, 3])
