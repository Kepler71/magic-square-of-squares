from sage.all import *
from pathlib import Path
import json,hashlib

out=Path(__file__).resolve().parent
source=out/'vendor'/'ctp_quartic.sage'
load(str(source))
random.seed(20260913)
E=EllipticCurve([0,1,0,-68110640,204180681300])
P=PolynomialRing(QQ,'x');x=P.gen()
covers=E.pari_curve().ell2cover()
gs=[[QQ(P(c[0])[4-j]) for j in range(5)] for c in covers]
I=quartic_I(gs[0]);J=quartic_J(gs[0])
assert all(quartic_I(g)==I and quartic_J(g)==J for g in gs)
F=FisherCTP(QQ,I,J,verbose=True)
g1=F._make_z_unit(gs[2]);g2=F._make_z_unit(gs[3])
print('I,J',I,J,flush=True)
print('g1,g2',g1,g2,flush=True)
g3=F.quartic_from_delta(F.z_inv(g1)*F.z_inv(g2),tag='sum_sha_candidates')
print('g3',g3,flush=True)
assert F.E.is_square(F.z_inv(g1)*F.z_inv(g2)*F.z_inv(g3))
record={'curve_ainvs':list(map(str,E.ainvs())),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'I':str(I),'J':str(J),'quartics':[[str(c) for c in g] for g in [g1,g2,g3]],
        'covers_origin':'PARI ell2cover; pairing recomputed with the separate project Fisher implementation',
        'runs':[]}
from sage.libs.eclib.interface import mwrank_EllipticCurve
ec=mwrank_EllipticCurve(list(map(int,E.ainvs())),verbose=False);ec.two_descent(verbose=False)
record['eclib_selmer_dimension']=int(ec.selmer_rank())
record['rational_two_torsion_dimension']=2
assert record['eclib_selmer_dimension']==4
EI=EllipticCurve([0,0,0,-27*I,-27*J])
iso=EI.isomorphism_to(E)
record['quartic_jacobian_model']=list(map(str,EI.ainvs()))
record['jacobian_to_original_isomorphism']=list(map(str,iso.tuple()))
orig_place=F.places_for;orig_point=F.local_point
orig_gamma=F.gamma1
active=None
def gamma(left,right,third):
    gam,m=orig_gamma(left,right,third)
    active.update(raw_gamma=list(map(str,gam)),m_coefficients=list(map(str,F.E.coeffs(m))))
    return gam,m
def places(g,gam,a,extra_norm=16):
    pls=orig_place(g,gam,a,extra_norm=extra_norm)
    active.update(gamma=list(map(str,gam)),a=str(a),places=[pl.name for pl in pls],local=[])
    return pls
def point(g,gam,pl,tries=60000):
    xx,zz=orig_point(g,gam,pl,tries=tries)
    gv=gam[0]*xx*xx+gam[1]*xx*zz+gam[2]*zz*zz
    val=quartic_eval(g,xx,zz)
    if isinstance(pl,RealPlace):assert val>0
    else:assert pl.is_sq(val)
    assert gv!=0
    hs=hilb(QQ,QQ(active['a']),gv,pl)
    active['local'].append({'place':pl.name,'x':str(xx),'z':str(zz),'quartic_value':str(val),
                            'gamma_value':str(gv),'hilbert':int(hs)})
    return xx,zz
F.places_for=places;F.local_point=point;F.gamma1=gamma
for left,right,label in [(g1,g2,'forward'),(g2,g1,'reverse')]:
    active={'direction':label};record['runs'].append(active)
    F._lpcache={}
    bit,n=F.pair(left,right,g3,reps=1,verbose=True)
    active.update(pairing_bit=int(bit),number_of_places=n)
    (out/'pairing_40_43.json').write_text(json.dumps(record,indent=2)+'\n')
    print(label,bit,n,flush=True)
assert all(r['pairing_bit']==1 for r in record['runs'])
record['nontrivial_pairing_both_directions']=True
record['everywhere_local_solubility']=[]
for index,g in enumerate([g1,g2,g3]):
    loc=[];F._lpcache={}
    for pl in orig_place(g,[QQ(1),QQ(0),QQ(0)],g[0],extra_norm=16):
        xx,zz=orig_point(g,[QQ(1),QQ(0),QQ(0)],pl)
        val=quartic_eval(g,xx,zz)
        assert val>0 if isinstance(pl,RealPlace) else pl.is_sq(val)
        loc.append({'place':pl.name,'x':str(xx),'z':str(zz),'quartic_value':str(val)})
    record['everywhere_local_solubility'].append({'quartic_index':index,'local_points':loc})
record['rank_upper_bound_from_pairing']=record['eclib_selmer_dimension']-2-2
assert record['rank_upper_bound_from_pairing']==0
(out/'pairing_40_43.json').write_text(json.dumps(record,indent=2)+'\n')
