"""Run the published Bianchi--Padurariu code on a reference or G1 fiber.
Downloaded vendor code is left unchanged. Full rational-point status requires
auditing precision and excluding all extra p-adic points.
"""
from sage.all import *
from sage.env import SAGE_VERSION
from pathlib import Path
import sys, json, time, hashlib, traceback

base=Path(__file__).resolve().parent
mode=sys.argv[1] if len(sys.argv)>1 else 'reference'
p=ZZ(sys.argv[2]) if len(sys.argv)>2 else ZZ(5)
prec=ZZ(sys.argv[3]) if len(sys.argv)>3 else ZZ(20)
tag=sys.argv[4] if len(sys.argv)>4 else ''
profile=sys.argv[5] if len(sys.argv)>5 else 'full_H'
assert profile in ['full_H','full_lifting_13']
if profile=='full_lifting_13':
    assert mode=='41_1' and p==13
    def _qc_disc_filter(P):
        # Proven in local_full_lifting_41_1.md: all other discs fail a square cell.
        return P[2]==0 or P[0]==0
prefix=f'qc_{mode}_p{p}_n{prec}'+('_'+tag if tag else '')
R=PolynomialRing(QQ,'x'); x=R.gen()
if mode=='reference':
    f=x**6+22*x**4-19*x**2+4
elif mode=='mixed_15_8':
    f=(49+529*x*x)*(83521*x**4+63358*x*x+83521)
elif mode=='41_1':
    f=1681*x**6+2827443*x**4+2827443*x**2+1681
else:
    raise ValueError(mode)
code=base/'vendor/qc_g2_mixed_full_lift.sage'
assert mode=='mixed_15_8' and p==11
profile='full_five_square_lift_only'
def _qc_disc_filter(P):
    return P[2]==0 or P[0]==0
out={'mode':mode,'polynomial':str(f),'p':int(p),'precision':int(prec),
     'sage_version':SAGE_VERSION,'vendor_url':'https://github.com/bianchifrancesca/QC_bielliptic',
     'vendor_sha256':hashlib.sha256(code.read_bytes()).hexdigest(),
     'profile':profile,'status':'running'}
dest=base/(prefix+'.json')
dest.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2),flush=True)
start=time.monotonic()
diagnostics=[]
def _qc_diagnostic_hook(entry):
    diagnostics.append(entry)
    (base/(prefix+'_functions.json')).write_text(json.dumps(diagnostics,indent=2)+'\n')
def _qc_disc_audit_hook(entry):
    out['disc_coverage']=entry
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print('QC disc coverage: '+json.dumps(entry),flush=True)
try:
    load(str(code))
    E1=EllipticCurve([0,f[4],0,f[2]*f[6],f[0]*f[6]**2])
    E2=EllipticCurve([0,f[2],0,f[0]*f[4],f[0]**2*f[6]])
    out['ellrank_intervals']=[str(pari(E).ellrank()) for E in [E1,E2]]
    assert all(E.rank(algorithm='pari',proof=True,use_database=False)==1 for E in [E1,E2])
    assert f.discriminant().valuation(p)==0
    assert all(E.has_good_reduction(p) and E.is_ordinary(p) for E in [E1,E2])
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print('Calling quadratic_chabauty_bielliptic',flush=True)
    rat,extra=quadratic_chabauty_bielliptic(f,p,prec,omega_info=True)
    assert not isinstance(rat,str) and not isinstance(extra,str)
    H=HyperellipticCurve(f)
    for row in rat:
        for P in row:
            if P[2] != 0: assert P[1]**2==f(P[0])
    out['rational_points_by_omega']=[[str(P) for P in row] for row in rat]
    out['extra_points_by_omega']=[[str(P) for P in row] for row in extra]
    out['rational_affine_coordinates']=[list(map(str,P[:2])) for row in rat for P in row if P[2]!=0]
    out['extra_counts']=[len(row) for row in extra]
    out['status']='completed; full certificate still requires audit'
except Exception as ex:
    out['status']='failed'
    out['error']=repr(ex)
    out['traceback']=traceback.format_exc()
    print(out['traceback'],flush=True)
out['elapsed_seconds']=time.monotonic()-start
dest.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2),flush=True)
