# Сетка по семье E_{m,n}: y^2 = x(x+m^2)(x+n^2), gcd(m,n)=1, 1<=m<n<=N.
# Для каждой пары: ellrank [r1,r2,s], знак функционального уравнения, кондуктор; C=dim Sel_2=r2+2+s.
import json,sys,time
N=int(sys.argv[1]); OUT=sys.argv[2]
res={}; t0=time.time(); cnt=0
for n in range(2,N+1):
    for m in range(1,n):
        if gcd(m,n)!=1: continue
        E=EllipticCurve([0,m*m+n*n,0,m*m*n*n,0])
        M=E.minimal_model()
        try:
            alarm(30); rk=M.pari_curve().ellrank(); cancel_alarm(); rk=[int(rk[0]),int(rk[1]),int(rk[2])]
        except Exception as e:
            cancel_alarm(); rk=None
        res[f'{m},{n}']={'rk':rk,'w':int(M.root_number()),'N':int(M.conductor())}
        cnt+=1
    if n%20==0:
        print(n,cnt,round(time.time()-t0),flush=True); json.dump(res,open(OUT,'w'))
json.dump(res,open(OUT,'w')); print('done',cnt,round(time.time()-t0))
