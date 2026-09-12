# M2 и M3 для H0 при заданном p (аргумент)
import sys, time
exec(preparse(open('genus_h0.sage').read().split('# ---------------- контрольные')[0]))
p = int(sys.argv[1])
t0 = time.time(); r2 = M2(H0, p); print(f"p={p} M2 (Hess genus) = {r2} ({time.time()-t0:.1f}s)", flush=True)
t0 = time.time(); r3 = M3(H0, p); print(f"p={p} M3 (Riemann-Hurwitz) = {r3} ({time.time()-t0:.1f}s)", flush=True)
