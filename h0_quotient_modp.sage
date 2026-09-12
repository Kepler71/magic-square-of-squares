import functools, time
print = functools.partial(print, flush=True)
src = open('h0_quotient.sage').read()
exec(preparse(src.split('t0 = time.time(); g1 = Curve(G)')[0]))
print("H0 at fixed points of iota (±1,±1):", {(a_, b_): H0(a_, b_) for a_ in (1, -1) for b_ in (1, -1)})
p = 101
Rp = PolynomialRing(GF(p), 'k,u'); kp, up = Rp.gens()
t0 = time.time(); gH = Curve(Rp(H0(kp, up))).geometric_genus(); print(f"control: genus(H0 mod {p}) = {gH} (expected 22) ({time.time()-t0:.1f}s)")
Sp = PolynomialRing(GF(p), 's,w'); sp, wp = Sp.gens()
Mp = Sp(G(sp, wp))
print("M mod p irreducible:", Mp.is_irreducible())
t0 = time.time(); gM = Curve(Mp).geometric_genus(); print(f"genus(M mod {p}) = {gM}  => r = 46 - 4g' = {46 - 4*gM} ({time.time()-t0:.1f}s)")
