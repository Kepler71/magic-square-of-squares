import sys
load('/home/kep/magicKube/fable_z2z6/three_isog.sage')
for (a,b) in [(6,17),(10,21),(19,35),(40,41),(27,28)]:
    E, Ep, phi, alpha, beta = curve_data(a,b)
    r = selmer_3isog(a,b)
    Em = E.minimal_model(); Epm = Ep.minimal_model()
    print((a,b), "bound3=", r['bound'], "dims", r['dim_hat'], r['dim_phi'], "rank E=", pari(Em).ellrank(), "rank E'=", pari(Epm).ellrank())
    try:
        print("   Sha_an(E)=", Em.sha().an(), " Sha_an(E')=", Epm.sha().an(), " tors E'=", Epm.torsion_order(), "cond=", Em.conductor().factor())
    except Exception as e:
        print("   sha error", e)
    print("   local dims E:", r['imE'] and {p: len(v) for p,v in r['imE'].items()}, " E':", {p: len(v) for p,v in r['imEp'].items()})
    print("   sel_hat:", r['sel_hat'])
