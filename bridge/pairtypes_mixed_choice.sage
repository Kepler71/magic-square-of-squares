# Control: only "fix a pair + parametrise the OTHER pair OF THE SAME TYPE" gives the shape
#   A1*(t^2+1)^2 +- B*t(t^2-1),  A2*(t^2+1)^2 +- B*t(t^2-1)   with a COMMON shift B.
R.<T>=QQ[]; var('sg')
p=(T^2+1)^2; q=T*(T^2-1)
Bf2=(1+2*T-T^2)^2/(1+T^2)^2; Hf2=(1-2*T-T^2)^2/(1+T^2)^2
def free_scaled(kind):
    # centre 1; sg = step of the FIXED pair
    if kind=='edge-edge':      # fixed EDGE (step sg = X+Y), parametrised EDGE (step X-Y)
        X=(sg+ (Bf2-1))/1  # not used; use direct cell formulas below
        tau=(Bf2-Hf2)/4    # (X-Y)/2 ... handled explicitly instead
    return None
# explicit cell lists, centre 1, scaled by (1+t^2)^2
rho=(Bf2-1)*(1+T^2)^2       # = -4 t(t^2-1)
res={}
# (a) fix EDGE (X+Y=sg), parametrise the other EDGE  -> free = 4 corners
res['edge + edge (free: corners)'] = [ ((1-sg/2)*p + q*2), ((1-sg/2)*p - q*2), ((1+sg/2)*p + q*2), ((1+sg/2)*p - q*2) ]
# (b) fix CORNER (X=sg), parametrise the other CORNER -> free = 4 edges
res['corner + corner (free: edges)'] = [ ((1+sg)*p - 4*q), ((1-sg)*p + 4*q), ((1-sg)*p - 4*q), ((1+sg)*p + 4*q) ]
# (c) fix EDGE (X+Y=sg), parametrise a CORNER (a,i)  -> free = one corner pair + one edge pair
res['edge + corner (mixed)'] = [ ((1+sg)*p + 4*q), ((1-sg)*p - 4*q), ((1+sg)*p + 8*q), ((1-sg)*p - 8*q) ]
for k,v in res.items():
    shifts=set()
    for f in v:
        g=(f).expand()
        shifts.add(abs(g.coefficient(T,3)))
    print(f"  {k:34s} shift magnitudes |coeff of t^3| present: {sorted(shifts)}  -> common shift: {len(shifts)==1}")
# independent check of (c) by direct cell arithmetic
X=rho/(1+T^2)^2
Y=sg-X
c=(1+Y); g_=(1-Y); d=(1-X+Y); f=(1+X-Y)
for nm,cell in [('c',c),('g',g_),('d',d),('f',f)]:
    val=((cell*(1+T^2)^2).expand())
    print(f"   direct: {nm} * (1+t^2)^2 =", val)
