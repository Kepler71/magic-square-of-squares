F.<sg>=FunctionField(QQ); R.<T>=PolynomialRing(F); FR=R.fraction_field()
p=(T^2+1)^2; q=T*(T^2-1)
Bf2=FR((1+2*T-T^2)^2)/FR((1+T^2)^2); Hf2=FR((1-2*T-T^2)^2)/FR((1+T^2)^2)
one=FR(1)
def report(tag, cells):
    print(" ",tag)
    for nm,cell in cells:
        g=R((cell*FR(p)).numerator()/ (cell*FR(p)).denominator()) if (cell*FR(p)).denominator().degree()==0 else None
        v=(cell*FR(p))
        assert v.denominator().degree()==0, (nm, v.denominator())
        poly=R(v.numerator()/v.denominator())
        print(f"    {nm}*(1+t^2)^2 = {poly}     [coeff t^3 = {poly[3]}]")
# (a) fix EDGE pair (step X+Y=sg); parametrise the OTHER EDGE pair (d,f)=Bf2,Hf2 ; free = 4 corners
# d=1-X+Y, f=1+X-Y  =>  X-Y = 1-d = 1-Bf2 ; X+Y=sg
XmY=one-Bf2; X=(sg+XmY)/2; Y=(sg-XmY)/2
report("(a) fix EDGE + parametrise EDGE  -> free cells = 4 CORNERS:",
       [('a',one+X),('i',one-X),('c',one+Y),('g',one-Y)])
# (b) fix CORNER pair (a,i) step X=sg ; parametrise OTHER CORNER (c,g)=Bf2,Hf2 => Y=Bf2-1 ; free = 4 edges
X=FR(sg); Y=Bf2-one
report("(b) fix CORNER + parametrise CORNER -> free cells = 4 EDGES:",
       [('h',one+X+Y),('b',one-X-Y),('d',one-X+Y),('f',one+X-Y)])
# (c) MIXED: fix EDGE pair (X+Y=sg) ; parametrise CORNER pair (a,i)=Bf2,Hf2 => X=Bf2-1 ; free = c,g,d,f
X=Bf2-one; Y=FR(sg)-X
report("(c) MIXED fix EDGE + parametrise CORNER -> free cells = 1 corner pair + 1 edge pair:",
       [('c',one+Y),('g',one-Y),('d',one-X+Y),('f',one+X-Y)])
