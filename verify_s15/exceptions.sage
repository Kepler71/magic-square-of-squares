# П.2: исключения в выводе тождества (1) (Codex, NO_LIFT_MOD7.md §2).
import functools
print = functools.partial(print, flush=True)
R.<t, X, v0> = QQ[]
a = 109*(t^2+1)^2; b = 338*t*(t^2-1); c = 229*(t^2+1)^2
e1 = v0^2 - 4*a*v0 + 4*b^2
e2 = X^4*v0^2 - 4*c*X^2*v0 + 4*b^2
res = e1.resultant(e2, v0)
codex = 16*b^2*(b^2*(1 - X^4)^2 - 4*X^2*(c - a*X^2)*(a - c*X^2))
print("(a) Res_{v0}(e1,e2) == 16 b^2 {b^2(1-X^4)^2 - 4X^2(c-aX^2)(a-cX^2)}  [Codex]:", res == codex)

S.<T> = QQ[]
A_ = 109*(T^2+1)^2; B_ = 338*T*(T^2-1); C_ = 229*(T^2+1)^2
print("(b) rational roots of b: ", sorted(B_.roots(QQ, multiplicities=False)))
for t0 in B_.roots(QQ, multiplicities=False):
    vals = [A_(t0) + B_(t0), A_(t0) - B_(t0), C_(t0) + B_(t0), C_(t0) - B_(t0)]
    print(f"    t = {t0}: Y_i^2 = {vals}; squares: {[QQ(v).is_square() for v in vals]}")
print("(c) t = oo: leading coefficients of a±b, c±b:", [(A_ + B_).leading_coefficient(), (A_ - B_).leading_coefficient(),
      (C_ + B_).leading_coefficient(), (C_ - B_).leading_coefficient()], " 109, 229 squares:", QQ(109).is_square(), QQ(229).is_square())
print("(d) Y1 = Y2 => Y1^2 - Y2^2 = 2b = 0 (тождество):", (A_ + B_) - (A_ - B_) == 2*B_, "; аналогично Y3,Y4:", (C_ + B_) - (C_ - B_) == 2*B_)
print("(e) w = 0 <=> 109p^2 = 229q^2: 229/109 square in Q:", (QQ(229)/109).is_square())
g4 = lambda x: (229 - 109*x^2)*(109 - 229*x^2)
print("(f) X = ±1: RHS of (1) =", g4(1), g4(-1), " (<0, а LHS — квадрат)")
print("(g) t^2+1 = 0 over Q: no rational roots:", (T^2 + 1).roots(QQ) == [])
print("(h) rational t with some Y_i = 0 (a±b = 0 or c±b = 0):", {nm: P.roots(QQ, multiplicities=False) for nm, P in
      [("a+b", A_ + B_), ("a-b", A_ - B_), ("c+b", C_ + B_), ("c-b", C_ - B_)]})
# (i) X = 0 или oo: X = (Y3-Y4)/(Y1-Y2); X = 0 <=> Y3 = Y4 => b = 0; X = oo <=> Y1 = Y2 => b = 0 (см. (d))
print("(i) X = 0 / oo reduce to (d): OK by construction")
