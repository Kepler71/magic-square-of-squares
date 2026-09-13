SetClassGroupBounds("GRH");
t := Cputime();
H := ReducedMinimalWeierstrassModel(HyperellipticCurve(Polynomial([0, 2164947841, 25215436800, 22212184318, 25215436800, 2164947841])));
print "MODEL:", H;
J := Jacobian(H);
lo, hi := RankBounds(J);
print "RANK_LO:", lo, "RANK_HI:", hi, "t:", Cputime(t);
