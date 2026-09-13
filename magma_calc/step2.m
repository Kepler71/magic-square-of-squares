SetClassGroupBounds("GRH");
t := Cputime();
H := ReducedMinimalWeierstrassModel(HyperellipticCurve(Polynomial([0, 4092529, 51379680, 84896062, 81791520, 44182609])));
J := Jacobian(H);
lo, hi := RankBounds(J);
print "RANK_LO:", lo, "RANK_HI:", hi, "t:", Cputime(t);
pts, ok := RationalPointsGenus2(H);
print "POINTS:", pts;
print "COMPLETE:", ok, "t:", Cputime(t);
