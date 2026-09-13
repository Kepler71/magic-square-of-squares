H := HyperellipticCurve(Polynomial([0, 2164947841, 25215436800, 22212184318, 25215436800, 2164947841]));
pts, ok := RationalPointsGenus2(H);
print "POINTS:", pts;
print "COMPLETE:", ok;
