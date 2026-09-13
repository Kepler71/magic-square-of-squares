H := HyperellipticCurve(Polynomial([4092529, 0, 47287151, 0, 37608911, 0, 44182609]));
pts, ok := RationalPointsGenus2(H);
print "POINTS:", pts;
print "COMPLETE:", ok;
