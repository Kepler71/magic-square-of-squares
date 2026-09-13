H := HyperellipticCurve(Polynomial([0, 4092529, 51379680, 84896062, 81791520, 44182609]));
pts, ok := RationalPointsGenus2(H);
print "POINTS:", pts;
print "COMPLETE:", ok;
