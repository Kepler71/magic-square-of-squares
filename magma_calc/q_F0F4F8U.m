H := HyperellipticCurve(Polynomial([0, 44182609, 81791520, 84896062, 51379680, 4092529]));
pts, ok := RationalPointsGenus2(H);
print "POINTS:", pts;
print "COMPLETE:", ok;
