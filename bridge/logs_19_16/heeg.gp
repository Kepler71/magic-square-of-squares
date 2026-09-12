\p 300
default(parisize, 4000000000);
E = ellinit([0,0,0,-1613212693425900,7222574158576522830000]);
print("cond=", ellglobalred(E)[1]);
gettime();
P = ellheegner(E);
print("TIME_MS=", gettime());
print("HEEGNER_POINT=", P);
print("on curve: ", ellisoncurve(E,P));
print("height=", ellheight(E,P));
quit
