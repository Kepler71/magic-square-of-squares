\p 40
default(parisize, 8000000000);
E = ellinit([0,0,0,-1613212693425900,7222574158576522830000]);
print("start p40");
gettime();
P = ellheegner(E);
print("TIME_MS=", gettime());
print("HEEGNER_POINT=", P);
print("oncurve=", ellisoncurve(E,P));
quit
