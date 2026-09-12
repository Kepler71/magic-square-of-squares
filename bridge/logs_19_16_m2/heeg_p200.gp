\p 200
default(parisize, 12000000000);
E = ellinit([0,0,0,-1613212693425900,7222574158576522830000]);
print("start p200");
gettime();
P = ellheegner(E);
print("TIME_MS=", gettime());
print("HEEGNER_POINT=", P);
print("oncurve=", ellisoncurve(E,P));
print("height=", ellheight(E,P));
quit
