default(parisize, 6000000000);
E = ellinit([0,0,0,-1613212693425900,7222574158576522830000]);
for(eff=1,20, gettime(); r = ellrank(E, eff); print("effort=",eff," result=",r," time_ms=",gettime()); if(#r[3]>0, print("POINTS FOUND: ", r[3]); break));
quit
