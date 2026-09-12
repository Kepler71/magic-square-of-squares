default(parisize, 6000000000);
E = ellinit([0,0,0,-1613212693425900,7222574158576522830000]);
for(eff=1,25, gettime(); r = ellrank(E, eff); t=gettime(); \
  print("effort=",eff," bounds=[",r[1],",",r[2],"] sha2dim=",r[3]," npts=",#r[4]," time_ms=",t); \
  if(#r[4]>0, print("POINTS FOUND: ", r[4]); break));
quit
