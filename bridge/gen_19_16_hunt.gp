default(parisize, 6000000000);
/* четыре кривые класса изогении (19,16) */
V = [ [0,0,0,-1613212693425900, 7222574158576522830000],
      [0,0,0,-1269386570673900, 17388781634016192078000],
      [0,0,0, 6144364201166100, 56474886831958509390000],
      [0,0,0,-14872007552049900,-692667016942944295602000] ];
{
for(i=1,4,
  E = ellinit(V[i]);
  rk = ellrankinit(E);
  for(eff=0,14,
    t = getabstime();
    r = ellrank(rk, eff);
    printf("curve %d effort %d -> [%d,%d,%d] npts=%d  %.1fs\n", i, eff, r[1],r[2],r[3], #r[4], (getabstime()-t)/1000.);
    if(#r[4] > 0,
      print("### POINT on curve ", i, " : ", r[4]);
      write("/home/kep/magicKube/bridge/logs_19_16_claude3/found_pts.txt", [i, r[4]]);
      break(1);
    );
  );
);
}
quit;
