default(parisize, 4000000000);
m=19; n=16; s=(m^2+n^2)/2; b=s*m^2*n^2;
A=4*b; B=4*s*m^4; C=4*s*n^4;
E = ellinit([0, A+B+C, 0, A*B+A*C+B*C, A*B*C]);
Em = ellminimalmodel(E);
print("ainvs ", [Em.a1,Em.a2,Em.a3,Em.a4,Em.a6]);
{
for(eff=1,60,
  t=getabstime();
  r=ellrank(Em,eff);
  print("effort=",eff," -> lo=",r[1]," hi=",r[2]," pts=",r[4],"  time=",(getabstime()-t)/1000.,"s");
);
}
quit;
