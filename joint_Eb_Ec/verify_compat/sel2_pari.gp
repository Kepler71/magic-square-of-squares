\\ dim_F2 Sel_2(E_n) по PARI ell2cover (для контроля формулы Монски)
{
for (n = 1, 600, if (issquarefree(n),
  E = ellinit([0,0,0,-n^2,0]);
  print(n, " ", #ell2cover(E))));
foreach([34,3434,374], n, E = ellinit([0,0,0,-n^2,0]); print(n, " ", #ell2cover(E)));
}
