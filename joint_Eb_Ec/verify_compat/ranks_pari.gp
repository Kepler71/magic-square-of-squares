\\ Независимая проверка рангов (PARI/GP, без eclib/mwrank).
\\ E_t : y^2 = x^3 - t^2 x.  Контроль Codex: b=34, c=3400, b+c=3434, c-b=3366.
default(parisize, 10^9);
default(realprecision, 60);
sqf(t) = core(t);
{
L = [34, 3400, 3434, 3366, 374];
for (k = 1, #L,
  n = L[k];
  E = ellinit([0,0,0,-n^2,0]);
  print("==== n = ", n, "   core(n) = ", core(n), "   N(E) = ", ellglobalred(E)[1], "  factor(N) = ", factor(ellglobalred(E)[1]));
  print("  torsion: ", elltors(E)[1..2]);
  print("  root number w = ", ellrootno(E));
  R = ellrank(E, 2);
  print("  ellrank(effort=2): r1 = ", R[1], "  r2 = ", R[2], "  s = ", R[3]);
  P = R[4];
  print("  ellrank points: ", P);
  if (#P > 0, print("  det heightmatrix(ellrank pts) = ", matdet(ellheightmatrix(E, P))));
  AR = ellanalyticrank(E);
  print("  ellanalyticrank: order = ", AR[1], "  L^(r)(1) = ", AR[2]);
);
}
