\\ Проверка образующих Codex (rank_compatibility.json): на кривой, регулятор,
\\ сравнение с точками ellrank, насыщение: 2-насыщение элементарно (отображение Куммера),
\\ нечётные простые < B через ellsaturation (PARI, не eclib).
default(parisize, 10^9);
default(realprecision, 60);
{C = [[34,   [[-16,120],[-2,48]]],
     [3400, [[-1600,120000],[-200,48000]]],
     [3434, [[243049/64,51166005/512],[166617/4,-67779459/8]]],
     [3366, [[-1683/841,116132049/24389]]]];}
\\ squarefree class of a nonzero rational
sq(q) = core(numerator(q)*denominator(q));
\\ Kummer map for y^2 = x(x-n)(x+n): e1=0,e2=n,e3=-n; returns F2-vector over primes {-1,2,p|n}
kum(n, R) = {
  my(x = R[1], a, b);
  if (R == [0], return([1,1]));
  if (x == 0,  a = sq(-n^2), a = sq(x));          \\ x - e1, fix at T=(0,0): (e1-e2)(e1-e3) = -n^2
  if (x == n,  b = sq(2*n^2), b = sq(x - n));     \\ x - e2, fix at T=(n,0): (e2-e1)(e2-e3) = 2n^2
  [a, b];
}
vecF2(n, pair) = {
  my(ps = concat([-1, 2], factor(core(n))[,1]~), v = []);
  ps = Set(ps);
  foreach(pair, d,
    foreach(ps, p, v = concat(v, if (p == -1, sign(d) < 0, valuation(d, p) % 2))));
  Mod(v~, 2);
}
{
for (k = 1, #C,
  n = C[k][1]; V = C[k][2];
  E = ellinit([0,0,0,-n^2,0]);
  print("==== n = ", n);
  for (i = 1, #V, if (!ellisoncurve(E, V[i]), error("point not on curve")));
  print("  all points on curve: yes");
  H = ellheightmatrix(E, V);
  print("  heights (PARI ellheight): ", vector(#V, i, H[i,i]));
  Reg = matdet(H);
  print("  regulator of Codex gens = ", Reg);
  R = ellrank(E, 2);
  RegR = matdet(ellheightmatrix(E, R[4]));
  print("  regulator of ellrank pts = ", RegR, "   ratio = ", Reg/RegR);
  \\ 2-saturation: image of <V> + E[2] in E(Q)/2E(Q) must have dim r+2
  T = [[0,0],[n,0],[-n,0]];
  gens = concat(V, T);
  M = matconcat(vector(#gens, i, vecF2(n, kum(n, gens[i]))));
  print("  Kummer image rank (need ", #V+2, ") = ", matrank(M));
  \\ odd primes up to B via ellsaturation
  forstep (B = 100, 2000, 1900,
    W = ellsaturation(E, V, B);
    RegW = matdet(ellheightmatrix(E, W));
    print("  ellsaturation(B=", B, "): Reg ratio old/new = ", Reg/RegW));
);
}
