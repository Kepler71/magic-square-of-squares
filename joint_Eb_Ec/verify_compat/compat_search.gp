\\ Дополнительно: пары (b,c) с НЕизоморфными E_b, E_c (b/c не квадрат), где
\\ доказанные нижние границы ellrank дают rank E_b>=2, rank E_c>=2, rank E_{b+c}>=1, rank E_{c-b}>=1,
\\ и верхние границы совпадают (ранг точен).
default(parisize, 10^9);
cache = Map();
rk(t) = {my(k = core(t), R); if (mapisdefined(cache, k, &R), return(R)); R = ellrank(ellinit([0,0,0,-k^2,0]), 1)[1..2]; mapput(cache, k, R); R;}
{
cnt = 0;
for (c = 2, 400, for (b = 1, c - 1,
  if (gcd(b, c) > 1 || issquare(b*c), next);
  Rb = rk(b); if (Rb[1] < 2, next);
  Rc = rk(c); if (Rc[1] < 2, next);
  Rp = rk(b + c); if (Rp[1] < 1, next);
  Rm = rk(c - b); if (Rm[1] < 1, next);
  cnt++;
  if (cnt <= 12, print("b=", b, " c=", c, "  ranks [lo,hi]: E_b ", Rb, " E_c ", Rc, " E_{b+c} ", Rp, " E_{c-b} ", Rm));
  if (cnt >= 12, break(2))));
print("found ", cnt);
}
