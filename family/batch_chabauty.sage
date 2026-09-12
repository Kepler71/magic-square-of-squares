# Для списка сечений: класс с известной точкой -> выбор расщеплённого p с малым N -> chabauty.sage;
# плюс строка Magma для Richelot-двойственной C1 (ранговая граница — у Codex).
import subprocess, sys, re
def sqf(n):
    n = ZZ(n); return sign(n) * prod(pr^(e % 2) for pr, e in n.abs().factor())
# (b,h,n, d3 положительного класса, X0 его точки) — из family/sections_batch.log
secs = [(7,1,5,15,2), (23,7,17,34,QQ(5)/3), (31,17,25,7,QQ(4)/3), (41,1,29,609,QQ(5)/2),
        (47,23,37,1295,6), (49,31,41,41,QQ(5)/4), (73,17,53,265,QQ(7)/2)]
magma = []
for (b, h, n, D, X0) in secs:
    A, C = ZZ((h^2 + n^2)/2), ZZ((b^2 + n^2)/2)
    S2 = ZZ((A^2 + C^2)/2) if (A^2 + C^2) % 2 == 0 else None
    magma.append(f"f:={D}*(t^2-1)*(t^2+1)*({S2}*t^2-{2*A*C}*t+{S2}); print \"({b},{h},{n})\", RankBounds(Jacobian(HyperellipticCurve(f)));")
    t0 = QQ(X0)^2
    dk = sqf((C - A)*(C + A))
    # выбор p: расщеплено в Q(sqrt(dk)), минимальный N — пробуем по возрастанию, берём первый с N <= 12
    chosen = None
    for p in primes(7, 200):
        if kronecker(dk, p) != 1 or (2*D*A*C*(C^2 - A^2)) % p == 0: continue
        out = subprocess.run(['sage', 'chabauty.sage', str(b), str(h), str(n), str(D), str(t0), str(p)],
                             capture_output=True, text=True, timeout=1800).stdout
        m = re.search(r"N = (\d+)", out)
        if not m: continue
        N = int(m.group(1))
        if N > 12: continue
        chosen = (p, N, out); break
    if chosen is None:
        print(f"({b},{h},{n}): no suitable p found"); continue
    p, N, out = chosen
    lines = [l for l in out.splitlines() if re.search(r"local saturation|classes:|rigorous|TOTAL|known P|exact zero|Strassmann bound", l)]
    print(f"######## ({b},{h},{n}) D={D} X0={X0} p={p} N={N}")
    print("\n".join(lines))
print("\n==== Magma lines (Q:=Rationals(); R<t>:=PolynomialRing(Q);) ====")
print("\n".join(magma))
