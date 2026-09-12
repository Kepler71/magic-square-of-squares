# Генераторы над квадратичными полями: P_sqrt6 над Q(√6), P_i над Q(i). Точная проверка и индекс.
import functools
print = functools.partial(print, flush=True)
src = open('mw_galois.sage').read()
exec(preparse(src.split('# --- минимальные сечения')[0]))          # K, E, A4, A6
exec(preparse(src.split('# --- высоты (как в heights_modp.sage)')[1].split('T  = E(')[0]))   # h, pair
T  = E(Fs(36*s^2 - 240*s + 144), Fs(0))
S1 = E(Fs(72*s^2 - 96*s + 288), Fs(432*(s+2)*(s^2-12)))
X6 = 36*s^2 - 240*s + 1680; Y6 = 1536*S3*S2*(20 - 3*s)
Xi = 36*s^2 - 240*s + 144 - 576*I_; Yi = -(3456 + 6912*I_)*s - 6912 + 13824*I_
print("P_sqrt6 on surface:", Y6^2 == X6^3 + A4*X6 + A6, "  P_i on surface:", Yi^2 == Xi^3 + A4*Xi + A6)
P6, Pi = E(Fs(X6), Fs(Y6)), E(Fs(Xi), Fs(Yi))
G = matrix(QQ, [[pair(x, y) for y in (S1, Pi, P6)] for x in (S1, Pi, P6)])
print("Gram <S1, P_i, P_sqrt6>:\n", G, "\ndet =", G.det(), "  index in MWL (det 1/4): sqrt(det/(1/4)) =", sqrt(G.det()/(QQ(1)/4)))
print("h(T+P_i) =", h(T + Pi), " h(P_i + P_sqrt6) =", h(Pi + P6), " h(S1 - P_i) =", h(S1 - Pi))
