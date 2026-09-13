# положительный контроль отображения: z=0 даёт точку (y=1) — должна быть среди точек кривой для любого S;
# отрицательный: для известного решения-заглушки проверка nine_ok на z из Бремнера-подобного примера не нужна; проверяем, что z=0 попадает в zvals при ранге 0
from g1census import *
import random
ok=0; tot=0
for sl in ['19/60','41/60','2/7','3/5']:
    r,s=map(int,sl.split('/')); vals=[r,s-r,s,s+r]; C=sorted(set(vals+[-v for v in vals]))
    for S in itertools.combinations(C,4):
        E,L,c0=weier([QQ(c) for c in S])
        if E.discriminant()==0: continue
        if E.pari_curve().ellrank()[1]==0:
            tot+=1; zs=zvals(E,L,c0)
            if 0 in zs: ok+=1
print('кривые ранга 0 с z=0 среди точек:',ok,'из',tot)
