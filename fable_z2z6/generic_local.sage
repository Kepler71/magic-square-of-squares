# Общие леммы при p=2,3 для y^2 = x^3 + (alpha x + beta)^2 с произвольными alpha, beta.
import random
load('/home/kep/magicKube/fable_z2z6/three_isog.sage')
random.seed(int(7))
def certified_image(alpha, beta, p):
    E = EllipticCurve(QQ, [0, alpha**2, 0, 2*alpha*beta, beta**2])
    T = E(0, beta); assert 3*T == E(0)
    Ep = E.isogeny(T).codomain()
    tang = tangent_Tprime(Ep)
    n = h1_mu3_dim(p); extra = 0
    while True:
        G = local_image_E(E, alpha, beta, p, n, extra)
        Gp = local_image_Ep(Ep, tang, p, n, extra)
        if G.dim() + Gp.dim() == n: return G
        assert G.dim() + Gp.dim() < n
        extra += 2000
        if extra > 12000: return None
from collections import Counter, defaultdict
tab = defaultdict(Counter)
for trial in range(300):
    # p = 2: alpha нечётное, beta = 2^k * нечётное
    k = int(random.choice([1,2,3,4,5,6,7,8]))
    alpha = int(random.randrange(-200, 200)) | 1
    beta = (2**k) * (int(random.randrange(-200, 200)) | 1)
    G = certified_image(alpha, beta, 2)
    tab[('p=2', 'v2(beta)=%d' % k, 'alpha mod 8=%d' % (alpha % 8))][G.dim() if G else 'uncert'] += 1
for trial in range(300):
    # p = 3: (i) 3 ∤ alpha, v3(beta)=k>=1 ; (ii) v3(alpha)=j>=1, 3 ∤ beta
    if random.random() < 0.5:
        k = int(random.choice([1,2,3,4,5,6]))
        alpha = int(random.randrange(-200, 200)); 
        while alpha % 3 == 0: alpha += 1
        beta = 3**k * int(random.choice([x for x in range(-60,60) if x % 3]))
        G = certified_image(alpha, beta, 3)
        tab[('p=3', 'v3(beta)=%d' % k, 'alpha mod 3=%d' % (alpha % 3))][G.dim() if G else 'uncert'] += 1
    else:
        j = int(random.choice([1,1,2,3]))
        alpha = 3**j * int(random.choice([x for x in range(-60,60) if x % 3]))
        beta = int(random.choice([x for x in range(-200,200) if x % 3]))
        G = certified_image(alpha, beta, 3)
        tab[('p=3', 'v3(alpha)=%d' % j, 'beta mod 9=%d' % (beta % 9))][G.dim() if G else 'uncert'] += 1
for key in sorted(tab, key=str): print(key, dict(tab[key]))
