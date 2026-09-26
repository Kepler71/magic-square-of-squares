# Однородная форма: центр - тоже неизвестный класс. Проверка: 8 линий => центр - квадрат (ранг 7).
from itertools import product
cells=[(i,j) for i in (-1,0,1) for j in (-1,0,1)]
lines=[[(i,-1),(i,0),(i,1)] for i in (-1,0,1)]+[[(-1,j),(0,j),(1,j)] for j in (-1,0,1)]
lines+=[[(-1,-1),(0,0),(1,1)],[(-1,1),(0,0),(1,-1)]]
sols=[x for x in (dict(zip(cells,b)) for b in product((0,1),repeat=9))
      if all(sum(x[c] for c in ln)%2==0 for ln in lines)]
print('решений:',len(sols),'=> ранг',9-(len(sols).bit_length()-1),'; центр всегда 0:',all(x[(0,0)]==0 for x in sols))
