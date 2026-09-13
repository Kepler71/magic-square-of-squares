import json
pr=json.load(open('g1_odd.json'))['m_le_20']
def EB(m,n): return EllipticCurve(QQ,[0,-((m^2-n^2)^2+(m^2+n^2)^2),0,(m^2-n^2)^2*(m^2+n^2)^2,0])
bad=[(m,n) for m,n in pr if int(pari(EB(m,n).minimal_model().a_invariants()).ellinit().ellrank()[1])!=0]
print("контроль: пары m,n ≤ 20, где критерий доказал ранг 0, а PARI для B видит ранг > 0:", bad if bad else "нет")
