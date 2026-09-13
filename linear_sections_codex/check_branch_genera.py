"""Riemann-Hurwitz from parity vectors; independent canonical-degree check."""
import json
from pathlib import Path

def genus_from_parities(r, branch_vectors):
    degree=2**r
    # In characteristic zero over an algebraic closure, each nonzero valuation
    # parity vector gives order-two inertia. Several odd coordinates do NOT
    # increase its order: the inertia element is their single combined vector.
    contributions=[]
    for v in branch_vectors:
        assert len(v)==r and all(x in (0,1) for x in v)
        inertia=2 if any(v) else 1
        contributions.append(degree-degree//inertia)
    two_g_minus_two=-2*degree+sum(contributions)
    assert two_g_minus_two%2==0
    return {'degree':degree,'parity_vectors':branch_vectors,'ramification_contributions':contributions,
            'total_ramification':sum(contributions),'genus':1+two_g_minus_two//2}

unit=lambda r,i:[int(j==i) for j in range(r)]
full=genus_from_parities(8,[unit(8,i) for i in range(8)]+[[1]*8])
pairs=genus_from_parities(4,[unit(4,i) for i in range(4) for _ in range(2)]+[[0]*4])
product=genus_from_parities(1,[[1]]*8+[[0]])
three=genus_from_parities(3,[unit(3,i) for i in range(3)]+[[1]*3])
assert [v['genus'] for v in (full,pairs,product,three)]==[321,17,3,1]
# Smooth complete intersection of seven quadrics in P^8: degree 128,
# omega=O(14-9)=O(5), independently giving degree(K)=640 and genus 321.
assert 2*full['genus']-2==(2*7-9)*2**7
out={'assumptions':'Characteristic zero, eight distinct nonzero linear coefficients; branch places and normalization taken geometrically.',
     'eight_individual_roots':full,'four_pair_products':pairs,'total_product':product,'three_individual_roots':three,
     'canonical_degree_check_full':{'degree':128,'canonical_O_degree':5,'canonical_divisor_degree':640,'genus':321}}
Path(__file__).with_name('branch_genera_checked.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v['genus'] for k,v in out.items() if isinstance(v,dict) and 'genus' in v})
