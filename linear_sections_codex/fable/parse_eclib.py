import re,sys
from criterion import criterion
txt=open('validate_eclib.log').read()
chunks=re.split(r'###PAIR ',txt)[1:]
checked=0; bad=[]; noblock=0
for ch in chunks:
    head=ch.split('\n',1)[0]; m,n=map(int,head.split()[:2]); ep=head.split('EPRIME')[1].strip()
    blocks=re.split(r'\* Using 2-isogeny number \d+ \*',ch)[1:]
    found=None
    for b in blocks:
        mm=re.search(r'minimal model \[([-\d,]+)\]',b)
        if mm and mm.group(1)==ep:
            a=re.search(r"rk\(S\^\{phi\}\(E'\)\)=\s*(\d+)",b); c=re.search(r"rk\(S\^\{phi'\}\(E\)\)=\s*(\d+)",b)
            if a and c: found=(int(a.group(1)),int(c.group(1)))
    if found is None: noblock+=1; continue
    cr=criterion(m,n); mine=(len(cr['S_phi']).bit_length()-1,len(cr['S_phihat']).bit_length()-1)
    checked+=1
    if set(found)!=set(mine): bad.append((m,n,'eclib',found,'criterion',mine))
print('eclib vs критерий: проверено',checked,'без блока',noblock,'расхождений',len(bad)); print(bad[:10])
