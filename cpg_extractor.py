with open('/data/dsun/ref/mouseigenome/mm10.fa') as f:
    lines = f.readlines()
chr=''
ans=[]
pos=0
mark=False
for line in lines:
    temp = line.strip()
    if mark and (temp[0]=='g' or temp[0]=='G'):
        t = chr+'\t'+str(pos)+'\t'+str(pos+1)+'\t1\n'
        ans.append(t)
    if temp[0]=='>':
        chr= temp[1:]
        pos=0
        continue
    mark=False
    for i in range(len(temp)-1):
        pos+=1
        if (temp[i]=='c' or temp[i]=='C') and (temp[i+1]=='g' or temp[i+1]=='G'):
            t = chr+'\t'+str(pos)+'\t'+str(pos+1)+'\t1\n'
            ans.append(t)
    if temp[-1]=='c' or temp[-1]=='C':
        mark=True
    pos+=1

with open('mm10_cpg.bed','w') as f:
    f.writelines(ans)
