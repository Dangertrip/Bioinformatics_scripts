conversion_ratio=0.96
with open('/data/dsun/ref/mouseigenome/mm10.fa') as f:
    lines = f.readlines()
chr=''
ans={}
pos=0
mark=False
s=''
length=0
l=[]
chr_order=[]
for line in lines:
    temp = line.strip()
    if temp[0]=='>':
        if chr!='':
            ans[chr]=s
            ll=len(s)
            length+=ll
            l.append(ll)
            s=''
        chr= temp[1:]
        chr_order.append(chr)
        continue
    s+=temp
ans[chr]=s
ll=len(s)
l.append(ll)
length+=ll
import random
#print(ans)
def getcg(read):
    num=0
    if random.random()>0.5:
        for i in range(1,len(read)):
            if ('G'==read[i] or 'g'==read[i]) and ('C'==read[i-1] or 'c'==read[i-1]):
                num+=1
    else:
        for i in range(1,len(read)):
            if ('C'==read[i] or 'c'==read[i]) and ('G'==read[i-1] or 'g'==read[i-1]):
                num+=1

    return num

num=[]
print(l)
zero=0
for i in range(100000000):
    pos = random.randint(0,length)
    chr=0
    while pos>l[chr]:
        pos-=l[chr]
        chr+=1

    read=ans[chr_order[chr]][max(1,pos-1):min(pos+76,l[chr])]
    g=getcg(read)
    num.append(g)
    if g==0: zero+=1
print(zero/float(100000000))
print((100000000-zero)/float(100000000))

#import seaborn as sns
#from matplotlib import pyplot as plt
#import numpy as np
#n=np.array_split(num,10)
#for nn in n:
#    sns.distplot(nn)
#plt.savefig('125.pdf')


#
