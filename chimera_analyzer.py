import sys
filename = sys.argv[1]
with open(filename) as f:
    lines = f.readlines()
distance=[]
for i in range(0,len(lines),2):
    line = lines[i]
    t = line.strip().split()
    achr,astart,aend,_,_,astrand = t
    line = lines[i+1]
    t = line.strip().split()
    bchr,bstart,bend,_,_,bstrand = t
    if achr!=bchr: continue
    if astart<bstart:
        dis = bstart-aend
    else:
        dis = astart-bend
    distance.append(dis)

from matplotlib import pyplot as plt
import seaborn as sns
sns.kdeplot(distance)
plt.show()


