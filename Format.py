#! /home/yyin/software/python/bin/python3
import subprocess
import argparse

def processline(line,dic,order,cov):
    if line[0]=='#': return
    t = line.strip().split()
    key=(t[0],t[1])
    if len(t)>4 and int(t[4])<cov: return
    if not key in dic:
        array=[]
        for i in range(order):
            array.append('-')
        array.append(t[3])
        dic[key]=array
    else:
        dic[key].append(t[3])

def formdata(files,cov=0,bedfile=''):
    dic={}
    i=0
    for file in files:
        if bedfile=='':
            with open(file) as f:
                for line in f:
                    processline(line,dic,i,cov)
        else:
            p = subprocess.Popen('bedtools intersect -a %s -b %s' %(file,bedfile),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            for line in p.stdout:
                processline(line,dic,i,cov)
        i+=1
    result=[]
    for key in dic:
        chr,x=key
        if len(dic[key])==len(files):
            a=[chr,x]
            a.extend(dic[key])
            result.append(a)
    return result

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--inputfile',nargs='+',help="input file name", metavar="FILE")
    parser.add_argument('-l','--label',nargs='+')
    parser.add_argument('-c','--cov',default=0)
    args = parser.parse_args()
    data = formdata(args.inputfile,cov=int(args.cov))
    header = 'chr\tpos'
    for l in args.label:
        header=header+'\t'+l
    print(header)
    import sys
#    sys.exit()
    ans = [header+'\n']
#    print(data[:10])
    for d in data:
        s=''
        for i in range(len(d)-1):
            s=s+str(d[i])+'\t'
        s=s+str(d[-1])+'\n'
        ans.append(s)
    with open('meth.bed','w') as f:
        f.writelines(ans)


