import sys
import multiprocessing

class Pshell():

    def __init__(self,commend):
        self.commend = commend
        self.out = ""

    def change(self,commend):
        self.commend = commend

    def getcommend(self,commend):
        return commend

    def process(self):
        import subprocess
        print(self.commend)
        t = subprocess.Popen(self.commend,shell=True,stdout = subprocess.PIPE,stderr = subprocess.PIPE)
        self.out = t.stderr.read().decode()
        print(self.out)

    def get_result(self):
        return self.out

class reads():

    def __init__(self,fqname,chromosome,startpos,order,mismatch):
        self.fqname = fqname
        self.chromosome = chromosome
        self.startpos = startpos
        self.sum = 0
        self.order = order
        self.chain = 0
        self.mismatch = mismatch

    def equal(self,rreads):
        if (rreads.getChrom()==self.chromosome)and(rreads.getStartPos()==self.startpos): return True
        return False

    def getChrom(self):
        return self.chromosome

    def getOrder(self):
        return self.order

    def getSum(self):
        return self.sum

    def getStartPos(self):
        return self.startpos

    def join(self,tail,mis):
        self.mismatch = self.mismatch + mis
        self.sum = tail.getOrder()-self.getOrder()

    def getChain(self):
        return self.chain

    def getMismatch(self):
        return self.mismatch

    def canjoin(self,tail,step,xx,binsize):
        if (tail.getChrom()==self.chromosome):
            if (tail.getStartPos() == self.startpos + step*(self.order-tail.getOrder()) or tail.getStartPos() == self.startpos - step*(self.order-tail.getOrder())):
            #if tail.getStartPos()>self.startpos: self.chain = 0
            #else: self.chain = 1
                return True
            if xx!=binsize:
                if xx % step == abs(self.startpos-tail.getStartPos())%step:
                    return True
        return False


def combine(readsname,fastqset,bamitem,startfileorder,step):
    fastq = fastqset[readsname] #0:info 1:reads 2:info 3:quality
    order = bamitem.getOrder()
    readstr=fastq[1][order]
    qualitystr=fastq[3][order]
    for i in range(order+1,order+bamitem.getSum()+1):
        readstrtail = fastq[1][i][(-1)*step:]
        qualitytail = fastq[3][i][(-1)*step:]
        if i==order+bamitem.getSum():
            readstrtail = fastq[1][i][len(fastq[1][0])-step:]
            qualitytail = fastq[3][i][len(fastq[3][0])-step:]
        readstr=readstr+readstrtail
        qualitystr = qualitystr+qualitytail
        #print(readstr,qualitystr)
 #   print(readsname)
    return fastq[0]+'\n'+readstr+'\n'+fastq[2]+'\n'+qualitystr+'\n'

def overlap(a,b):
    #if a.getChrom()!=b.getChrom(): return False
    ordera = a.getOrder()
    orderb = b.getOrder()
    if (orderb<ordera):
        k=a
        a=b
        b=k
    ordera = a.getOrder()
    orderb = b.getOrder()
    sa = a.getStartPos()
    sb = b.getStartPos()
    suma = a.getSum()
    sumb = b.getSum()
    #if sa+suma*5+30>sb: return True
    if orderb-ordera-suma<=5: return True
    return False

def cut(step,bin,reads):
    start = 0
    end = 0
    result=[]
    l = len(reads)
    while end<l:
        end = start+bin
        if end>l: end = l
        result.append(reads[start:end])
        start += step
    return result



def do_process(l):
    #print(l+'in')
    temp = l.strip().split()
    length = len(temp)
    if length<=0 or length>2:
        print("Parameter error in "+l)
        sys.exit()
    outputname = temp[0][:len(temp[0])-6]
    #print(outputname)
    if length==2 :
        commend='bsmap -a '+temp[0]+' -b '+temp[1]+'  -d /data/dsun/ref/humanigenome/hg19.fa -o '+outputname+'.bam -n 1 -q 3 -r 0'
    else:
        commend='bsmap -a '+temp[0]+' -d /data/dsun/ref/humanigenome/hg19.fa  -o '+outputname+'.bam -n 1 -q 3 -r 0'
    First_try = Pshell(commend)
    First_try.process()

#Test1 done
    inputfileinfo=l.strip().split()
    commend = 'samtools view '+outputname+'.bam > '+outputname+'.sam'
    BamFileReader = Pshell(commend)
    BamFileReader.process()
    with open(outputname+".sam") as sam:
    #second column in sam file: 64, mate 1; 128, mate 2;
        samlines = sam.readlines()
    set_sam = {}
    for line in samlines:
        temp = line.strip().split()
        m1 = (int(temp[1]) & 64)
        m2 = (int(temp[1]) & 128)
#        print(temp[1],m1,m2)
        if m1>0: mate = 1
        elif m2>0: mate = 2
        else: mate = 0
        if temp[0] in set_sam: set_sam[temp[0]]=3
        else: set_sam[temp[0]]=mate
#        print(mate)
#    for k in set_sam:
#        print(k)
#        break

    UnmappedReads = {}
    o=0
    for filename in inputfileinfo:
        o+=1
        with open(filename) as f:
            fastqlines = f.readlines()
        for i in range(len(fastqlines)):
            fastqlines[i] = fastqlines[i].strip()
        for i in range(0,len(fastqlines),4):
            line1 = fastqlines[i].split()
            line1[0] = line1[0]
 #           print(line1[0][1:])
            if (line1[0][1:] in set_sam):
                #print(line1[0][1:],set_sam[line1[0][1:]],int(line1[1][0]))
                if set_sam[line1[0][1:]]==0 or set_sam[line1[0][1:]]==3 : continue
                if set_sam[line1[0][1:]]==int(line1[1][0]) : continue
            fastqlines[i] = line1[0]+'_'+line1[1][0]+' '+line1[1]
            UnmappedReads[line1[0]+'_'+line1[1][0]]=fastqlines[i:i+4]

#We've got a dictionary named UnmappedReads = {readsname:[line1,line2,line3,line4]}

    step = 5
    length_bin = 30#30
    max_length = 0
    for readsname in UnmappedReads:
        result = cut(step,length_bin,UnmappedReads[readsname][1])
        if len(result)>max_length: max_length = len(result)
        UnmappedReads[readsname][1]=result
        UnmappedReads[readsname][3]=cut(step,length_bin,UnmappedReads[readsname][3])

    Part_Fastq_Filename = []
    for i in range(max_length):
        filecontent = []
        for readsname in UnmappedReads:
            link = UnmappedReads[readsname]
            if len(link[1])<=i: continue
            filecontent.append(link[0]+'\n'+link[1][i]+'\n'+link[2]+'\n'+link[3][i]+'\n')
        name = outputname+'.part'+str(i+1)+'.fastq'
        Part_Fastq_Filename.append(name)
        with open(name,'w') as f:
            f.writelines(filecontent)

    #We've got the splited fastq file, filename is stored in Part_Fastq_Filename
   # p = multiprocessing.Pool(processes=7)
    for i in range(len(Part_Fastq_Filename)):
        commend = 'bsmap -a '+Part_Fastq_Filename[i]+'  -d /data/dsun/ref/humanigenome/hg19.fa  -o '+Part_Fastq_Filename[i]+'.bam -n 1 -q 3 -r 0 -R'
        Batch_try = Pshell(commend)
        Batch_try.process()
   # p.close()
   # p.join
   #run bsmap and get bam files named as Part_Fastq_Filename[i].bam

    result = {}
    file_order = 0
    for name in Part_Fastq_Filename:
        commend = 'samtools view '+name+'.bam> '+name+'.sam'
        SamFileMaker = Pshell(commend)
        SamFileMaker.process()
        #print(name)
        with open(name+'.sam') as f:
            partsamlines = f.readlines()
        for line in partsamlines:
            s = line.strip().split('\t')
            mismatch = int(s[11][s[11].rfind(':')+1:])
            #s[0] = s[0]
            if not(s[0] in result) or (len(result[s[0]])==0):
                result[s[0]]=[reads(s[0],s[2],int(s[3]),file_order,mismatch)]

            else:
               # print(1)
                temp = reads(s[0],s[2],int(s[3]),file_order,mismatch)
                #reads from cliped mapped bam
                join_or_not=False
                i=0

                for ss in result[s[0]]:
                    #distance=length_bin-abs(ss.getOrder()+ss.getSum()-file_order)*step
                    if True:#distance>step:
                        if ss.canjoin(temp,step,len(s[9]),length_bin):
                            refseq = s[12][s[12].rfind(':')+1:][2:-2][-1*step:].upper()
                            readsseq = s[9][-step:].upper()
                            mis=0
                            for ppp in range(step):
                                if refseq[ppp]!=readsseq[ppp]: mis+=1
                            ss.join(temp,mis)
                            join_or_not=True
                            break
                        #else:
                            #if distance<=step, we just need to pick the smallest mismatch
                            #new reads should be uniquely mapped
                        #    if mismatch<=ss.getMismatch():
                        #        if ss.getSum()==0: result[s[0]][i]=temp
                        #        join_or_not=True
                        #        break
                    i+=1

                if not join_or_not:
                    #head reads should be uniquely mapped
                    result[s[0]].append(temp)

        file_order+=1

    #join done
    for name in result:
        nonjoin_num=0
        for item in result[name]:
            if item.getSum()==0 and item.getMismatch()>1:
                result[name].remove(item)
    for name in result:
        num = len(result[name])
        del_mark = [0 for i in range(num)]
        for i in range(num):
            for j in range(i+1,num):
                #if name=='NS500669:63:H2NW2BGX2:1:13302:19288:9216_1':
                #    print(overlap(result[name][i],result[name][j]))
                if overlap(result[name][i],result[name][j]):
                    sss = result[name][i].getSum()-result[name][j].getSum()
                    if sss>0: del_mark[j]=1
                    elif sss<0: del_mark[i]=1
                    else:
                        mis = result[name][i].getMismatch()-result[name][j].getMismatch()
                        if mis>0: del_mark[i]=1
                        if mis<=0: del_mark[j]=1
       # if name=='NS500669:63:H2NW2BGX2:1:13302:19288:9216_1':
       #     print(del_mark)
       #     print(result[name])
       #     print(result[name][0].getSum(),result[name][1].getSum())
       #     print(result[name][0].getOrder(),result[name][1].getOrder())
        for i in range(num-1,-1,-1):
            if del_mark[i]==1:
                del result[name][i]

    finalfastq = []
    for name in result:
        for read in result[name]:
            finalfastq.append(combine('@'+name,UnmappedReads,read,read.getOrder(),step))

    #form a final new fastq file and do the final mapping
    with open(outputname+'_finalfastq.fastq','w') as f:
        f.writelines(finalfastq)

    commend = 'bsmap -a '+outputname+'_finalfastq.fastq -d /data/dsun/ref/humanigenome/hg19.fa  -o '+outputname+'_split.bam -n 1 -q 3 -r 0'
    Bam = Pshell(commend)
    Bam.process()
    m=Pshell('samtools merge '+outputname+'_combine.bam '+outputname+'.bam '+outputname+'_split.bam')
    m.process()
    print("Merge done!\nCreated final bam file called "+outputname+'_combine.bam')

if __name__=="__main__":
    with open("config.txt") as f:
        lines = f.readlines()
    import multiprocessing
    pool = multiprocessing.Pool(processes=1)
    for l in lines:
        #pool.apply_async(do_process,(l,))
        do_process(l) #pass file name to do_process
    pool.close()
    pool.join()

