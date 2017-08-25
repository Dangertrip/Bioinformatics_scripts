import multiprocessing
from Comb_fastq import combine
from utils import *
import sys
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
    #First_try.process()

#Test1 done
    inputfileinfo=l.strip().split()
    commend = 'samtools view '+outputname+'.bam > '+outputname+'.sam'
    BamFileReader = Pshell(commend)
    #BamFileReader.process()
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
        '''
        with open(filename) as f:
            while 1:
                line1 = f.readline()
                if not line1:
                    break
                line1 = line1.strip().split()
                line2 = f.readline().strip()
                line3 = f.readline()
                line4 = f.readline().strip()
                line1[0] = line1[0]
 #           print(line1[0][1:])
                if (line1[0][1:] in set_sam):
                #print(line1[0][1:],set_sam[line1[0][1:]],int(line1[1][0]))
                    if set_sam[line1[0][1:]]==0 or set_sam[line1[0][1:]]==3 : continue
                    if set_sam[line1[0][1:]]==int(line1[1][0]) : continue

                temp = line1[0]
                if length>1: temp+='_'+line1[1][0]
                #Maybe the mate search method is buggy. Cuz there are different structures of reads name generated by different sequencing machine.
                #fastqlines[i] = line1[0]+'_'+line1[1][0]+' '+line1[1]
                UnmappedReads[temp]=[line2,line4]
'''
#We've got a dictionary named UnmappedReads = {readsname:[line1,line2,line3,line4]}

    step = 5
    length_bin = 30#30
    max_length = 10
#Change cut funtion into cut(setp,length_bin,string,fileorder), return Available(T/F), reads_fraction

    Part_Fastq_Filename = []
    for i in range(max_length):
        filecontent = []
        #for readsname in UnmappedReads:
        #    link = UnmappedReads[readsname]
            #if len(link[1])<=i: continue
        #    mark,cutreads = cut(step,length_bin,link[0],i)
        #    if not mark: continue
        #    _,cutquality = cut(step,length_bin,link[1],i)
        #    filecontent.append(readsname+'\n'+cutreads+'\n+\n'+cutquality+'\n')
        #if len(filecontent)==0: break
        name = outputname+'.part'+str(i+1)+'.fastq'
        Part_Fastq_Filename.append(name)
        #with open(name,'w') as f:
        #    f.writelines(filecontent)
    del UnmappedReads
    #We've got the splited fastq file, filename is stored in Part_Fastq_Filename
   # p = multiprocessing.Pool(processes=7)
    for i in range(len(Part_Fastq_Filename)):
        commend = 'bsmap -a '+Part_Fastq_Filename[i]+'  -d /data/dsun/ref/humanigenome/hg19.fa  -o '+Part_Fastq_Filename[i]+'.bam -n 1 -q 3 -r 0 -R'
        Batch_try = Pshell(commend)
        #Batch_try.process()

   #run bsmap and get bam files named as Part_Fastq_Filename[i].bam
    #import combine to generate the finalfastq
    combine(outputname,Part_Fastq_Filename,step,length_bin)

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

