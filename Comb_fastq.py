from utils import overlap,Pshell
import utils

def GetFastqList(joined_reads,Part_Fastq_Filename,step,length_bin):

    nameset={}
    contentset={}
    filenum = len(Part_Fastq_Filename)
    #Generate a dictionary which contains readsname, start file order and extend fraction number
    for name in joined_reads:
        reads_list = joined_reads[name]
        if len(reads_list)==0: continue
        #n = reads_list[0].fqname
        n = name
        nameset[n]=[[read.order,read.sum] for read in reads_list]
        contentset[n]=[['',''] for i in range(len(nameset[n]))]#read_content,read_quality


    file_order=0
    for name in Part_Fastq_Filename:
        with open(name) as f:
            lines = f.readlines()
        for i in range(0,len(lines),4):
            fqname = lines[i].strip().split()[0]
            if not fqname in nameset: continue
            read = lines[i+1].strip()
            quality = lines[i+3].strip()
            fraction_num=0
            for order,sum in nameset[fqname]:
                if file_order>=order and file_order<=order+sum:
                    add_length = ((len(read)-1)%step)+1
                    contentset[fqname][fraction][0]+=read[-1*add_length:]
                    contentset[fqname][fraction][1]+=quality[-1*add_length:]
                fraction_num+=1
    return contentset


def combine(outputname,Part_Fastq_Filename,step,length_bin):
    cache_length=3
    result={}
    file_order=0
    for name in Part_Fastq_Filename:
        command = 'samtools view '+name+'.bam'
        SamFileMaker = Pshell(command)
        SamFileMaker.process()
        partsamlines = SamFileMaker.out.split('\n')[:-1]
        del SamFileMaker
        #print(partsamlines[-1])
        for line in partsamlines:
            #print(line)
            s = line.strip().split('\t')
            mismatch = int(s[11][s[11].rfind(':')+1:])
            if not(s[0] in result) or (len(result[s[0]])==0):
                result[s[0]]=[utils.reads(s[2],int(s[3]),file_order,mismatch)]
            else:
                temp = utils.reads(s[2],int(s[3]),file_order,mismatch)
                #reads from cliped mapped bam
                join_or_not=False
                read_length = len(s[9])
                tail_length = ((read_length-1)%step)+1
                refseq = s[12][-2-tail_length:-2]
                readsseq = s[9][-step:]
                strand = s[13][-2:]
                for reads in result[s[0]]:#Try to join existing seeds
                    if reads.canjoin(temp,step,read_length,length_bin):
                        mis=0
                        for ppp in range(step):
                            if refseq[ppp]!=readsseq[ppp]:
                                #Here ++/+-/-+/-- should be considered. C/T or A/G match should be identified.
                                if strand[0]=='+':
                                    if (refseq[ppp]=='C' and readsseq=='T'): continue
                                else:
                                    if (refseq[ppp]=='G' and readsseq=='A'): continue
                                mis+=1
                        reads.join(temp,mis)
                        join_or_not=True
                        break

                if not join_or_not: #temp reads haven't join any exist reads
                    result[s[0]].append(temp) #add temp reads to array as new seed
                    if len(result[s[0]]>=s):
                        for read in result[s[0]]:
                            if read
                    print(len(result[s[0]]))
        file_order+=1
    #join done
    #filter results: filter1
    for name in result:
        nonjoin_num=0
        reads_list=result[name]
        for i in range(len(reads_list)-1,-1,-1):
            if reads_list[i].getSum()==0 and reads_list[i].getMismatch()>1:
                reads_list.pop(i)#Remove all reads which have more than 1 mistake and never be joined
    #filter results: filter2
    for name in result:
        reads_list = result[name]
        num = len(reads_list)
        del_mark = [0 for i in range(num)]
        for i in range(num):
            for j in range(i+1,num):
                if overlap(result[name][i],result[name][j]):
                    sss = result[name][i].getSum()-result[name][j].getSum()
                    if sss>0: del_mark[j]=1
                    elif sss<0: del_mark[i]=1
                    else:
                        mis = result[name][i].getMismatch()-result[name][j].getMismatch()
                        if mis>0: del_mark[i]=1
                        else: del_mark[j]=1
        #Only keep the best read which has the most extends and the least mismatches.
        for i in range(num-1,-1,-1):
            if del_mark[i]==1:
                reads_list.pop(i)

    fastq_dic = GetFastqList(result,Part_Fastq_Filename,step,length_bin)
    with open(outputname+'_finalfastq.fastq','w') as f:
        for name in fastq_dic:
            read,quality = fastq_dic[name]
            f.write(name+'\n')
            f.write(read+'\n')
            f.write('+\n')
            f.write(quality+'\n')







