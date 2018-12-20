
import sys,re,csv,itertools
from operator import *
from itertools import *
import numpy
from Bio import SeqIO

args=sys.argv

if len(args)<2:
    snpsFile="/media/Second3TB/Work/WorkVLA/Data/Martin/WGS_Results/20130606_RD0057__Project_RD0057/AF2122/Mapped/12-0023-05/12-0023-05.pileup_SN.csv"
    gbkFile="/media/Second3TB/Work/WorkVLA/Software/references/AF2122.gbk"
    fnaFile="/media/Second3TB/Work/WorkVLA/Software/references/AF2122.fna"
else:    
    snpsFile=sys.argv[1]
    gbkFile=sys.argv[2]
    fnaFile=sys.argv[3]
    
def readTable(fname,ch,dig=0):
    infile=open(fname,"rb")
    data = csv.reader(infile, delimiter=ch)
    if dig!=0:
        dataOut = [map(strToint,row) for row in data]
    else:
        dataOut = [row for row in data]
    infile.close()
    print str(len(dataOut))+" X "+str(len(dataOut[0]))+" Matrix in "+fname
    return(dataOut)
    
def strToint(s):
    if all(s[i] in ".0123456789" for i in range(len(s))):
        if "." in s:
            return(float(s))
        if s.isdigit():
            return int(s)
    else:
        return s

def readfnaFileSeveralSequences(fname):
    fileIn = open(fname, 'rb')
    lines= fileIn.readlines()
    ids=[]
    seqs=[]
    seq=""
    for line in lines:
        if line[0]==">":
            if len(ids)>0:
                seqs.append(seq)
                seq=""
            ids.append(line[1:-1].strip())
        else:
            seq=seq+line[:-1].strip()
    seqs.append(seq)
    print str(len(ids))+" sequences in "+fname
    return ids,seqs

def writeCSV(fname,matrix):
    with open(fname, "wb") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print "file "+fname+" saved."

def isACGT(c):
    return c in ["a","A","c","C","g","G","t","T"]
    
def complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
    letters = list(s)
    new=[]
    for l in letters:
        if isACGT(l):
            new=new+[basecomplement[l]]
        else:
            new=new+[l]
    return ''.join(new)

def revcom(s):
    return complement(s[::-1])

def translate_dna(sequence,pos=-1):

    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    #proteinsequence = ''
    #start = sequence.find('ATG')
    #sequencestart = sequence[int(start):]
    #stop = sequencestart.find('TAA')
    #cds = str(sequencestart[:int(stop)+3])

    #for n in range(0,len(cds),3):
    #    if cds[n:n+3] in codontable:
    #        proteinsequence += codontable[cds[n:n+3]]
    #        print proteinsequence
    #    sequence = ''
    codonInfo="NA"
    proteinsequence = ''
    for n in range(0,len(sequence),3):
        if sequence[n:n+3] in codontable:
            proteinsequence += codontable[sequence[n:n+3]]
            if pos>-1 and pos in range(n,n+3):
                codonInfo=[sequence[n:pos]+"("+sequence[pos]+")"+sequence[pos+1:n+3],codontable[sequence[n:n+3]]]
        #else:
        #    print "not found "+sequence[n:n+3]+" at position "+str(n)
    return proteinsequence,codonInfo

    
def gsAnnot(x,seqref):
    refSeq=seqref
    geneSeq=seqref[:(x[ppos]-1)]+x[palt]+seqref[x[ppos]:]
    
    
    refPro,refCodon=translate_dna(refSeq,x[ppos]-1)
    genePro,geneCodon=translate_dna(geneSeq,x[ppos]-1)
    if refPro==genePro:
        annot="syn"
    else:
        annot="non"    
    return x+[refCodon,geneCodon,annot]


refid,refseq=readfnaFileSeveralSequences(fnaFile)

snps=readTable(snpsFile,",",1)
prefName=0
ppos=1
pref=2
palt=3

for gbk_record in SeqIO.parse(open(gbkFile,"r"), "genbank"):
        print "Name %s, %i features" % (gbk_record.name, len(gbk_record.features))

newtab=[]

for snp in snps:
    fea=[feature for feature in gbk_record.features if snp[ppos]>=feature.location.start.position and snp[ppos]<=feature.location.end.position and feature.type=='CDS']
    if len(fea)>0:    
        fea=fea[0]    
        gene="none"
        if "gene" in fea.qualifiers:
            gene=fea.qualifiers['gene'][0]
        db_xref="none"
        if 'db_xref' in fea.qualifiers:
            db_xref=fea.qualifiers['db_xref'][0]
        product="none"
        if 'product' in fea.qualifiers:
            product=fea.qualifiers['product'][0]
            
        start=fea.location._start.position
        end=fea.location._end.position
        strand=fea.location._strand
        seqref=refseq[0][start:end]
        seqSnp=refseq[0][start:snp[ppos]-1]+snp[palt]+refseq[0][snp[ppos]:end]
        other=refseq[0][start:snp[ppos]-1]+snp[pref]+refseq[0][snp[ppos]:end]
        if other!=seqref:
            print "error no match between seqref and other"
            break
        posAnnot=snp[ppos]-1-start
        if strand==-1:
            seqref=revcom(seqref)
            seqSnp=revcom(seqSnp)
            posAnnot=end-snp[ppos]
        transRef,coRef=translate_dna(seqref,pos=posAnnot)
        transSnp,coSnp=translate_dna(seqSnp,pos=posAnnot)
        if coRef[1]==coSnp[1]:
            syn="Syn"
        else:
            syn="non-syn"
        newtab.append(snp+[gene,db_xref,strand,coRef,coSnp,syn,product])
    else:
        newtab.append(snp+["intergene","NA","NA","NA","NA","NA","NA"])
    


writeCSV(snpsFile[:-4]+"_Annotation.csv",newtab)





























