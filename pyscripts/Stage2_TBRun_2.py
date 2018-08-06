#from Bio import SeqIO
import os, sys, os.path
import numpy
import pickle, csv, os, time, pp
#from operator import itemgetter


#python /media/Big/Work/WorkVLA/Projects/RD0070/pipeline/code/Woking/Stage2_TBRun_2.py /media/Big/Work/WorkVLA/Data/TB_Runs/TB_fastqs /media/Big/Work/WorkVLA/Data/TB_Runs/TB_Results /media/Big/Work/WorkVLA/Projects/RD0070/pipeline/stage2 AF2122 20140530_SB4020__WGS-TB-20140521 1 8 false 2
#########

thHitLen=23
thHitSim=85
wl="7"

#   parsing the arguments
args=sys.argv
    
if len(args)<2:
    pathTBfastqs="/media/Second3TB/Work/WorkVLA/Data/TB_Runs/TB_fastqs"
    pathTBRunsResutls="/media/Second3TB/Work/WorkVLA/Data/TB_Runs/TB_Results"
    pathPatterns="/media/Second3TB/Work/WorkVLA/Projects/RD0070/pipeline/stage2"
    TBRun ="20150817_TBFT5168__WGS-TB-20150810"
    ncpus=1
    qth=8
    stage1="false"
    spCountTh=1.5 # times sd of spCountTh for cut off of false space appareances
    spth=1 # if spth is != 0 spCountTh is ignored
else:    
    pathTBfastqs=sys.argv[1]
    pathTBRunsResutls=sys.argv[2]
    pathPatterns=sys.argv[3]
    TBRun=sys.argv[4]
    ncpus=int(sys.argv[5])
    qth=float(sys.argv[6])
    stage1=sys.argv[7]
    spCountTh=float(sys.argv[8])
    spth=float(sys.argv[9])


def listT(matrix):
    return map(list, zip(*matrix))
    
def readTable(fname,ch):
    infile=open(fname,"rb")
    data = csv.reader(infile, delimiter=ch)
    dataOut = [row for row in data]
    infile.close()
    return(dataOut)

def writeCSV(fname,matrix):
    with open(fname, "wb") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print "file "+fname+" saved."
def writeTable(fname,matrix):
    with open(fname, "wb") as fileOut:
        writer = csv.writer(fileOut,delimiter='\t')
        writer.writerows(matrix)
        print "file "+fname+" saved."
def writeListToTextFile(lista,fileName):
    f = open(fileName, 'w')
    f.write("\n".join(lista))
    f.close()
    print "file "+fileName+" saved."

def fastqToFasta(fastqFile):
    fileIn = open(fastqFile, 'rb')
    fileOut = open(fastqFile[:-1]+"a", 'wb')
    line = fileIn.readline()
    while line: 
        if line[0]=="@":
            fileOut.write(">"+line[1:])
            fileOut.write(fileIn.readline())
            fileIn.readline()
            fileIn.readline()
            line = fileIn.readline()
    fileIn.close()
    fileOut.close()

def extractFromFasta(fastaFile,seqsId,outFastaFile):
    f=open(fastaFile,"r")
    seqExtracted=[]
    line=f.readline()
    while line:
        flag=0
        if line.split()[0][1:] in seqsId:
            flag=1
            seqExtracted=seqExtracted+[line[:-1]]
            line=f.readline()
            while line[0]!=">":
                seqExtracted=seqExtracted+[line[:-1]]
                line=f.readline()
        if flag==0:
            line=f.readline()
    writeListToTextFile(seqExtracted,outFastaFile)
    f.close()
    
def fRev(lis):
    if lis[0]<lis[-1]:
        return lis
    else:
        lis.reverse()
        lis[1]=lis[1]+"-RC"
        return lis

def getMotifs(crunchPFile,crunchSFile,fileOutName,thHitLen,thHitSim):    
    crunchP=readTable(crunchPFile,'\t')
    crunchS=readTable(crunchSFile,'\t')
    names=list(set([x[1] for x in crunchP]+[x[1] for x in crunchS]))
    seqPat=[]
    for seqid in names:
        seqsP=[[int(x[8]),x[0]+"-"+str(int(x[7])-int(x[6])+1),int(x[9])] for x in crunchP if x[1]==seqid and int(x[3])>16]
        seqsP=[fRev(x) for x in seqsP]
        seqsS=[[int(x[8]),x[0]+"-"+str(int(x[7])-int(x[6])+1),int(x[9])] for x in crunchS if ((x[1]==seqid and int(x[3])>thHitLen) or (x[1]==seqid and int(x[3])>16 and x[0]=="Spacer02")) and (float(x[2])>thHitSim)]
        seqsS=[fRev(x) for x in seqsS]
        seqs=seqsP+seqsS
        seqPat=seqPat+[[seqid]+sorted(seqs,key=lambda x: x[0])]
    writeCSV(fileOutName+"_motifs.csv",seqPat)
    return spacersInMotifs(seqPat)   

def spacersInMotifs(motifs):
    motifs=[[motif[0],x] for motif in motifs for y in motif[1:] for x in y]
    sps=["01","02","03","04","05","06","07","08","09"]+[str(x) for x in range(10,44)]
    sps=["Spacer"+x for x in sps]
    spacersFound=[]
    for sp in sps:
        spacersSp=[x for x in motifs if sp in str(x[1])]
        if len(spacersSp)>0:
            spacersFound=spacersFound+spacersSp
    return spacersFound        

def countSpacersInMotifs(spacersFound):
    sps=["01","02","03","04","05","06","07","08","09"]+[str(x) for x in range(10,44)]
    sps=["Spacer"+x for x in sps]
    counts=[]
    for sp in sps:
        spCount=[x[0] for x in spacersFound if sp in str(x[1])]
        spCount=len(list(set(spCount)))
        counts=counts+[spCount]
    return counts        

def transfCounts(lis):
    out="r_"
    for li in lis:
        if li<10:
            out=out+str(li)
        if li>=10 and li<20:
            out=out+"a"
        if li>=20 and li<30:
            out=out+"b"
        if li>=30 and li<40:
            out=out+"c"
        if li>=40:
            out=out+"d"
    return out

def checkDifs(par1,par2):
    if par1=="NA" or par1==par2:
        return "NA"
    par1=list(par1)
    par2=list(par2)
    mismatch=[]
    for i in range(len(par1)):
        if par1[i]!=par2[i]:
            if par1[i]==1:
                mismatch=mismatch+[str(i+1)+"p"]
            else:
                mismatch=mismatch+[str(i+1)+"o"] 
    return "-".join(mismatch)
                    
     
def typeOneStrain(strainsDetailsHeader,strainDetails,spoPatterns,pathTBfastqs,pathPatterns,pathAux,wl,pathCrunchs,ppredType,thHitLen,thHitSim,qth,spth,spCountTh,stage1):
    TBRun=strainDetails[strainsDetailsHeader.index('runName')]
    meanCov=float(strainDetails[strainsDetailsHeader.index('meanCov')])
    print strainsDetailsHeader
    print strainDetails
    fileName=strainDetails[strainsDetailsHeader.index('fileName')]
    flag=strainDetails[strainsDetailsHeader.index('flag')]
    if meanCov>=qth:
        if flag not in ["BritishbTB","nonBritishbTB"] or stage1=="false":
            fastaSpacersFile=os.path.join(pathPatterns,"spacers.fasta")
            strainFastqPath=os.path.join(pathTBfastqs)
            strainFilesGZ=[x for x in os.listdir(strainFastqPath) if fileName+"_" in x]
            crunchesSpacers=[]
            for strainFile in strainFilesGZ:
                print "Processing files "+strainFile
                cmd = "gunzip -c "+os.path.join(strainFastqPath,strainFile)+ " > "+ os.path.join(pathAux,strainFile[:-3])
                os.system(cmd)
                outFastaFile=os.path.join(pathAux,strainFile[:-9]+".fasta")
                fastqToFasta(os.path.join(pathAux,strainFile[:-9]+".fastq"))        
                cmd = "formatdb -p F -i " + outFastaFile
                os.system(cmd)
                crunchSpacersFile=os.path.join(pathAux,strainFile[:-9]+"_crunchspacers_"+wl+".txt")
                cmd = "blastall -p blastn -d "+ outFastaFile +" -i "+ fastaSpacersFile + " -o " + crunchSpacersFile + " -q -1 -m 8 -e 5 -W "+wl
                os.system(cmd)
                crunchSpacers=readTable(crunchSpacersFile,'\t')
                if "_R1_" in strainFile:
                    crunchSpacers=[x+["R1"] for x in crunchSpacers]
                else:
                    crunchSpacers=[x+["R2"] for x in crunchSpacers]
                crunchesSpacers=crunchesSpacers+crunchSpacers
                cmd="rm "+os.path.join(pathAux,fileName+"*")
                os.system(cmd)
            
            print "Before Filtering "+str(len(crunchesSpacers))
            writeCSV(os.path.join(pathCrunchs,fileName+"_Crunch.csv"),crunchesSpacers)
            crunchesSpacers=[x for x in crunchesSpacers if ((int(x[3])>thHitLen) or (int(x[3])>16 and x[0]=="Spacer02")) and (float(x[2])>thHitSim)]
            writeCSV(os.path.join(pathCrunchs,fileName+"_CrunchFiltered_"+str(thHitLen)+"_"+str(thHitSim)+".csv"),crunchesSpacers)
            print "After Filtering "+str(len(crunchesSpacers))
            complete = [[x[0][-2:],x[1]] for x in crunchesSpacers]
            complete=map(list,list(set(map(tuple, complete))))
            print "After Removing Same read "+str(len(complete))
            complete = [x[0] for x in complete]
            spCounts=[]
            spotype=['0']*43
            for i in range(len(spotype)):
                if i+1<10:
                    sp = "0"+str(i+1)
                else:
                    sp = str(i+1)
                spCount = complete.count(sp)
                spCounts = spCounts + [spCount]
            if spth==0:
                spCountsMean = numpy.mean([x for x in spCounts if x>0])
                print [x for x in spCounts if x>0]
                print "spCountsMean="+str(spCountsMean)
                spCountSD = numpy.std([x for x in spCounts if x>0])
                print "spCountSD="+str(spCountSD)
                print "spCountTh="+str(spCountTh)
                
                spth = spCountsMean-(spCountTh*spCountSD)
                print "spth="+str(spth)
                if spth < 1:
                    spth = 1
            print "spth="+str(spth)
            for i in range(len(spotype)):
                if spCounts[i] > spth:
                    spotype[i]='1'    
            spotype="".join(spotype)
            match=[x for x in spoPatterns if x[1]==spotype]
            writeListToTextFile(map(str,spCounts),os.path.join(pathCrunchs,fileName+"_Counts.txt"))
            spCounts=transfCounts(spCounts)
            if len(match)==0:
                return [strainDetails+["NA","p_"+spotype,round(spth,2)]+[spCounts]]
            else:
                return [strainDetails+[match[0][0],"p_"+match[0][1],round(spth,2)]+[spCounts]]
        else:
            if stage1=="true":
                print "Stage1_OK"
                return [strainDetails+["Stage1_OK","","",""]]
    else:
        print "LowCoverage"
        return [strainDetails+[""]*(ppredType-len(strainDetails))+["LowCoverage","","",""]]
        
print TBRun
spoPatterns=readTable(os.path.join(pathPatterns,"typePatterns.csv"),",")
strainDetailsFile=os.path.join(pathTBRunsResutls,TBRun,"Stage1",TBRun+"_stage1.csv")
if stage1=="true":
    pathResutls=os.path.join(pathTBRunsResutls,TBRun,"Stage1_Stage2")
else:
    pathResutls=os.path.join(pathTBRunsResutls,TBRun,"Stage2")

if not os.path.exists(pathResutls):
    os.makedirs(pathResutls)
    print "Created directory " + pathResutls

pathCrunchs=os.path.join(pathResutls,"crunches")
pathAux=os.path.join(pathResutls,"extrafiles")
if not os.path.exists(pathCrunchs): os.system("mkdir "+ pathCrunchs)
if not os.path.exists(pathAux): os.system("mkdir "+ pathAux)

strainsDetails=readTable(strainDetailsFile,",")
prunName=strainsDetails[0].index('runName')
pmeanCov=strainsDetails[0].index('meanCov')
pfileName=strainsDetails[0].index('fileName')
pgroup=strainsDetails[0].index('group')
pflag=strainsDetails[0].index('flag')
newstrainsDetails=[strainsDetails[0]+["predType","predictedPattern","th_#copies","copiesFoundperSpacer(a=10-19,b=20-29..)"]]
ppredType=newstrainsDetails[0].index('predType')

print "Stage2. Processing "+str(len(strainsDetails))+" samples"
if stage1=="true":
    outFileName=os.path.join(pathResutls,TBRun+"_stage1_stage2.csv")
else:
    outFileName=os.path.join(pathResutls,TBRun+"_stage2.csv")

if ncpus==1:
    for strainDetails in strainsDetails[1:]:
        print "Processing sample "+str(strainDetails)
        newstrainsDetails=newstrainsDetails+typeOneStrain(strainsDetails[0],strainDetails,spoPatterns,pathTBfastqs,pathPatterns,pathAux,wl,pathCrunchs,ppredType,thHitLen,thHitSim,qth,spth,spCountTh,stage1)
        writeCSV(outFileName,newstrainsDetails)

else:
    ppservers = ()
    job_server = pp.Server(ncpus, ppservers=ppservers)
    print "Starting pp with", job_server.get_ncpus(), "workers"
    jobs = [(strainDetails, job_server.submit(typeOneStrain,(strainsDetails[0],strainDetails,spoPatterns,pathTBfastqs,pathPatterns,pathAux,wl,pathCrunchs,ppredType,thHitLen,thHitSim,qth,spth,spCountTh,stage1),(transfCounts,countSpacersInMotifs,spacersInMotifs,fRev,writeCSV,readTable,writeListToTextFile,getMotifs,extractFromFasta,fastqToFasta,),("csv","os","os.path","numpy",))) for strainDetails in strainsDetails[1:]]
    cont=0
    for strainDetails, job in jobs:
        print "Processing sample "+str(strainDetails)
        newstrainsDetails=newstrainsDetails+job()
        writeCSV(outFileName,newstrainsDetails)
os.system("rm -rf "+ pathAux)