
import os
import csv
import sys
import fileinput

args=sys.argv
if len(args)>2:
    fnameI=sys.argv[1]
    thmin=int(sys.argv[2])
    thprop=float(sys.argv[3])
    thqual=int(sys.argv[4])
    
else:
    thqual=150
    thprop=0.2
    thmin=2
    fnameI="/media/Second3TB/Work/WorkVLA/Data/WGS_Analysis/results/20150903_Lorraine/AF2122/Mapped/C12/C12.pileup.NonEqual.vcf"

fnameO=fnameI[:-4]+"_SN.csv"
fnameODUO=fnameI[:-4]+"_DUO.csv"
fnameOIND=fnameI[:-4]+"_INDEL.csv"

fOut=open(fnameO,"w")
writer = csv.writer(fOut, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
fOutDUO=open(fnameODUO,"w")
writerDUO = csv.writer(fOutDUO, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
fOutIND=open(fnameOIND,"w")
writerIND = csv.writer(fOutIND, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
for line in fileinput.input(fnameI):
        if line[0]!="#":
            if "INDEL" not in line and "DP4" in line: 
                line=line.split()
                refGenome=line[0]
                pos=int(line[1])
                ref=line[3]
                alt=line[4]
                qual=float(line[5])
                det=line[7].split(";")
                cov=int([s for s in det if "DP=" in s][0].split("DP=")[1])
                gcov=map(int,[s for s in det if "DP4=" in s][0].split("DP4=")[1].split(","))
                if len(alt)>1:
                    writerDUO.writerow(line)
                else:
                    if ref!="n" and alt!="n" and len(alt)==1 and qual>=thqual and gcov[0]<=thprop*gcov[2] and gcov[1]<=thprop*gcov[3] and (gcov[2]>thmin or gcov[3]>thmin):
                        writer.writerow([refGenome,pos,ref,alt,qual,cov]+gcov)
            else:
                lineS=line.split()
                if "DP4=" in lineS[7]:
                    det=lineS[7].split(";")
                    gcov=map(int,[s for s in det if "DP4=" in s][0].split("DP4=")[1].split(","))
                    if gcov[0]<=thprop*gcov[2] and gcov[1]<=thprop*gcov[3] and (gcov[2]>thmin or gcov[3]>thmin):
                        writerIND.writerow(lineS)
fOut.close()
fOutDUO.close()
fOutIND.close()
