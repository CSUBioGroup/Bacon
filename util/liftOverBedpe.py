import argparse
import os

#####################################################################################################
#																									#			
#											Functions												#
#																									#			
#####################################################################################################

""" Define a function that splits up bedpe files into 2 temp files """
def splitBedpe(bedpe,tmp1="tmp1.bed",tmp2="tmp2.bed",header="F",verbose="F"):
	
	t1 = open("tmp1.bed",'w')
	t2 = open("tmp2.bed",'w')
	bedpein = open(bedpe,'r')
	
	if header == "T":
		headerline = bedpein.readline()
	
	for l in bedpein:
		e = l.strip().split("\t")
		t1.write("\t".join(e[0:3] + e[6:9])+"\n")
		t2.write("\t".join(e[3:9])+"\n")
		
	t1.close()
	t2.close()
	bedpein.close()
	

""" Define a function that implements liftOver """
def doliftOver(liftOver,chain,infile,verbose="T"):
	
	cmd = " ".join([liftOver,infile,chain,infile + ".success",infile + ".failure"])
	#if verbose == "T":
	print(cmd)
	os.system(cmd)
	

""" Define a function that merges liftOver """
def mergeliftOver(f1,f2,outputfile,verbose="F"):
	
	o = open(outputfile,'w')
	
	# read in file 1 and make dictionary
	readdict = dict()
	f = open(f1,'r')
	for l in f:
		e = l.strip().split('\t')
		readdict[e[3]] = e[0:3]
	f.close()
	
	# read in file2 and print out matches
	f = open(f2,'r')
	n=1
	for l in f:
		e = l.strip().split('\t')
		if (e[3] in readdict) and n<10001:
			r1 = readdict[e[3]]
			#r2 = e
			o.write(r1[0]+"\t"+r1[1]+"\t"+r1[2]+"\t"+e[0]+"\t"+e[1]+"\t"+e[2]+"\t"+e[3]+"\n")
			n=n+1
	f.close()
	o.close()


#####################################################################################################
#																									#			
#										Parse arguments												#
#																									#			
#####################################################################################################


parser = argparse.ArgumentParser(description='wrapper for liftOver to accomodate bedpe files')

# required arguments
parser.add_argument('--lift',        dest='liftOver', 	help='path to liftOver')
parser.add_argument('--chain', 	     dest='chain', 	    help='(e.g. hg19ToHg18.over.chain)')
parser.add_argument('--i', 	         dest='infile', 	    help='input file in bedpe format')
parser.add_argument('--o', 	         dest='outfile', 	    help='output file')
parser.add_argument('--v', 	         dest='verbose', 	    help='verbose' , default = "F")
parser.add_argument('--h', 	         dest='header', 	    help='T /  F if there is a header line', default = "F")

# parse arguments
args = parser.parse_args()

# read in args
LO       = args.liftOver
chain    = args.chain
bedpeIN  = args.infile
bedpeOUT = args.outfile
tmp1     = "tmp1.bed"
tmp2     = "tmp2.bed"
header	 = args.header
verbose	 = args.verbose

#####################################################################################################
#																									#			
#										Run the Code  												#
#																									#			
#####################################################################################################

# break up the files
splitBedpe(bedpeIN,tmp1,tmp2,header,verbose)

# perform liftOver
doliftOver(LO,chain,tmp1,verbose)
doliftOver(LO,chain,tmp2,verbose)

# merge liftOvered files
mergeliftOver(tmp1+".success",tmp2+".success",bedpeOUT,verbose)

# remove tmp files
os.remove(tmp1)
os.remove(tmp2)
os.remove(tmp1+".success")
os.remove(tmp1+".failure")
os.remove(tmp2+".success")
os.remove(tmp2+".failure")





















