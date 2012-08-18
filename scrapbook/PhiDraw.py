### PhiDraw : a Python-based sequence map and alignment drawing tool for phages ###
### Version 1.0 by Geraldine A. Van der Auwera, AUG 2011 ###
### Feedback welcome at Geraldine_VdAuwera@hms.harvard.edu ###

# Required Python software/libraries :
# BioPython (http://www.biopython.org/)
# ReportLab (http://www.reportlab.org/)

print "\n*** Running PhiDraw ***"

#### CONFIGURATION SECTION ####

# This config section contains the filename specifications for the input files, as well as the basic definition of all graphics-related variables used to draw plasmid maps and alignments using PlasmiG. It is located at the start of the main script to make the definitions available throughout the program and to make it easy to find for customization. The values can be modified within this section (recommended), or they can be overridden by defining alternative values within the main script (not recommended). For ease of use, it is convenient to save specific customized versions of this file as separate copies, preferably under an explicit name. In any case, the default values indicated between brackets should be left for future reference. 

print "Loading configuration section..."

### Load all modules and dependencies first ###

## Import general modules ##
import re
import sys
from sys import argv
from math import log
from Bio import GenBank, SeqFeature, Seq, SeqIO
from Bio.Seq import Seq
from reportlab.pdfgen import *
from reportlab.pdfgen import canvas
from reportlab.lib import *
from reportlab.lib.units import *
from reportlab.lib.colors import *

print "    all modules imported"

### Input filenames and paths ###

## Input files ##
aln_seqlist 	= argv[1]			# list of plasmids
aln_outfile 	= argv[2]			# pdf output file

### Base units, measurements and colors ###

## General proportions ##
u 		= 0.025			# conversion factor for sequence length-related values (0.0125, increase for short sequences)
hmar 	= 1*cm			# horizontal margin to canvas (1*cm)
vmar 	= 2*cm			# vertical margin to canvas (2*cm)
pNsize 	= 4*cm			# width to set aside for plasmid names (2.5*cm)
di 		= 0.18*cm		# half-length of interruption/frameshift tick marks (0.18*cm)
doL 	= 0.7*cm		# distance of ORF labels from their ORF (0.7*cm)
minL 	= 400			# minimum ORF size (in base pairs) for 'full arrows' (700, decrease for short seq)
w 		= 0.3*cm     	# half-width of the tail (0.3*cm)
h 		= 0.175*cm		# distance from the side tips of the head to the neck (0.175*cm)

## Typeface, fonts, sizes ##	
rFont 	= "Helvetica"			# regular type ("Helvetica")
bFont 	= "Helvetica-Bold"		# bold type ("Helvetica-Bold")
LfSize 	= 14					# large font size (14) 
NfSize 	= 12					# normal font size (12)
SfSize 	= 10					# small font size (10)

## special feature settings ##
SFX 	= 'off'			# allow drawing special sequence features beside CDS/ORF ('on' or 'off')
osym 	= '*'			# symbol for origin of sequence (*)
snp		= '*'			# symbol for SNP (*)

## Alignment settings ##
clump 	= 'off'			# OPTION to clump or not to clump, that is the question ('on' or 'off') #DEV# disabled for now
Llimit 	= 200			# threshold distance (in bases) for clumping segments together (100)
multipass = 1			# set how many clumping iterations should be passed (3)
mLS 	= 30			# minimum sequence length of segments to qualify for drawing (3)
heat 	= 'on'			# OPTION to color shaded areas according to similarity ('on' or 'off')
Dlimit 	= 0.2 			# maximum ratio of asymmetry within segment pairs for shading (0.2)
chopper	= 'on'			# OPTION to chop long segments into smaller ones
chopX	= 1000			# size of chopped up segments
idpt 	= [60,50,40,30]	# array of 4 threshold values for heatmap shading of alignment + 1 minimum cutoff (95,85,70,50)
dBL	 	= 3*cm			# distance between plasmid baselines (4*cm)
da 		= 0.4*cm		# distance of alignment tick marks from the corresponding plasmid baseline (0.4*cm)
idpdiff	= 1 			# max difference in segment pair identity for glomping
seg_space = 100			# max distance between consecutive segments for glomping
tm 		= 0.2*cm		# size of alignment tick marks (0.2*cm)

## Color definitions for alignment shading cues ##
shtop 	= HexColor('#444444')		# top similarity class (HexColor('#444444'))
shupmid = HexColor('#777777')		# upper middle class (HexColor('#777777'))
shlomid = HexColor('#BBBBBB')		# lower middle class (HexColor('#BBBBBB'))
shlow 	= HexColor('#DDDDDD')		# low similarity class (HexColor('#DDDDDD'))
shblank = HexColor('#FFFFFF')		# lower than cutoff (HexColor('#FFFFFF'))
shdef	= HexColor('#DDDDDD')		#  default shading color (HexColor('#DDDDDD'))
idpc 	= (shtop,shupmid,shlomid,shlow)
# examples of other suitable gradients using Reportlab colors : 
# dimgray / lightslategray / lightsteelblue / lavender		('atmosphere' theme)
# darkkhaki / khaki / palegoldenrod / lightgoldenrodyellow	('desertstorm' theme)

## Color definitions for functional categories of CDS ##
# functions specified by tags in the GenBank files (more can be added ad libitum)
# color names as specified in the ReportLab library 
# general conjugation/transfer terms
transfer 	= mediumslateblue	# (mediumslateblue)
T4SS 		= mediumslateblue	# (mediumslateblue)
T4SLS 		= mediumorchid		# (mediumorchid)
Tra 		= steelblue			# (steelblue)
Trb 		= mediumslateblue	# (mediumslateblue)
#  control functions
replication = darkorange	# main replication protein (darkorange)
inheritance = darkorange	# partitioning, addiction etc. (gold)
regulation 	= yellow		# transcriptional regulators (yellow)
# miscellaneous
patho 		= red 			# pathogenicity determinants (red)
mobile 		= lime			# mobility, transposases etc. (lime)
other 		= deepskyblue	# other functions (deepskyblue)
accessory 	= deepskyblue	# accessory genes (deepskyblue)
# default for any and all non-specified stuff
default 	= seashell		# (seashell)

## Automated labeling derived from GB locus tags ##		##DEV## needs more work
# uses pattern recgnition to label features of interest automatically
# defaults are the most general Gram-negative plasmid genes
# (Tra|Trb|Mob|RepA|TrfA|ParA)
#ORF_Lloci = '(Tra|Trb|Mob|RepA|TrfA|ParA)'


### Legend elements ###

## Scale of sequence, in legend ##
scX 	= 0*u			# starting point on the X axis(0)
incrT 	= 1000*u		# increment size in basepairs (1000, signifying 1 kb)
incrN 	= 5				# increment number (5, or 1 for short sequences)
dip 	= -2*cm			# starting point on Y axis (-2*cm) 
dop 	= 0.3*cm		# size of vertical tick

## Color key settings ##
ckX 		= 5*cm		# starting point on X axis (5*cm)
ckY 		= dip		# starting point on Y axis (dip, same as legend scale, but can take any numerical value)
ck_hsp 		= 170		# horizontal spacing between columns (170)
ck_vsp 		= 17		# vertical spacing between items (17)
ck_boxX 	= 12		# horizontal side size of box (10)
ck_boxY 	= 12		# vertical side size of box (10)
ck_htxoff 	= 20		# horizontal offset for text from box origin (1) 
ck_vtxoff 	= 1			# vertical offset for text from box origin (20)
lay_MAX 	= 2			# maximum items in a column (4, can also be a fraction of len(functions)) 

## Heatmap key settings ##
hk_u		= 1.5		# length conversion factor (1.5)
hkY			= dip		# starting point on Y axis (dip)
hk_boxY		= ck_boxY	# horizontal side size of box (ck_boxY)
hk_htxoff 	= -100		# horizontal offset for text from box origin (1) 
hk_vtxoff 	= 1			# vertical offset for text from box origin (20)
aligner		= 'ClustalW nt'	# text string for legend indicating basis for heatmap shading

## Array of names for function types ##
# specifies the label text for the legend color key (more can be added ad libitum)
dfuncts = [(replication, "Replication"),		# ("Replication")
		#(regulation, "Regulation"),				# ("Regulation")
		#(inheritance, "Inheritance"),			# ("Inheritance", could also be "Stability")
		#(Tra, "Tra genes"),					# ("Tra genes")
		#(Trb, "Trb genes"),					# ("Trb genes")
		#(T4SS, "T4SS"),							# ("T4SS", could also be "Type IV Secretion")
		#(T4SLS, "T4SLS"),						# ("T4SLS", could also be "Type IV-Like Secretion")
		#(transfer, "Transfer"),				# ("Transfer")
		#(patho, "Pathogenicity"),				# ("Pathogenicity", could also be "Virulence")
		#(mobile, "Transposition/MGE"),			# ("MGE", could also be "Mobile Genetic Elements")
		#(other, "Other functions"),			# ("Other functions")
		(accessory, "Accessory genes")]			# ("Accessory genes")
		#(duo, "Unknown function"),				# ("Unknown function")
		#(trio, "Unknown function"),			# ("Unknown function")
		#(none, "No known homologs")			# ("No known homologs")

print "    OK"

#### END OF CONFIGURATION SECTION ####


#### LIBRARY OF FUNCTIONS ####

print "Loading library of functions... "

### Data input processing functions ###

## EnsureDir : Creates directory if it doesn't exist ##
def EnsureDir(directory) :
    dirhandle 	= os.path.dirname(directory)
    if not os.path.exists(dirhandle):
        os.makedirs(dirhandle)
        
## Strippa : Cleans up raw data from file read-in ##
def Strippa(line) :
	stripline 	= line.strip("\n")
	stripline 	= stripline.strip("\r")
	line_array 	= stripline.split("\t")
	return line_array

## Cleppie : Cleans up single lines and arrays plasmids ##
def Cleppie(line_array) :
	plasmids 	= []
	pnames		= []
	for rawline in line_array :	
		line 	= Strippa(rawline)	
		plasmids.append(line)
		pnames.append(line[0])
	pN = len(plasmids)		# length of the array is the number of plasmids
	return pN,pnames,plasmids

## Off7 : Calculates coordinate changes according to plasmid-specific offset ##
def Off7(coord,pLen,offset) :
	# calculate linear coordinate change
	coft = (abs(coord)-abs(offset)+1)
	if coft < 0 :
		coff 	= pLen + coft
	else :
		coff 	= coft
	# invert sequence if offset is negative
	if offset < 0 :
		cofi 	= pLen - coff
		sign7 	= -1
	else :
		cofi 	= coff
		sign7	= 1
	# recall original sign 
	if coord < 0 :
		sign0 	= -1
	else :
		sign0 	= 1
	# set final sign
	coff7 	= sign0 * cofi
	return coff7,sign7

## Quixote : Does a quick'n'dirty identification of file format based on file extension ##
def Quixote(filename) :
	gb_ext 	= re.compile(r"\.gb")
	fa_ext 	= re.compile(r"\.fa")
	sq_ext 	= re.compile(r"\.seq")
	genbank = gb_ext.search(filename)
	fasta 	= fa_ext.search(filename)
	seq 	= sq_ext.search(filename)
	if genbank :
		format = 'genbank'
	elif fasta : 
		format = 'fasta'
	elif seq : 
		format = 'seq'
	else :
		format = 'unidentified'
	return format
	
## pLonk : Parses sequence files to get their length ##
def pLonk(plasmids) :
	pLenks = []
	pLasmids = []
	for (pName,seq_infile,offset,order) in plasmids :
		fhandle = open(seq_infile, 'r')	# load plasmid sequence file
		# evaluate file name to detect format using Quixote [ filename ]
		format =Quixote(seq_infile)
		if format == 'genbank' :
			parser = GenBank.FeatureParser()
			gb_entry = parser.parse(fhandle)	
			pLen = len(gb_entry.seq)		# read in length of plasmid sequence
			print pName, pLen
		elif format == 'fasta' or format == 'seq' :
			for fa_entry in SeqIO.parse(fhandle, "fasta") :
				pLen = len(fa_entry.seq)		# read in length of plasmid sequence
		else :		
			print "TERMINAL ERROR : file format not recognized for " + pName + " !!!"
			break;
		fhandle.close()					# close sequence file (to free up memory)
		pLenks.append(pLen)
		pLasmids.append((pName,seq_infile,int(pLen),int(offset)))
		pLenks.sort()					
		pLenks.reverse()
		pLen_MAX = pLenks[0]
	return pLenks,pLen_MAX,pLasmids		

### Plasmid alignment map ###

## Flipper : Flip negative coordinates ##
def Flipper(axr,bxr) :
	if axr < 0 :
		ax,bx = abs(bxr),abs(axr)
	else :
		ax,bx = axr,bxr
	return ax,bx 

## RCshift : Performs reverse-complementation on non-Seq ojects ##
def RCshift(seqslice) :
	tempseq	= Seq(seqslice)
	rc_seq	= str(tempseq.reverse_complement())
	return rc_seq

### Sequence features : processing and  drawing tools ###

## Canvasser : Initialize canvas ##
def Canvasser(hCan,vCan,transX,transY,outfile) :
	canvasN = canvas.Canvas(outfile, pagesize=(hCan,vCan)) 
	canvasN.translate(transX,transY)	  
	canvasN.setStrokeColor(black)	 
	canvasN.setFillColor(white) 
	canvasN.setLineWidth(1) 
	canvasN.setLineJoin(1) 
	canvasN.setLineCap(0) 
	return canvasN
	
## BaseDraw : draws plasmids baselines and features ##
def BaseDraw(plasmids) :
	ordN = 0
	for (pName,seq_infile,pLen,offset) in plasmids :
		# set Y axis once and for all for the plasmid being processed
		y0 = (pNs-ordN)*dBL				# starts from the top
		pLeni = int(pLen)
		print 'offset', offset
		offset = int(offset)
		# draw plasmid baseline
		BaseL(ordN,pName,pLeni,y0,canvas_main)
		# label the baseline with plasmid name and size
		LabeL(ordN,pName,pLeni,y0,canvas_main)
		# evaluate file name to detect format using Quixote [filename]
		format = Quixote(seq_infile)
		# mark up sequence origin if there is an offset
		if offset < -1 or offset > 1 :
			Zs,dir 		= Off7(1,pLeni,offset)
			xs 			= Zs*u
			canvas_main.setFont(bFont,NfSize) 
			canvas_main.drawString(xs,y0+da/2,osym) 
		# filter and draw annotation features 
		if format == 'genbank' :
			# load GB file to filter features
			parser = GenBank.FeatureParser()
			fhandle = open(seq_infile, 'r')			# load GenBank file
			gb_entry = parser.parse(fhandle)
			ORFcnt = 0
			for feature in gb_entry.features :
				if feature.type == 'CDS' or feature.type == 'cds':		# draw CDS using ORFeus
					ORFcnt += 1
					ORFeus(feature,pLeni,offset,y0,ORFcnt)
				elif SFX == 'on' :
					if feature.type == 'SNP' :	
						Snippit(feature,pLeni,offset,y0)		# draw asterisk at feature location
					if feature.type == 'IR' :
						IRFlag(feature,pLeni,offset,y0)			# draw flag at feature location
				# need other functions for other features ( with conditional, default switch off)
			fhandle.close()
			print "    got a GenBank-style file for " + pName + " with " + str(ORFcnt) + " ORFs"
		else : 					
			# no features so just skip this step
			print "    got a non-genbank-style file for " + pName + "; no features to draw"
		# increment plasmid ordinal count
		ordN = ordN + 1
		print "    " + pName + " (" + str(pLeni) + " bp) drawn with " + str(ORFcnt) + " ORFs"
	print "    OK"
					
## BaseL : Draws plasmid baselines ##
def BaseL(ordN,pName,pLen,y0,canvas_def) :
	canvas_def.setLineWidth(3) 
	x0 = 0						# all sequences are aligned on the left
	x1 = pLen*u					# convert sequence length to workspace dimensions
	pBL = canvas_def.beginPath()			
	pBL.moveTo(x0,y0)
	pBL.lineTo(x1,y0)
	canvas_def.drawPath(pBL, stroke=1, fill=0) 
	pBL.close()

## LabeL : Labels baselines with plasmid name and size ##
def LabeL(ordN,pName,pLen,y0,canvas_def) :
	canvas_def.setFillColor(black)
	x0 = -pNsize				# retreat into the left margin to write out name and size
	y1 = y0 + ck_vsp/10			# bump name up a bit from BL level
	y2 = y0 - ck_vsp			# bump size down a bit from BL level
	pLenStr = str(float(pLen)/1000)		# converts number to suitable form
	canvas_def.setFont(bFont,LfSize) 
	canvas_def.drawString(x0,y1,pName) 			
	canvas_def.setFont(rFont,NfSize) 
	canvas_def.drawString(x0,y2,pLenStr+' kb') 	

## ORFeus : Draws plasmid ORFs ##
def ORFeus(feature,pLen,offset,y0,count) :
	canvas_main.setLineWidth(1) 
	# extract relevant properties of the CDS feature
	locus 	= str(feature.qualifiers.get('locus_tag'))
	fct 	= feature.qualifiers.get('fct_class')
	hits 	= feature.qualifiers.get('db_hits')
	shape 	= feature.qualifiers.get('shape')
	# evaluate what strand the ORF is on (determines direction of arrow)
	cstrand = feature.strand 
	if cstrand == None:
		cstrand = 1
	# take start and end points
	location 	= feature.location		
	Zs 			= location.nofuzzy_start
	Ze 			= location.nofuzzy_end
	featL 		= Ze - Zs
	# run offsetting function
	mode = 'feature'
	Zs7,dir = Off7(Zs,pLen,offset)
	Ze7,dir = Off7(Ze,pLen,offset)
	if dir == -1 :
		cstrand = -cstrand
		Zs,Ze 	= Ze7,Zs7
	else :
		Zs,Ze 	= Zs7,Ze7
	# calculate X axis coordinates (expr of cstrand has changed!)
	if cstrand == -1 :	# reverse orientation
		xs,xe 	= Ze*u,Zs*u		# start and end
		xn 		= xe+minL*u		# neck of arrow
	else :				# forward orientation
		xs,xe 	= Zs*u,Ze*u		# start and end
		xn 		= xe-minL*u		# neck of arrow
	midLZ 	= ((Zs+Ze)/2)*u	# middle of ORF for optional label
	# set Y axis coordinates
	yt,yb,ynt,ynb = y0+w,y0-w,y0+h,y0-h
	# assign colors to function tags
	if fct == ['transfer'] : 
		canvas_main.setFillColor(transfer) 
	elif fct == ['T4SS'] :	
		canvas_main.setFillColor(T4SS) 
	elif fct == ['T4SLS'] :	
		canvas_main.setFillColor(T4SLS) 
	elif fct == ['tra'] :	
		canvas_main.setFillColor(transfer) 
	elif fct == ['trb'] :	
		canvas_main.setFillColor(Trb)
	elif fct == ['replication'] : 
		canvas_main.setFillColor(replication) 
	elif fct == ['inheritance'] : 
		canvas_main.setFillColor(inheritance) 
	elif fct == ['regulation'] : 
		canvas_main.setFillColor(regulation) 
	elif fct == ['patho'] : 
		canvas_main.setFillColor(patho) 
	elif fct == ['AT'] : 
		canvas_main.setFillColor(patho)
	elif fct == ['mobile'] : 
		canvas_main.setFillColor(mobile)
	elif fct == ['ice'] : 
		canvas_main.setFillColor(mobile)
	elif fct == ['other'] : 
		canvas_main.setFillColor(other)
	elif fct == ['accessory'] : 
		canvas_main.setFillColor(accessory)
	else : 
		canvas_main.setFillColor(default) 
	# initialize path
	pORF = canvas_main.beginPath() 
	# draw square-shaped ORFS
	if shape == ['square'] : 	
		pORF.moveTo(xs,ynt)
		pORF.lineTo(xe,ynt)
		pORF.lineTo(xe,ynb)
		pORF.lineTo(xs,ynb)		
		pORF.lineTo(xs,ynt)	
	# draw triangle-shaped ORFS		
	elif featL <= minL :	
		pORF.moveTo(xs,yt)
		pORF.lineTo(xe,y0)
		pORF.lineTo(xs,yb)
		pORF.lineTo(xs,yt)	
	# draw arrow-shaped ORFS
	elif featL > minL :
		pORF.moveTo(xs,ynt)
		pORF.lineTo(xn,ynt)
		pORF.lineTo(xn,yt)
		pORF.lineTo(xe,y0)
		pORF.lineTo(xn,yb)
		pORF.lineTo(xn,ynb)
		pORF.lineTo(xs,ynb)
		pORF.lineTo(xs,ynt)
	# finalize object path
	canvas_main.drawPath(pORF, stroke=1, fill=1) 
	pORF.close()	
	canvas_main.setFillColor(black)
	canvas_main.setFont(rFont, SfSize) 
	canvas_main.drawCentredString(midLZ,y0-doL,str(count))
	canvas_main.setFont(rFont, NfSize) 

## Snippit : Draws an asterisk at the location of a SNP ##
def Snippit(feature,pLen,offset,y0) :
	location 	= feature.location		
	Zs 			= location.nofuzzy_start
	# run offsetting function
	mode 		= 'feature'
	Zs,dir 		= Off7(Zs,pLen,offset)
	# draw asterisk
	xs 			= Zs*u
	canvas_main.setFont(bFont,NfSize) 
	canvas_main.drawString(xs,y0+da/2,snp) 
	
## IRflag : Draws a flag symbol (default for IR of transposons and/or insertion sequences) ##
def IRFlag(feature,pLen,offset,y0) :
	location	= feature.location
	Zs 			= location.nofuzzy_start
	Ze 			= location.nofuzzy_end
	vdir 		= feature.qualifiers.get('vdir')
	# evaluate what strand the ORF is on (determines direction of flag)
	cstrand = feature.strand 
	# run offsetting function
	mode = 'feature'
	Zs7,dir = Off7(Zs,pLen,offset)
	Ze7,dir = Off7(Ze,pLen,offset)
	if dir == -1 :
		Zs,Ze 	= Ze7,Zs7
		cstrand	= -cstrand
	else :
		Zs,Ze 	= Zs7,Ze7
	# convert sequence coordinates to canvas coordinates
	Cs,Ce = Zs*u,Ze*u
	# size and direction settings for flag symbol
	if vdir == 'up' :		#to make the flag point upwards
		Ch 		= y0+15
		Chtip 	= Ch+16
		Chmip 	= Ch+12
		Chmid 	= Ch+8
	elif vdir == 'dn' :		#to make the flag point downwards
		Ch 		= y0-15
		Chtip 	= Ch-16
		Chmip 	= Ch-12
		Chmid 	= Ch-8
	if cstrand == -1 :		#to make flag point left (for right IR)
		Chick 	= -4
	elif cstrand == 1 :		#to make flag point right (for left IR)
		Chick 	= 4
	# draw flag symbol
	ir = canvas_main.beginPath()
	canvas_main.setLineWidth(1)
	canvas_main.setFillColor(mobile)
	ir.moveTo(Ce,Ch)
	ir.lineTo(Ce,Chtip)
	ir.lineTo(Cs+Chick,Chmip)
	ir.lineTo(Cs,Chmid)
	canvas_main.drawPath(ir, stroke=1, fill=1)	#close symbol path


### Legend elements ###

## Draws the sequence scale bar ##
def SeqScale(scX,incrT,incrN,dip,dop) :
	canvas_main.setLineWidth(1.2)
	canvas_main.setFillColor(black)
	incrCNT = 0							# initialize count of increments
	psc = canvas_main.beginPath()
	psc.moveTo(scX,dip-dop)				# start at beginning (duh!)
	psc.lineTo(scX+incrT*incrN,dip-dop)	# draw the scale baseline
	# draw ticks until the max number of increments is reached
	while incrCNT <= incrN :
		psc.moveTo(scX+incrT*incrCNT,dip-dop) 
		psc.lineTo(scX+incrT*incrCNT,dip)	
		incrCNT = incrCNT + 1			# annoyingly, there's no increment function in python
	canvas_main.drawPath(psc, stroke=1, fill=0) 
	psc.close()
	# write out scale extremities values (needs hand-fix if not using kbs)
	canvas_main.setFont(rFont,NfSize) 
	canvas_main.drawRightString(scX,dip+dop,'0')
	canvas_main.drawString(scX+incrT*incrN,dip+dop,str(incrN)+' kb')

## CKeyOne : Draws a color key for functions, with corresponding text (consider doing symbol key) ##
def CKeyOne(functions,lay_MAX) :
	canvas_main.setLineWidth(1)
	canvas_main.setLineCap(0)
	canvas_main.setFillColor(black)
	canvas_main.setFont(rFont,LfSize) 
	# use counters to organize items into columns of limited length
	stk_CNT = 0							# initialize counter of stacks
	lay_CNT = 0							# initialize counter of layers
	# requires items (function color + corresponding text) to be loaded in specific order in the main script!
	for (type,name) in functions :
		if lay_CNT == lay_MAX :
			stk_CNT = stk_CNT + 1
			lay_CNT = 0
		canvas_main.setFillColor(type)
		# draw the colored rectangle in the right column (stack) and at the right level (layer)
		rectX,rectY = ckX+ck_hsp*stk_CNT,ckY-ck_vsp*lay_CNT
		canvas_main.rect(rectX,rectY,ck_boxX,ck_boxY,fill=1)
		# add the corresponding text next to the colored box
		canvas_main.setFillColor(black)
		txtX,txtY = ckX+ck_hsp*stk_CNT+ck_htxoff,ckY-ck_vsp*lay_CNT+ck_vtxoff
		canvas_main.drawString(txtX,txtY,name)
		# finish processing each item by incrementing the layer counter to go to the next 'line'
		lay_CNT = lay_CNT + 1

	
print "    OK"
			
#### END OF LIBRARY OF FUNCTIONS ####


#### LINEAR ALIGNMENT SCRIPT ####

# This is the main script for drawing linear alignments of plasmid sequences. Explanations regarding the required input files, their formatting constraints and how the program uses the data are provided at the end of this document. 

### Actual program ###

## Data processing ##
# Load in basic plasmid info from the plasmid list
print "Processing input file \"" + aln_seqlist + "\" ..."
inlist 				= open(aln_seqlist, 'r') 
headlineP 			= inlist.readline()		# set aside the first line 
rawlinesP 			= inlist.readlines()	# read in the rest of the lines
inlist.close()
pN,pnames,plasmids 	= Cleppie(rawlinesP)	# prepare array of plasmids and count them
pNs = pN-1			# adjusted pN for evaluating array indices starting from 0
print "    found " + str(pN) + " plasmids to align :"
print pnames

# Get plasmid lengths from their sequence files (GenBank or Fasta)
pLenks,pLen_MAX,pLasmids = pLonk(plasmids) 	# returns plasmid lengths, max and augmented array pLasmids

print "    OK"

## Graphical output ##

# Calculate main canvas dimensions
print "Calculating main canvas dimensions..."
hCan 	= pLen_MAX*u + hmar*2 + pNsize			# set horizontal size
vCan 	= pNs*dBL + lay_MAX*ck_vsp + vmar*2		# set vertical size
transX 	= hmar + pNsize							# set zero of X axis 
transY 	= vmar + lay_MAX*ck_vsp					# set zero of Y axis
print "    dimensions are " + str(hCan/cm) + " cm wide and " + str(vCan/cm) + " cm high"
print "    origin translated by " + str(transX/cm) + " right and " + str(transY/cm) + " up"

# Set up main canvas
canvas_main 	= Canvasser(hCan,vCan,transX,transY,aln_outfile) 

print "    OK"
	
# Process features of plasmids in array pLasmids :
print "Drawing plasmids baselines and features..."

# Draw plasmid baselines and features
BaseDraw(pLasmids)

# Draw legend (scale and color key)
print "Drawing legend (scale and color keys)..."
SeqScale(scX,incrT,incrN,dip,dop)	# draw scale
CKeyOne(dfuncts,lay_MAX)			# draw color key for functions 
print "    OK"

# Write to file and finalize the figure
canvas_main.showPage()	 
canvas_main.save()

print "The program run is finished. Figure name is \"" + aln_outfile + "\".\n*** END ***\n"

### END OF THE LINEAR ALIGNMENT SCRIPT ###

