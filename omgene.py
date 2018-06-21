#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2017 Michael Dunne
#
# OMGene version 0.1.0
#
# This program (OMGene) is distributed under the terms of the GNU General Public License v3
#
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# For any enquiries send an email to Michael Dunne
# michael.dunne@worc.ox.ac.uk

__author__ = "Michael Dunne"
__credits__ = "Michael Dunne, Steve Kelly"

########################################################
################ Safely import packages ################
########################################################

errors = []
libnames = ["csv", "re", "os", "sys", "itertools", "copy", "subprocess", "multiprocessing", \
	"commands", "datetime", "tempfile", "pickle", "string", "pyfaidx", "argparse", \
	"scipy", "math", "warnings", "random"]

for libname in libnames:
    try:
        lib = __import__(libname)
    except ImportError as e:
	errors.append(e)
    else:
        globals()[libname] = lib

try:
	import networkx as nx
except ImportError as e:
	errors.append(e)
try:
	from shutil import copyfile
except ImportError as e:
        errors.append(e)
try:
	from collections import Counter
except ImportError as e:
        errors.append(e)
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.SubsMat import MatrixInfo as matlist
except ImportError as e:
        errors.append(e)
try:
	from scipy.special import binom as bn
except ImportError as e:
        errors.append(e)
try:
	import numpy as np
except ImportError as e:
	errors.append(e)

if errors:
        print("Missing modules :(\nThe following module errors need to be resolved before running OMGene:")
        for error in errors: print("-- " + str(error))
        sys.exit()

##########################################
# Check programs
##########################################

def canRunCommand(command, qAllowStderr = False):
#        sys.stdout.write("Test can run \"%s\"" % command)
        capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = [x for x in capture.stdout]
        stderr = [x for x in capture.stderr]
        return check(len(stdout) > 0 and (qAllowStderr or len(stderr) == 0))

def canRunSpecific(line, lineFormatted):
        return True if canRunCommand(line) else outputBool(False, "ERROR: Cannot run " + lineFormatted, \
                "    Please check "+ lineFormatted +" is installed and that the executables are in the system path\n")

def canRunMinusH(package, packageFormatted):
        return canRunSpecific(package + " -h", packageFormatted)

def canRunAwk():
        return returnsExpected("awk '{print 4 + 5}'", "awk '{print 4 + 5}'", "a", "9\n")

def canRunMan(package, packageFormatted):
        return canRunSpecific("man " + package, packageFormatted)

def check(boolVal, trueMsg=" - ok", falseMsg=" - failed", trueExtra="", falseExtra=""):
        return outputBool(True, trueMsg, trueExtra) if boolVal else outputBool(False, falseMsg, falseExtra)

def returnsExpected(message, command, contentsIn, expectedOut):
       # sys.stdout.write("Test can run \""+message+"\"\n")
        path_tf1 = tempfile.mktemp(); path_tf2 = tempfile.mktemp()
        write(contentsIn, path_tf1)
        callFunction(command + " " + path_tf1 + " > " + path_tf2)
        res = read(path_tf2)
        return res == expectedOut

def outputBool(boolVal, msg, msgExtra):
#        print msg
        if msgExtra: print msgExtra
        return boolVal


def checkShell():
        """ Run quick checks to ensure the relevant packages are installed, and that they do
                the things we need them to do (applies in particular to bedtools and to R).
                - awk           - exonerate   - bedtools
                - echo          - sed
                - cut           - grep
                - sort          -tac/cat
                """
	sprint("Checking installed programs...")
        checks = []
        for program in ["sed", "echo", "cut", "sort", "grep", "uniq", "tac", "cat", "mktemp"]:
                checks.append(canRunMan(program, program))
        for program in ["exonerate", "bedtools"]:
                checks.append(canRunMinusH(program, program))
        #Check presence of orthofinder
        checks += [canRunAwk()]
        if not all(checks):
                print("\nSome programs required to run omgene are not installed or are not callable.\nPlease ensure all of the above programs are installed and in the system path.")
                sys.exit()

##########################################
# Define blosum matrices
##########################################

static_aa_basic = ['A', 'C', 'B', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y', 'X', 'Z']
static_aa = static_aa_basic + ['-']

def goodBlosum(singleGap=-1, doubleGap=0.5):
	blos = copy.deepcopy(matlist.blosum62)
	aa = static_aa_basic
	for a in aa:
		for b in aa:
			if (a,b) in blos:
				blos[(b,a)] = blos[(a,b)]
	for a in aa:
		blos[(a, "-")] = singleGap
		blos[("-", a)] = singleGap
	blos["-","-"] = doubleGap
	return blos

def blosMat(blosdict):
	blosmat = np.zeros([len(static_aa), len(static_aa)])
	for i,e in enumerate(static_aa):
		for j,f in enumerate(static_aa):
			blosmat[i,j] = blosdict[(e,f)]
	return blosmat

static_blosdict = goodBlosum(-1, 0)
static_blosmat = blosMat(static_blosdict)
static_blospos = dict((e,i) for i,e in enumerate(static_aa))

# Save time later by storing calculated binomial values
static_binoms = {}


########################################
# Generic Utilities
########################################

def lprint(l):
	"""Print a list piece by piece
	"""
	for i in l: print i

def dprint(d):
	"""Print a dict piece by piece
	"""
	for k in d: print k, d[k]

def printSeqs(seqs):	
	for s in seqs: 	sprint(str(s.seq))

def sprint(string):
	print string

def read(path_file, tag="r"):
        with open(path_file, tag) as f:
                return f.read()

def readCsv(path_file, ignoreBlank = True, ignoreHash = False):
	with open(path_file, "r") as f:
		data = [line for line in csv.reader(f, delimiter="\t") if \
				(not ignoreBlank or ''.join(line).strip()) and (not ignoreHash or not line[0].startswith('#'))]
	return data

def writeCsv(rows, path_file):
	with open(path_file, "w") as f:
		writer = csv.writer(f, delimiter="\t", quoting = csv.QUOTE_NONE, quotechar='')
		writer.writerows(rows)

def write(msg, path_dest):
	with open(path_dest, "w") as f:
		f.write(str(msg))

def makeIfAbsent(path_dir):
	if not os.path.exists(path_dir): os.makedirs(path_dir)
	return path_dir

def concatFiles(fileslist, dest):
	with open(dest, "w") as o:
		for file in fileslist:
			with(open(file, "r")) as b: o.write(b.read())
	
def callFunction(str_function):
	"""Call a function in the shell
	"""
	subprocess.call([str_function], shell = True)

def grabLines(command):
	return commands.getstatusoutput(command)[1].split("\n")

def interpolate(l1, l2=[]):
	return range(min(l1+l2), max(l1+l2))

def deleteIfPresent(path):
	"""Delete a folder which may or may not exist
	"""
	try:
		os.remove(path)
	except OSError:
		pass

def async(pool, function, args):
	"""Run asynchronously
	"""
	pool.apply_async(function, args=args)
	
def pause():
	callFunction("sleep 10000")

def checkFileExists(path_file):
	return path_file if os.path.exists(path_file) else sys.exit("File does not exist: " + path_file)

def idict(keys, defaultvalue):
	return dict((k, copy.deepcopy(defaultvalue)) for k in keys)

def delIndices(l, indices):
	# Only delete indices that exist, and do them in the right order!
	k = l[:]
	indices_sorted = sorted(set(i % len(l) for i in indices if -len(l) <= i < len(l)), key = lambda x: -x)
	for i in indices_sorted: del k[i]
	return k

def flattenLol(lol):
	"""Flatten a list of lists
	"""
	res = []
	for l in lol: res += l
	return res

def deDup(dlist):
        nlist = []
        for i in dlist:
                if i not in nlist:
                        nlist.append(i)
        return nlist

def mergeDicts(dicts):
	# Note - identical keys are flattened.
	res = {}
	for d in dicts:
		for k in d:
			res[k] = d[k]
	return res

def wrapString(string, n):
	return [string[start:start+n] for start in xrange(0, len(string), n)]

def blankFile(path_file):
	dirname = os.path.dirname(path_file)
	makeIfAbsent(dirname)
	callFunction("echo -n \"\" >" + path_file)
	return path_file

def anyperm(listin, checker):
	return any([(list(i) in checker) for i in itertools.permutations(listin)])

########################################
# Feature validation
########################################

# These get set on entry.
static_do = ["gc", "gt"]
static_ac = ["ag"]
static_sc = ["atg"]

def isStartCodon(c):
	return isFeature(c, 3, static_sc)

def isStopCodon(c):
	return isFeature(c, 3, ["tga", "taa", "tag"])

def isDonor(c):
	return isFeature(c, 2, static_do)

def isAcceptor(c):
	return isFeature(c, 2, static_ac)

def isDonorSite(p, cds):
	return isDonor(cds[p:p+2])

def isAcceptorSite(p, cds):
	return isAcceptor(cds[p-3:p-1])

def isFeature(c, length, choices):
	if len(c) != length: return False
	return c.lower() in choices


########################################
# Binary algebra
########################################

def binAnd(a, b):
	return a & b

def binAndMulti(l):
	r = l[0]
	for i in l: r = binAnd(r,i)
	return r

def binCompat(bin1, bin2, checkpoints):
	anded = bin1 & bin2
	return not any(anded & c == 0 for c in checkpoints)

def support(a, bases):
        c  = Counter(bases)
        return sum(c[j] for j in set(bases) if binSub(a,j))

def binSub(a, b):
	# Takes binary numbers as input.
        return (a | b) == b


#########################################
# Sequence I/O
#########################################

def readSeqs(path_in, tolist=True):
	seqs = SeqIO.parse(path_in, "fasta")
	if tolist: return list(seqs)
	else: return seqs

def writeSeqs(seqs, path_out):
	SeqIO.write(seqs, path_out, "fasta")

def makeSeq(label, aa):
	return SeqRecord(id=label, name="", description="",  seq=Seq(aa))


#########################################
# Live gtf operations
#########################################

def absFrame(gtfline):
	return (int(gtfline[3]) + int(gtfline[7])) % 3
	
def gtfCompatibleFrame(gtfline1, gtfline2):
	return absFrame(gtfline1) == absFrame(gtfline2)

def getProtoexon(gtf):
	return [int(gtf[3]), int(gtf[4]), int(gtf[7])]

def gtfLength(line):
	return int(line[4]) - int(line[3]) + 1

def gtfLinesOverlap(i,j):
	return (int(i[3]) <= int(j[3]) and int(j[3]) <= int(i[4])) or (int(j[3]) <= int(i[3]) and int(i[3]) <= int(j[4]))

def sameFrame(i, j):
	return (int(i[3]) + int(i[7])) % 3 == (int(j[3]) + int(j[7])) % 3

def overlapInFrame(i,j):
	return gtfLinesOverlap(i,j) and sameFrame(i,j)
	
def gtfLineEquals(gtf1, gtf2):
	return int(gtf1[3])==int(gtf2[3]) and int(gtf1[4])==int(gtf2[4]) and gtf1[0]==gtf2[0]

def safeGtf(gtflines):
	return sorted(deDup([a for a in gtflines if a[2].lower()=="cds"]), key = lambda x: (math.pow(-1,x[6]=="-"))*int(x[3]))

def cleanGtfLive(gtf):
	return [line for line in gtf if not isBadGtf([line])]

def isBadGtf(gtf):
	return any(int(line[3]) >= int(line[4]) for line in gtf)

def getCdsLive(gtf, cds):
	return string.join([cds[int(line[3]) -1: int(line[4])] for line in sorted(gtf, key =  lambda x : int(x[3]))], "")

def translateGtfLive(gtf, cds):
	return translatePart(getCdsLive(gtf, cds))

def translatePart(cdsPart):
	return str(Seq(cdsPart[0:len(cdsPart)-(len(cdsPart)%3)]).translate())

def getEndFrame(gtf):
	return (int(gtf[4]) - int(gtf[3]) - int(gtf[7]) + 1)%3

def getOverlaps(list_gtf):
	overlaps=[]; safe=list(list_gtf)
	for i in list_gtf:
		for j in list_gtf:
			if i == j: continue
			if gtfLinesOverlap(i,j):
				if not i in overlaps: overlaps.append(i)
				if i in safe: safe.remove(i)
	return safe, overlaps

def getNonOverlappingSubsets(overlaps):
	# For a set of possibly overlapping gtflines, returns 
	# all maximally non-overlapping subsets
	d_o={}
	for i, e in enumerate(overlaps): d_o[i] = e
	G=nx.Graph(); keys=d_o.keys()
	for i in keys: G.add_node(i)
	for i in itertools.product(keys, keys):
		if i[0] != i[1] and gtfLinesOverlap(d_o[i[0]], d_o[i[1]]):
			G.add_edge(i[0], i[1])
	H = nx.complement(G)
	C = nx.find_cliques(H)
	return [[d_o[i] for i in k] for k in C]


########################################
# Part group operations
########################################

def groupParts(dict_parts, ranges=True):
	partGraph = nx.Graph()
	for generegion in dict_parts:
		for key in dict_parts[generegion]:
			tag = generegion+"."+str(key)
			if not tag in partGraph.nodes(): partGraph.add_node(tag) 
			part = dict_parts[generegion][key]
			for ogeneregion in dict_parts:
				if generegion == ogeneregion: continue
				for okey in dict_parts[ogeneregion]:
					otag = ogeneregion + "." + str(okey)
					opart = dict_parts[ogeneregion][okey]
					f_o, b_o = overlapProportion(part, opart, ranges)
					if max(f_o, b_o) > 0.333:
						partGraph.add_edge(tag, otag)
	C = list(nx.find_cliques(partGraph))
	groups = [sorted(c) for c in sorted(C, key=lambda x: sorted(x)[0])]
	a = sortPartsList(C)
	return C, a

def sortPartsList(list_c):
	if len(list_c) == 1: return list_c
	G = nx.DiGraph()
	for i, c in enumerate(list_c):
		for j, d in enumerate(list_c):
			if lexLess(c,d): G.add_edge(i,j)
	return [list_c[i] for i in nx.topological_sort(G)]

def lexLess(item1, item2):
	if item1 == item2: return False
	d1 = {}; d2 = {}
	for i in item1:
		m1, o1 = parseId(i); d1[m1] = o1	
	for j in item2:
		m2, o2 = parseId(j); d2[m2] = o2
	isec = set(d1.keys()).intersection(set(d2.keys()))
	if all([d1[i] == d2[i] for i in isec]): return False
	return all([int(d1[i]) <= int(d2[i]) for i in isec])

def overlapProportion(a1, a2, ranges=True):
	a1r = range(min(a1), max(a1)+1) if ranges else a1
	a2r = range(min(a2), max(a2)+1) if ranges else a2
	if len(a1r)*len(a2r) == 0: return 0, 0
	isec = 1.0*len(set(a1r).intersection(a2r))
	return isec/len(a1r), isec/len(a2r)
	

########################################
# Gene part operations
########################################
			
def isDonorPart(gtfline, cds):
	if isTerminalPart(gtfline, cds): return False
	if isDonor(cds[int(gtfline[4]): int(gtfline[4])+2]): return True

def isAcceptorPart(gtfline, cds):
	if isAcceptor(cds[int(gtfline[3])-3: int(gtfline[3])-1]): return True

def isInitialPart(gtfline, cds):
	# Must be zero-framed
	frame = int(gtfline[7])
	if frame != 0: return False
	# If it's a micromicropart, can't tell at this stage if it's a start. return true.
	part = getCdsLive([gtfline], cds)
	if len(part) < 3: return True
	return isStartCodon(part[0:3])

def isTerminalPart(gtfline, cds):
	# Must be zero-ended
	endframe = getEndFrame(gtfline)
	if endframe != 0: return False
	# If it's a micromicropart, can't tell at this stage if it's a stop. return true.
	part = getCdsLive([gtfline], cds)
	if len(part) < 3: return True
	return isStopCodon(part[-3:])

def seqFromPartsList(plist, cds):
	return translateGtfLive([p["gtfline"] for p in plist], cds)

def getPartString(prefix, parts):
	return prefix + string.join([".b" + str(part["gtfline"][3])+"e"+ str(part["gtfline"][4]) for part in parts], "")


########################################
# Id parsing
########################################

def parseId(str_id):
	generegion = re.sub(r"([^.]*)\.(.*)", r"\1", str_id)
	option = int(re.sub(r"([^.]*)\.([0-9]*)", r"\2", str_id))
	return generegion, option

def getGr(str_in):
	return re.sub(r".*(generegion_[0-9]+).*", r"\1", str_in)


########################################
# Alignments
########################################

def align(f_in, f_out, double=False, safety = False):
	alignVanilla(f_in, f_out)
	if safety: safetyCheck(f_in, f_out, f_out)

def alignVanilla(f_in, f_out):
	callFunction("linsi --quiet " + f_in +" > "  + f_out)

def alignRef(f_in, f_out, p_ref, p_refOut):
	callFunction("sed -r \"s/>/>dummy./g\" " + p_ref + " > " + p_refOut)
	callFunction("linsi --quiet --add " + f_in +" " + p_refOut + " > "  + f_out)
	refseqs = [a for a in readSeqs(f_out) if "dummy" in a.id]
	writeSeqs(refseqs, p_refOut)
	cleanDummies(f_out)

def alignSeeded(list_f_in, f_out, prealigned = False, safety = False):
	# Use each of the input fasta files as seeds. Need to double-check
	# that they each have more than one entry.
	functionString = "linsi --quiet "
	f_all = tempfile.mktemp()
	for f_in in list_f_in:
		callFunction("cat " + f_in +" >> " + f_all)
		seqs = readSeqs(f_in)
		if len(seqs) == 0: continue
		if len(seqs) == 1:
			dummy = copy.deepcopy(seqs[0])
			dummy.id = seqs[0].id + ".dummy"
			seqs += [dummy]
		writeSeqs(seqs, f_in)
		if prealigned:
			functionString += " --seed " + f_in
		else:
			f_in_aln = re.sub(r"\fa", r"", f_in) + ".aln"
			align(f_in, f_in_aln)
			functionString += " --seed " + f_in_aln
	functionString += " /dev/null > " + f_out
	callFunction(functionString)
	# Remove the _seed_ prefix and remove any dummy entries
	cleanDummies(f_out)
	if safety: safetyCheck(f_all, f_out, f_out)

def cleanDummies(f_out):
	toclean = readSeqs(f_out)
	cleaned = []
	dummies = []
	for c in toclean:
		if "dummy" in c.id: continue
		s = copy.deepcopy(c)
		s.id = re.sub(r"_seed_", r"", c.id)
		cleaned.append(s)
	writeSeqs(cleaned, f_out)
	writeSeqs(dummies, f_out+".dummies")

def alignmentScore(sequences, omitEmpty = False, scaled = False):
	sequences = [s for s in sequences if not omitEmpty or list(set(str(s.seq))) != ["-"]]
	if not sequences: return 0
	slen   = len(sequences[0])
	if not slen: return 0
	dist   = dict((i, [s.seq[i] for s in sequences]) for i in range(0,slen))
	scores = dict((c,colScore(dist[c])) for c in dist)
	score  = sum(scores.values())
	return score if not scaled else score/(1.0*slen)

def colScore(column):
	counts = Counter(column)
	l = len(static_aa); cl = len(column)
	countsmatrix = np.zeros((l,l));
	for k in set(counts.values()):
		if not k in static_binoms: static_binoms[k] = bn(k,2)
	if not cl in static_binoms: static_binoms[cl] = bn(cl,2)
	# Else-else
	for i,j in itertools.combinations(counts.keys(), 2):
		ipos = static_blospos[i.upper()]
		jpos = static_blospos[j.upper()]
		countsmatrix[ipos, jpos] = counts[i]*counts[j]
	# Self-self
	for i in counts:
		ipos = static_blospos[i.upper()]
		countsmatrix[ipos, ipos] = static_binoms[counts[i]]
	# Don't count things twice.
	scoresmatrix = static_blosmat * countsmatrix
	score = np.sum(scoresmatrix) / static_binoms[cl]
	return score

def flattenAlignment(alignedseqs, gtfs, path_out = "", preSorted = False):
	# FlattenAlignment accepts a list of aligned Seqs, and a dict of loose leaf gtfs.
	# It returns a list of flattened Seqs, and a set of aa ands cds coordinates for the gtfs, and the sorted set of gtfs
	aa=["A","G"]
	flatseqs     = []; goodGtfs     = {}
	locationsAa  = {}; locationsCds = {}
	for s in alignedseqs:
		gtf = gtfs[s.id] if preSorted else safeGtf(gtfs[s.id])
		goodGtfs[s.id] = gtf
		if not gtf: continue
		# The stop codon won't be included in the alignment.
		seqExpanded = string.join([a*3 for a in str(s.seq)], "")
		seqExpanded += "@@@"
		seqExpanded = re.sub(r"(-*)@@@$", r"***\1", seqExpanded)
		# Furthermore, partial codons are possible, so tag some crap on the end for safety
		seqExpanded += "###"
		res = {}; pos = 0; flat = ""
		strand = gtf[0][6]
		for i, line in enumerate(gtf):
			cdslen = int(line[4]) - int(line[3]) + 1
			localpos = 0
			localres = []
			while localpos < cdslen:
				if len(seqExpanded) <= pos or seqExpanded[pos] != "-":
					localres.append(pos)
					flat     += aa[i%2]
					localpos += 1
				else:
					flat += "-"
				pos += 1
			res[i] = sorted(localres)
		outstr = ""
		for qpos in range(len(s)):
			pos=3*qpos
			triplet =flat[pos:pos+3]
			ender = "A" if triplet == "AAA" else ("G" if triplet == "GGG" else ("-" if triplet == "---" else "Q"))
			outstr += ender
		outstr = re.sub(r"Q([Q-]*$)", r"", outstr)
		outstr = outstr + (len(s) - len(outstr))*"-"
		outstr = re.sub(r"(-*)E", r"E\1", outstr + "E")
		outstr = re.sub(r"^(-*)[AG]", r"\1S", outstr)
		t = copy.deepcopy(s)
		t.seq = Seq(outstr)
		flatseqs.append(t)
		locationsCds[s.id] = res
		locationsAa[s.id] = dict((i, [a/3 for a in res[i]]) for i in res)
	if path_out: writeSeqs(flatseqs, path_out)
	return locationsAa, locationsCds, goodGtfs

def cdsToAaLocs(locs):
	return list(set([a/3 for a in locs] + [(a+((3-a)%3))/3 for a in locs]))

def chopAlignment(alnseqs, chopper, negative = False):
	"""chop an alignment based on a list of coordinates
	"""
	if not alnseqs: return []
	if negative: chopper = [a for a in range(len(alnseqs[0])) if not a in chopper]
	return list(chop(alnseqs, chopper))

def chop(alnseqs, chopper):
	res = []
	for a in alnseqs:
		s = copy.deepcopy(a)
		s.seq = Seq(string.join([a[i] for i in chopper if 0 <= i < len(a.seq)], ""))
		res.append(s)
	return res


########################################
# Feature finding
########################################

def findStarts(protoExon, cds, left=True, right=True, central=True, happyWithCentral=True):
	"""Find starts in either direction.
	"""
	pos, boundary, frame = protoExon
	c = cds[pos-1:pos+2]
	cds = cds.lower()
	out = []
	if central and frame == 0 and isStartCodon(c):
		out = [pos]
		if happyWithCentral: return out
	if right: out += walkTermini(cds, pos + frame, 3, lambda c: isStopCodon(c), lambda x: x + 1 > boundary, lambda c: isStartCodon(c), lambda x,c: c[x-1:x+2])
	if left: out += walkTermini(cds, pos + frame, -3, lambda c: isStopCodon(c), lambda x: x - 1 < 0, lambda c: isStartCodon(c), lambda x,c: c[x-1:x+2])
	return out

def findStops(protoExon, cds):
	""" Find the first stop codon in the right hand direction.
	"""
	left, right, frame = protoExon
	pos = right + ((3 - (1 + right - left-frame))%3) -3 
	cds = cds.lower()
	c   = cds[pos-3:pos]
	if isStopCodon(c): return [int(pos)]
	else: return walkTermini(cds, pos, 3, lambda c: False, lambda x: x + 1 > len(cds), lambda c: isStopCodon(c), lambda x,c: c[x-3:x])

def walkTermini(cds, pos, step, escCod, escPos, winFn, cdsSlice):
	while True:
		pos = pos + step
		c = cdsSlice(pos, cds)
		if escCod(c) or escPos(pos): break
		if winFn(c):
			return [int(pos)]
	return []

def containsStop(sequence, frame=0):
	return any(isStopCodon(sequence[i:i+3]) for i in range(frame, len(sequence), 3))

def walkSplice(pos, cds, step, escPos, escCod, nframe, items, direction, cdsSlice, codTest):
	framesDone = []
	while True:
		nframe = (nframe + direction*step) %3 
		pos = pos + step
		if escPos(pos): break
		c=cdsSlice(pos, cds)
		if codTest(c):
			if escCod(pos, nframe): break
			if not nframe in framesDone: items[nframe].append(int(pos))
			framesDone += [nframe]
			if len(set(framesDone)) ==3: break
	return items

def walkSpliceLeft(pos, cds, step, escPos, escCod, nframe, lefts):
	return walkSplice(pos, cds, step, escPos, escCod, nframe, lefts, -1, lambda x,c: c[x-3:x-1], lambda x: isAcceptor(x))

def walkSpliceRight(pos, cds, step, escPos, escCod, nframe, rights):
	return walkSplice(pos, cds, step, escPos, escCod, nframe, rights, 1, lambda x,c: c[x:x+2], lambda x: isDonor(x))

def findSpliceSites(protoExon, cds, left=True, right=True, happyWithCentral = True, first = True):
	# LHS splice sites need to lend the right frame.
	frame = protoExon[2]
	cds = cds.lower()
	lefts = idict([0,1,2], [])
	if left:
		pos, boundary = protoExon[0:2]
		c = cds[pos-3:pos-1]
		good = isAcceptor(c)
		if good: lefts[frame].append(int(pos))
		if (not good) or (good and not happyWithCentral) :
			lefts = walkSpliceLeft(pos, cds, -1, lambda x: x -3 < 0, lambda x, f: containsStop(cds[x-((3-f)%3)-1:x+frame], 0), frame, lefts)
			lefts = walkSpliceLeft(pos, cds, 1, lambda x: x >= boundary, lambda x, f: False, frame, lefts)
	# The frames here are donor frames.
	rights=idict([0,1,2], [])
	if right:
		boundary, pos, oframe = protoExon
		c = cds[pos:pos+2]
		good = isDonor(c)
		frame = (pos - boundary + 1 - protoExon[2]) % 3
		if good: rights[frame].append(int(pos))
		if (not good) or (good and not happyWithCentral):
			rights = walkSpliceRight(pos, cds, -1, lambda x: x < boundary, lambda x,f: False, frame, rights)
			rights = walkSpliceRight(pos, cds, 1, lambda x: x + 1 > len(cds) -1, lambda x,f: containsStop(cds[boundary -((3-oframe)%3)-1:x], 0), frame, rights)
	return lefts, rights

def splitAtStops(ko, string_context, minExon=20):
	happy=[]; waiting=ko
	original = ko[:]
	while waiting:
		waiting_next=[]
		for line in waiting:
			cds       =string_context[int(line[3])-1:int(line[4])]
			frame     =int(line[7])
			scanstart =frame
			scanend   =frame+3
			safe      = True
			while scanend <= len(cds):
				codon=cds[scanstart:scanend]
				if isStopCodon(codon):
					newline1     = line[:]
					newline2     = line[:]
					newline1[4]  = int(line[3]) + scanstart - 1
					newline2[3]  = int(line[3]) + scanend
					newline2[7]  = 0
					waiting_next += [newline1, newline2]
					safe = False
					break
				scanstart += 3
				scanend += 3
			if safe: happy.append(line)
		waiting = waiting_next
	# Throw out any bad exons. Also throw out any of the new split bits
	# if they are tiny: these bits are probably rubbish.
	happy = [i for i in happy if int(i[4]) >= int(i[3])]
	happy = [i for i in happy if [a for a in original if gtfLineEquals(a, i)] or int(i[4])-int(i[3]) >= minExon]
	return happy
	

##########################################
# File-based gtf/fasta operations
##########################################

def translateGtfToFile(gtf, cds, tag, path_out):
	aaParts = translateGtf(gtf, cds, tag)
	if aaParts: writeSeqs(aaParts, path_out)
	else: blankFile(path_out)

def translateGtf(gtf, cds, tag):
	for line in gtf:
		cdspart = cds[int(line[3])+int(line[7]) -1 : int(line[4])]
		aapart  = translatePart(cdspart)
		if aapart: yield makeSeq(tag + ".b" + str(line[3]) + "e" + str(line[4]), aapart)

def framify(path_gtf):
	"""Framify assumes that the first exon is in frame 0
	"""
	data = readCsv(path_gtf)
	geneNames = list(set([a[1] for a in data]))
	genes = [[a for a in data if a[1] == gene] for gene in geneNames]
	results = []
	for g in genes:
		gs=sorted([e for e in g if e[2].lower()=="cds"], key = lambda x: int(x[3]))
		lastframe=0
		for cds in gs:
			cdsc=cds[:]
			cdsc[7]=lastframe
			results.append(cdsc)
			nextframe=(3-(int(cds[4])-int(cds[3])+1-lastframe)) % 3
			lastframe=nextframe
	writeCsv(results, path_gtf)

def splitToAbsFrame(gtf):
	gtf_framed = idict([0,1,2], [])
	for line in gtf:
		start	 = int(line[3])
		frameRel = int(line[7])
		frameAbs = (start + frameRel) % 3
		gtf_framed[frameAbs].append(line)
	return gtf_framed

def mergeFriendly(path_gtf, path_out, sort=True):
	"""Merge elements of a gtf file, making sure only to merge
	   which are frame-compatible. Requires frame info.
	"""
	gtf = readCsv(path_gtf)
	gtf_framed = splitToAbsFrame(gtf)
	out = []
	for i in gtf_framed:
		if not any(gtf_framed[i]): continue
		path_gtftmp = path_gtf+".f"+str(i)+".tmp"
		writeCsv(gtf_framed[i], path_gtftmp)
		merged = [a.split("\t") for a in grabLines("sort -k4,4n " + path_gtftmp + " | bedtools merge -i -")]
		model = gtf_framed[i][0]
		for line in merged:
			newline = model[:]
			newline[1] = "."
			newline[8] = "source_id \"any\""
			newline[3] = int(line[1]) + 1
			newline[4] = line[2]
			newline[7] = (i - int(newline[3])) % 3
			out.append(newline)
	if sort: out = sorted(out, key = lambda x: int(x[3]))
	writeCsv(out, path_out)

def cleanGtf(path_gtfO, path_gtf):
	callFunction("sort -u " + path_gtfO + " | sort -k4,4n > " + path_gtf)

def replaceSeq(seqObject, newSeq):
	s = copy.deepcopy(seqObject)
	s.seq = newSeq
	return s

def fetchAa(path_gtfIn, path_genome, path_cdsFastaOut, path_aaFastaOut, int_translationTable):
	"""Fetch cds and translate to aa
	"""
	fetchCds(path_gtfIn, path_genome, path_cdsFastaOut, "CDS")
	protSequences = [replaceSeq(s, s.seq.translate(table=int_translationTable)) for s in readSeqs(path_cdsFastaOut, False)]
	writeSeqs(protSequences, path_aaFastaOut)

def fetchCds(path_gtfIn, path_genome, path_cdsFastaOut, str_token):
	"""Fetch the cds fasta sequences for the given gtfs.
	"""
	callFunction("infile=" +  path_gtfIn+ "; outfile=" + path_cdsFastaOut + "; genome=" + path_genome + "; token=" + str_token + """;
		tf=`mktemp -d`
		gtfCds="$tf/gtfCds"
		gtfBed="$tf/gtfBed"

		#Prepare the gtf
		grep -vP "^$" $infile | awk -v token="$token" '$3==token' > $gtfCds
		cut -f1-8 $gtfCds > $gtfBed.1
		sed -r "s/.*transcript_id[ =]\\"?([^\\";]*)\\"?;?.*/\\1/g" $gtfCds > $gtfBed.2
		paste $gtfBed.1 $gtfBed.2 | perl -ne 'chomp; @l=split; printf "$l[0]\\t%s\\t$l[4]\\t$l[8]\\t.\\t$l[6]\\n", $l[3]-1' | sort -u | sort -k1,1V -k2,2n > $gtfBed

		#Negative strand
		awk '$6=="-"' $gtfBed > $gtfBed.neg
		bedtools getfasta -name -s -fullHeader -fi $genome -fo $gtfBed.neg.tab -bed $gtfBed.neg -tab
		tac $gtfBed.neg.tab | awk '{a[$1]=a[$1]""$2} END {for (i in a) {print ">"i"\\n"a[i]}}' > $gtfBed.neg.fa

		#Then positive strand
		awk '$6=="+"' $gtfBed > $gtfBed.pos
		#cat $gtfBed.pos
		bedtools getfasta -name -s -fullHeader -fi $genome -fo $gtfBed.pos.tab -bed $gtfBed.pos -tab
		cat $gtfBed.pos.tab | awk '{a[$1]=a[$1]""$2} END {for (i in a) {print ">"i"\\n"a[i]}}' > $gtfBed.pos.fa

		cat $gtfBed.pos.fa $gtfBed.neg.fa | sed -r "s/^>(.*)$/£££>\\1###/g" | sed -r \"s/$/###/g\" | tr '\\n' ' ' | sed -r "s/£££/\\n/g" | sed -r "s/### //g" | grep -v XXX | grep -v "\*[A-Z]" | grep -v "###$" | sed -r "s/###/\\n/g" | grep -vP "^$" > $outfile

		rm -r $tf
		""")

def getBase(path_gtf, int_slop, dict_chrSizes, path_base, path_baseRel, path_baseTight, sequenceId):
	"""Get the base for a gtf file, with margins of a specified size
	"""
	gtf = readCsv(path_gtf)
	# Use the data from the first line, but using min and max location coordinates.
	geneEnd		= max([int(b[4]) for b in gtf])
	geneStart	= min([int(b[3]) for b in gtf])
	# Grab the gene and transcript names where possible
	str_geneId       = re.sub(r'.*gene_id \"([^\"]*)\".*', r'\1', [a[8] for a in gtf if "gene_id" in a[8]][0])
	str_transcriptId = re.sub(r'.*transcript_id \"([^\"]*)\".*', r'\1', [a[8] for a in gtf if "transcript_id" in a[8]][0])
	# Slop the entry, ensuring that the slopped coordinates are still viable within the genome
	entry	  = gtf[0][0:9]
	chrName   = entry[0]
	chrSize   = dict_chrSizes[chrName]
	descr     = "transcript_id \"" + str_transcriptId + ".t1\"; gene_id \"" + str_geneId + "\""
	entry     = makeGtfLine(entry, max(int(geneStart) - int_slop, 1), min(int(geneEnd) + int_slop, chrSize), chrName, "", descr, "base")
	# Construct both slopped and unslopped versions of the base.
	entry.append(sequenceId)
	entry_tight = makeGtfLine(entry, geneStart, geneEnd, "", "", "", "base_tight")
	writeCsv([entry], path_base)
	writeCsv([entry_tight], path_baseTight)
	# Construct the relative gtf for the slopped base
	entryRel = makeGtfLine(entry, 1, entry[4] - entry[3] +1, "relative", "", "", "")
	writeCsv([entryRel], path_baseRel)

def makeGtfLine(refline, pstart, pend, chrom, strand, descr, tag):
	line = refline[:]
	if chrom:  line[0] = chrom
	if tag:    line[2] = tag
	if pstart: line[3] = pstart
	if pend:   line[4] = pend
	if strand: line[6] = strand
	if descr:  line[8] = descr 
	return line

def toggleRelative(path_gtfIn, path_gtfOut, path_gtfBase, cdsOnly=False, tag="", dirn = "fwd"):
	# Read in the data
	data_gtfIn   = readCsv(path_gtfIn)
	data_gtfBase = readCsv(path_gtfBase)[0]
	if cdsOnly: data_gtfIn = [a for a in data_gtfIn if a[2] == "CDS"]
	strand  = data_gtfBase[6]
	gStartO = int(data_gtfBase[3])
	gEndO   = int(data_gtfBase[4])
	chrom   = data_gtfBase[0]
	descr   = data_gtfBase[8]
	# Go through each line and flip it.
	data_gtfFlipped = []
	for oline in data_gtfIn:
		if dirn == "fwd":
			lEnd = gEndO - int(oline[4]) + 1 if strand == "-" else int(oline[3]) - gStartO + 1
			rEnd = gEndO - int(oline[3]) + 1 if strand == "-" else int(oline[4]) - gStartO + 1
		else:
			lEnd = gEndO - int(oline[4]) + 1 if strand == "-" else gStartO + int(oline[3]) - 1
			rEnd = gEndO - int(oline[3]) + 1 if strand == "-" else gStartO + int(oline[4]) - 1
		line = makeGtfLine(oline, lEnd, rEnd, chrom, "+" if dirn =="fwd" else strand, descr, "")
		if tag != "": line[1]=tag
		data_gtfFlipped.append(line)
	writeCsv(data_gtfFlipped, path_gtfOut)


########################################
# Setup functions
########################################

def getGenomeInfo_indi(path_genome):
	fas = pyfaidx.Fasta(path_genome)
	return dict((record.name, len(record)) for record in fas)

def getGenomeInfo(dict_seqInfo):
	"""Get the sizes of all chromosomes in each genome
	"""
	genomeFiles = list(set([dict_seqInfo[b]["genome"] for b in dict_seqInfo]))
	genomeSizes = dict((path_genome, getGenomeInfo_indi(path_genome)) for path_genome in genomeFiles)
	return genomeSizes

def readInputLocations(path_inf, path_ref):
	"""Read CSV file containing the locations for the sequence files, etc.
	   Each row is a sequence. The columns are [gtf, genome].
	"""
	sprint("loading and checking input data locations...")
	data = readCsv(path_inf, ignoreHash = True)
	dict_seqInfo = readInputFile(readCsv(path_inf, ignoreHash = True), {}, False)
	if path_ref: dict_seqInfo = readInputFile(readCsv(path_ref, ignoreHash = True), dict_seqInfo, True)
	return dict_seqInfo

def processCodonOption(co, length):
	s_co = co.split(",")
	for a in s_co:
		a = re.sub(r" ", r"", a)
		if not len(a) == length or any(not s.lower() in "acgtu" for s in a):
			sys.exit("Invalid start codon/splice site choice: " + str(a))
		yield re.sub(r"u", r"t", a.lower())

def processCodonOptions(do, ac, sc):
	s_do = list(set(processCodonOption(do, 2)))
	s_ac = list(set(processCodonOption(ac, 2)))
	s_sc = list(set(processCodonOption(sc, 3)))
	return s_do, s_ac, s_sc

def readInputFile(data, dict_seqInfo, keep = False):
	nextIndex = len(dict_seqInfo)
	for i, line in enumerate(data):
		# Check files exist
		path_gtf    = checkFileExists(line[0])
		path_genome = checkFileExists(line[1])
		# TODO CUPCAKES:
		# need to check chromosomes, coordinates, and that each gtf refers to a single gene
		# Also need to check that each genome id is unique.
		seqId = "sequence_" + str(nextIndex+i)
		dict_seqInfo[seqId] = {"id": seqId, "gtf": path_gtf, "genome": path_genome, "gtfId": os.path.basename(path_gtf), "genomeId": os.path.basename(path_genome), "keep": keep}
	return dict_seqInfo

def prepareOutputFolder(path_outDir):
	path_outDir = makeIfAbsent(path_outDir if path_outDir else generateDir())
	return makeIfAbsent(path_outDir + "/results"), makeIfAbsent(path_outDir + "/working")

def generateDir():
	return tempfile.mkdtemp(prefix="OMGeneOut_" + datetime.datetime.now().strftime("%y%m%d") + "_RunId_", dir=".")

def prepareSequences(path_wDir, dict_seqInfo, dict_genomeInfo, int_numCores, int_slopAmount):
	pool = multiprocessing.Pool(int_numCores)
	path_sDir = makeIfAbsent(path_wDir+"/sequences")
	for seqId, ds in dict_seqInfo.items():
		# Get the base gtf and its CDS. Organise them nicely into folders.
		path_sDir_seqid    = makeIfAbsent(path_sDir +"/" + seqId)
		ds["gtfBase"]      = path_base	          = path_sDir_seqid + "/" + seqId + ".base"
		ds["gtfBaseTight"] = path_baseTight       = path_base + ".tight"
		ds["gtfBaseRel"]   = path_baseRel         = path_base + ".relative"
		ds["cdsBase"]      = path_cdsFastaOutBase = path_base + ".cds.fasta"
		# Get gtf files
		path_gtfO	= ds["gtf"]
		path_genome	= ds["genome"]
		dict_chrSizes	= dict_genomeInfo[path_genome]
		# Clean up the gtf file and move it to the working directory.
		ds["gtf"] = path_gtf = path_sDir_seqid + "/"+seqId+".gtf"
		cleanGtf(path_gtfO, path_gtf)
		async(pool, grabBases, args=(path_gtf, int_slopAmount, dict_chrSizes, path_base, path_baseRel, path_baseTight, path_genome, path_cdsFastaOutBase, seqId))
		# Get the cds and aa for the genes
		ds["aa"]  = path_aaFastaOut  = path_sDir_seqid + "/" + seqId + ".aa.fasta"
		ds["cds"] = path_cdsFastaOut = path_sDir_seqid + "/" + seqId + ".cds.fasta"
		async(pool, fetchAa,  args=(path_gtf, path_genome, path_cdsFastaOut, path_aaFastaOut, 1))
	pool.close()
	pool.join()
	for ds in dict_seqInfo.values():
		ds["sequence_aa"] = str(readSeqs(ds["aa"])[0].seq).strip("-")
		ds["sequence_cds"] = str(readSeqs(ds["cds"])[0].seq).strip("-")

def grabBases(path_gtf, int_slopAmount, dict_chrSizes, path_base, path_baseRel, path_baseTight, path_genome, path_cdsFastaOutBase, sequenceId):
	getBase(path_gtf, int_slopAmount, dict_chrSizes, path_base, path_baseRel, path_baseTight, sequenceId)
	fetchCds(path_base, path_genome, path_cdsFastaOutBase, "base")


def writeGeneRegionFile(line, generegion, path_generegion):
	nline = [line[0], "orthofxer", "generegion", line[1], line[2],\
			".", line[3], ".", "transcript_id \"" + generegion + "\"; gene_id \"" + generegion + "\"", generegion]
	writeCsv([nline], path_generegion)
	return path_generegion

def sequenceMatchesRegion(seqGtf, regionGtf):
	return [l for l in grabLines("bedtools intersect -a "+regionGtf+" -b "+seqGtf) if ''.join(l).strip()]

def prepareGeneregions(dict_seqInfo, dict_genomeInfo, path_wDir, int_numCores, int_slopAmount):
	dict_generegions = {}
	path_mDir = makeIfAbsent(path_wDir + "/generegions")
	path_gDir = makeIfAbsent(path_mDir + "/genomes")
	pool      = multiprocessing.Pool(int_numCores)
	for path_genome in list(set([ a["genome"] for a in dict_seqInfo.values()])):
		# Join all tight base files together and then perform a bedtools merge
		genomeId		= os.path.basename(path_genome)
		bases_tight		= [a["gtfBaseTight"] for a in dict_seqInfo.values() if a["genome"] == path_genome]
		path_allBaseLoci	= path_gDir + "/" + genomeId + ".bases"
		path_mergedBaseLoci	= path_gDir + "/" + genomeId + ".bases.merged"
		concatFiles(bases_tight, path_allBaseLoci)
		callFunction("sort -k1,1V -k4,4n " + path_allBaseLoci + " | bedtools merge -s -i - > " + path_mergedBaseLoci)
		sequences = [copy.deepcopy(a) for a in dict_seqInfo.values() if a["genome"] == path_genome]
		for line in readCsv(path_mergedBaseLoci):
			# Initialise the generegion in the holder dict
			generegion = "generegion_" + str(len(dict_generegions.keys()))
			dg = dict_generegions[generegion] = {}
			dg["mDirl"]  = path_mDir_l = makeIfAbsent(path_mDir + "/" + str(generegion))
			dg["genome"] = path_genome
			dg["generegion"] = generegion
			# The path of the generegion base descriptor
			path_generegion	= writeGeneRegionFile(line, generegion, path_mDir_l + "/" + str(generegion) + ".gtf")
			# Figure out which sequences correspond to this generegion.
			matchingseqs=[]
			for s in sequences:
				if sequenceMatchesRegion(s["gtf"], path_generegion):
					dict_seqInfo[s["id"]]["generegion"] = generegion
					matchingseqs.append(s["id"])
			dg["sequences"] = matchingseqs
			# Prepare all the lovely file names
			path_generegionStub = path_mDir_l + "/" + str(generegion)
			dg["base"]      = path_generegionSlopped      = path_generegionStub + ".base"
			dg["baseRel"]   = path_generegionSloppedRel   = path_generegionStub + ".baseRel"
			dg["baseTight"] = path_generegionSloppedTight = path_generegionStub + ".baseTight"
			dg["cdsBase"]   = path_generegionCds	      = path_generegionStub + ".base.cds.fasta"
			dict_chrSizes = dict_genomeInfo[path_genome]
			grabBases(path_generegion, int_slopAmount, dict_chrSizes, \
                                                      path_generegionSlopped, path_generegionSloppedRel, path_generegionSloppedTight, \
                                                      path_genome, path_generegionCds, generegion)
			# Sort out any "keep" gene models.
			keepseqs   = [sequence for sequence in dg["sequences"] if dict_seqInfo[sequence]["keep"]]
			changeseqs = [sequence for sequence in dg["sequences"] if not sequence in keepseqs]
			if len(keepseqs) == 0:
				dg["keep"] = False
				continue
			else:
				ccc = copy.deepcopy(dg)
				if len(changeseqs) == 0:
					del dict_generegions[generegion]
				else:
					dict_generegions[generegion]["sequences"] = changeseqs
					dict_generegions[generegion]["keep"]      = False
				for k in keepseqs:
					ngeneregion = "generegion_" + str(len(dict_generegions.keys()))
					dict_generegions[ngeneregion] = copy.deepcopy(ccc)
					dict_generegions[ngeneregion]["sequences"] = [k]
					dict_generegions[ngeneregion]["keep"]      = True
	pool.close()
	pool.join()
	# Add the cds string as a data item.
	for dg in dict_generegions.values():
		dg["cdsBase_data"] = str(readSeqs(dg["cdsBase"])[0].seq)
		dg["baselen"] = len(dg["cdsBase_data"])
	return dict_generegions

def writeOutput(dict_generegions, dict_seqInfo, res, path_resultsDir):
	"""Write out all the various bits of output
	"""
	path_gtfDir = makeIfAbsent(path_resultsDir+"/gtf/")
	path_aaDir  = makeIfAbsent(path_resultsDir+"/aa/")
	path_cdsDir = makeIfAbsent(path_resultsDir+"/cds/")
	allAa = []
	for generegion, dg in dict_generegions.items():
		for sequence in dg["sequences"]:
			gtfId  = dict_seqInfo[sequence]["gtfId"]
			gtfOut = path_gtfDir + gtfId + ".fixed.gtf"
			aaOut  = path_aaDir + gtfId + ".fixed.aa.fa"
			cdsOut = path_cdsDir + gtfId + ".fixed.cds.fa"
			# Prepare gtf
			gtf = [part["gtfline"] for  part in res[generegion]] if generegion in res else []
			writeCsv(gtf, gtfOut + ".tmp")
			toggleRelative(gtfOut+".tmp", gtfOut, dict_generegions[generegion]["base"], False, "", "rev")
			os.remove(gtfOut + ".tmp")
			dict_seqInfo[sequence]["res_gtf"] = readCsv(gtfOut)
			dict_seqInfo[sequence]["res_gtfpath"] = gtfOut
			# Prepare cds
			cdsString = string.join([dg["cdsBase_data"][int(line[3])-1: int(line[4])] for line in gtf], "")
			dict_seqInfo[sequence]["res_cds"] = cdsString
			writeSeqs([makeSeq(gtfId, cdsString)], cdsOut)
			# Prepare aa
			aaString = translatePart(cdsString)
			dict_seqInfo[sequence]["res_aa"] = aaString
			writeSeqs([makeSeq(gtfId, aaString)], aaOut)
			allAa.append(makeSeq(gtfId, aaString))
	writeSeqs(allAa, path_resultsDir+"/all.fa")
	align(path_resultsDir+"/all.fa", path_resultsDir+"/all.aln")
	return dict_generegions

def relativiseSequences(dict_generegions, dict_seqInfo):
	for generegion, dg in dict_generegions.items():
		path_gtfBase = dg["base"]
		path_mDir_l  = dg["mDirl"]
		dg["boundaries"] = path_boundaries = path_mDir_l + "/" + generegion + ".boundaries"
		with open(path_boundaries, "w") as o:
			for seqId in dg["sequences"]:
				path_gtf    = dict_seqInfo[seqId]["gtf"]
				dict_seqInfo[seqId]["gtfRel"] = path_gtfRel = path_mDir_l + "/" + seqId + ".rel." + generegion + ".gtf"
				toggleRelative(path_gtf, path_gtfRel, path_gtfBase, True, seqId, "fwd")
				# Can't guarantee that the gff will have frame information in it, so provide that information.
				framify(path_gtfRel)
				with(open(path_gtfRel, "r")) as b: o.write(b.read())


########################################
# Exonerate functions
########################################

def goExonerate(path_wDir, dict_generegions, dict_seqInfo, int_numCores):
	sprint("Performing first round exonerate...")
	# Prepare output directories
	path_eDir  = makeIfAbsent(path_wDir + "/exonerate")#qe
	path_eDir1 = makeIfAbsent(path_eDir + "/firstRound")#qe
	path_eDir2 = makeIfAbsent(path_eDir + "/secondRound")#qe
	path_eDir3 = makeIfAbsent(path_eDir + "/thirdRound")#qe
	# Perform first round exonerate
	exonerateFirstRound(dict_generegions, dict_seqInfo, int_numCores, path_eDir1)#qe
	# Run exonerate again using all the new potential gene parts.
	sprint("Performing second round exonerate...")
	exonerateSecondRound(dict_generegions, int_numCores, path_eDir2)#qe
	combineExonerate1and2(dict_generegions, path_eDir)#qe
	# Run a third time.
	sprint("Performing third round exonerate...")
	exonerateThirdRound(dict_generegions, dict_seqInfo, path_eDir3)#qe
	sprint("done getting options")

def exonerate(path_base, list_allAaBase, path_exonerateOut):
	exonerateLines = []
	for aa in list_allAaBase:
		aaName = os.path.basename(aa)
		capture = subprocess.Popen("exonerate --showalignment no --showtargetgff --frameshift -10000 --model protein2genome " \
				+ aa + " " + path_base + "  | grep \"exonerate:protein2genome:local\" | grep -P \"(cds)|(gene)\" " \
				+"| grep -v \"#\" | sed -r \"s/.*exonerate:protein2genome:local/relative\\t"+aaName+"/g\" " \
				+"| sed -r \"s/cds/CDS/g\"", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
		stdout = [x for x in capture.stdout]
		stderr = [x for x in capture.stderr]
		if not stderr: exonerateLines += [a.strip("\n").split("\t") for a in stdout]
	geneNum = 0; result = []
	for line in exonerateLines:
		if line[2]=="gene":
			geneNum += 1
			source = re.sub(r".*sequence ([^ ;]*) ;.*", r"\1", line[8])
		nline = line[:]
		nline[8] = source
		nline[1] = "gene_" + str(geneNum) + "."+line[1]
		if line[2].lower() == "cds": result.append(nline)
	writeCsv(result, path_exonerateOut)

def implementExonerate(path_cdsBase, list_allAa, path_exonerateOut):
	exonerate(path_cdsBase, list_allAa, path_exonerateOut)
	framify(path_exonerateOut)

def exonerateFirstRound(dict_generegions, dict_seqInfo, int_numCores, path_eDir):
	sprint("Exonerating and piling sequences against sequence regions...")
	for generegion, dg in dict_generegions.items():
		dg["exonerated"] = path_exonerateOut = path_eDir + "/" + generegion + ".exonerate1"
		# There should only be one sequence in each keep generegion.
		if dg["keep"]:
			copyfile(dict_seqInfo[dg["sequences"][0]]["gtfRel"], path_exonerateOut)
		else:
			list_allAa   = [dict_seqInfo[a]["aa"] for a in dict_seqInfo if not a in dg["sequences"]]
			path_cdsBase = dg["cdsBase"]
			path_baseRel = dg["baseRel"]
			implementExonerate(path_cdsBase, list_allAa, path_exonerateOut)
		# Merge the exonerate output
		dg["exonerated_merged"] = path_exonerateOutMerged = path_exonerateOut +".merged"
		mergeFriendly(path_exonerateOut, path_exonerateOutMerged)
		# Deal with any in-frame stops
		gtf = readCsv(path_exonerateOutMerged)
		cds = dg["cdsBase_data"]
		if not dg["keep"]: gtf = splitAtStops(gtf, cds)
		writeCsv(gtf, path_exonerateOutMerged)
		# Split out the individual aa parts for each gtf line.
		dg["path_genePartsFasta"] = path_parts = path_exonerateOutMerged + ".parts"
		translateGtfToFile(gtf, cds, generegion, path_parts)

def exonerateSecondRound(dict_generegions, int_numCores, path_eDir2):
	for generegion, dg in dict_generegions.items():
		dg["exonerated2"] = path_secondRound = path_eDir2 + "/" + generegion + ".exonerated2"
		path_cds = dg["cdsBase"]
		if dg["keep"]:
			copyfile(dg["exonerated"], path_secondRound)
		else:
			othergeneParts = [dict_generegions[ogeneregion]["path_genePartsFasta"] for ogeneregion in dict_generegions if ogeneregion != generegion]
			implementExonerate(path_cds, othergeneParts, path_secondRound)
			path_secondRoundMerged = dg["exonerated2_merged"] = path_secondRound +".merged"
			mergeFriendly(path_secondRound, path_secondRoundMerged)

def exonerateThirdRound(dict_generegions, dict_seqInfo, path_eDir3):
	"""This function aims to clean up any regions in the exonerated output
	   that overlap but are in incompatible frames.
	"""
	# Generate possible versions and choose between them.
	for generegion, dg in dict_generegions.items():
		# Get possible combinations of bits of gene from the exonerate output.
		options = getOptions(readCsv(dg["exoneratedAllMerged"]))
		# Pick the best of those options and write them to file.
		dg["path_options"] = makeIfAbsent(dg["mDirl"] + "/options/") + "/" + generegion + ".options.fasta"
		chooseOptions(dict_generegions[generegion], options)
	# Run exonerate on the cleaned up genes again.
	for generegion, dg in dict_generegions.items():
		options   = flattenLol([[o["fasta"] for o in dict_generegions[omet]["options"].values()] for omet in dict_generegions if not omet == generegion ])
		path_base = dict_generegions[generegion]["cdsBase"]
		path_ex3  = path_eDir3 +"/"+generegion+".exonerated3"
		# Implement exonerate with this new set of data
		if not dg["keep"]: implementExonerate(path_base, options, path_ex3)
		else: copyfile(dg["exonerated_merged"], path_ex3)
		gtf = readCsv(path_ex3)
		cds = dg["cdsBase_data"]
		# Be careful about stops.
		if not dg["keep"]: gtf = splitAtStops(gtf, cds)
		writeCsv(gtf, path_ex3)
		path_all = dg["exoneratedAllMerged"]
		# Add them all up.
		allgtf   = readCsv(path_all) + readCsv(path_ex3)
		writeCsv(allgtf, path_all)
		mergeFriendly(path_all, path_all)
	# Choose between the resulting options.
	for generegion in dict_generegions:
		# Add in any bits of the input that do not overlap the exonerate-found regions.
		gtf_o = safeGtf(flattenLol([readCsv(dict_seqInfo[oseq]["gtfRel"]) for oseq in dict_generegions[generegion]["sequences"]]))
		gtf   = readCsv(dict_generegions[generegion]["exoneratedAllMerged"])
		gtf   = addNonOverlappers(gtf, gtf_o)
		# Hack on lazy code
		for i in gtf: i[0] = "relative"
		options = getOptions(gtf)
		chooseOptions(dict_generegions[generegion], options)

def addNonOverlappers(list_gtfLines, gtf_o):
	res = copy.deepcopy(list_gtfLines)
	framed_r = splitToAbsFrame(res)
	framed_o = splitToAbsFrame(gtf_o)
	for frame in framed_o:
		for line in framed_o[frame]:
			if not any([gtfLinesOverlap(line, i) for i in framed_r[frame]]):
				res += [line]
	return safeGtf(res)

def combineExonerate1and2(dict_generegions, path_eDir):
	for generegion, dg in dict_generegions.items():
		# Set up files
		dg["exoneratedAll"]       = path_all       = path_eDir + "/"+generegion + ".all"
		dg["exoneratedAllMerged"] = path_allMerged = path_all + ".merged"
		# Perform the merge
		if dg["keep"]:
			copyfile(dg["exonerated_merged"], path_all)
		else:
			concatFiles([dg["exonerated_merged"], dg["exonerated2_merged"]], path_all)
		mergeFriendly(path_all, path_allMerged)
		# Split into parts and split at stops.
		gtf = readCsv(path_allMerged)
		cds = dg["cdsBase_data"]
		if not dict_generegions[generegion]["keep"]:
			gtf = splitAtStops(gtf, cds)
		writeCsv(gtf, path_allMerged)
		translateGtfToFile(gtf, cds, generegion, path_allMerged+".parts")

def getOptions(list_gtfLines):
	xsafe, goverlaps = getOverlaps(list_gtfLines)
	options = []
	if not goverlaps: return [xsafe]
	# The first set of options should involve simply removing one or other of the offending overlappers.
	goodFromOverlapping = getNonOverlappingSubsets(goverlaps)
	for suboption in goodFromOverlapping:
		options.append(cleanGtfLive(xsafe + suboption))
	# The second set of options should involve chopping lines along their join point, and throwing away parts that end up too tiny.
	xoverlaps = goverlaps[:]
	for pair in itertools.combinations(goverlaps, 2):
		l_safe, l_overlaps = getOverlaps(pair)
		if l_overlaps:
			firstFirst = int(l_overlaps[0][3]) <= int(l_overlaps[1][3])
			first = l_overlaps[int(not firstFirst)]
			second = l_overlaps[int(firstFirst)]
			# Only want to consider things of this form:
			# ----------------
			#	   -----------------
			if int(first[4]) >= int(second[4]): continue
			# Find the centroid
			overlapStart = int(second[3])
			overlapEnd   = int(first[4])
			breakpoint   = (overlapStart + overlapEnd)/2
			# Generate new gene parts
			newfirst = first[:]
			newfirst[4] = breakpoint - 1
			newsecond = second[:]
			newsecond[3] = breakpoint
			newsecond[7] = (3 - (breakpoint - int(second[3]) - int(second[7]))) % 3
			# Don't want to include really tiny ones!!
			if min(gtfLength(newfirst), gtfLength(newsecond)) < 20: continue
			xoverlaps+=[newfirst, newsecond]
	goodFromOverlapping = getNonOverlappingSubsets(xoverlaps)
	for suboption in goodFromOverlapping:
		options.append(cleanGtfLive(xsafe + suboption))
	# The final set of options should involve simply splitting up into chunks, and throwing away parts that end up too tiny.
#	xoverlaps=goverlaps[:]
#	print "==============4"
#	for pair in itertools.combinations(goverlaps, 2):
#		l_safe, l_overlaps = getOverlaps(pair)
#		if l_overlaps:
#			firstFirst = int(l_overlaps[0][3]) <= int(l_overlaps[1][3])
#			first = l_overlaps[int(not firstFirst)]
#			second = l_overlaps[int(firstFirst)]
#			firstchunk1 = first[:]
#			firstchunk1[4] = int(second[3]) - 1
#			firstchunk2 = first[:]
#			firstchunk2[3] = int(second[3])
#			firstchunk2[7] = (3 - (int(second[3]) - int(first[3]) - int(first[7]))) % 3
#			secondchunk1 = second[:]
#			secondchunk1[4] = int(first[4])
#			secondchunk2 = second[:]
#			secondchunk2[3] = int(first[4]) + 1
#			secondchunk2[7] = (3 - (int(first[4]) + 1  - int(second[3]) - int(second[7]))) % 3
#			for opt in [firstchunk1, firstchunk2, secondchunk1, secondchunk2]:
#				if gtfLength(firstchunk1) >= 20: xoverlaps += [opt]
#	goodFromOverlapping = getNonOverlappingSubsets(xoverlaps)
#	for suboption in goodFromOverlapping:
#		options.append(cleanGtfLive(xsafe + suboption))
	return options

def chooseOptions(dg, final):
	dg["options"] = {}
	for i,ens in enumerate(final):
		# Construct entries for each option.
		e = sorted(ens,key=lambda x: int(x[3]))
		path_optionOut   = dg["path_options"] + ".option_" + str(i)
		writeCsv(e, path_optionOut)
		# Grab the fasta sequences.
		path_optionFasta = path_optionOut+".aa.fasta"
		dg["options"][i] = {"gtf": path_optionOut, "fasta": path_optionFasta, "parts": {}}
		cds = dg["cdsBase_data"]
		with open(path_optionFasta, "w") as f:
			f.write(">" +dg["generegion"] + "."+ str(i) + "\n")
			for j, line in enumerate(e):
				cdspart = cds[int(line[3])+int(line[7]) -1 : int(line[4])]
				aapart = translatePart(cdspart)
				f.write(aapart)
				dg["options"][i]["parts"][j] = {"part": aapart, "partlength": len(aapart), "gtfline": line}
			f.write("\n")		

########################################
# Stats
########################################

def getStats(dict_generegions, dict_seqInfo, res, path_resultsDir, path_wDir):
	# Get a gene tree.
	path_sDir = makeIfAbsent(path_resultsDir+"/stats/")
	path_tDir = makeIfAbsent(path_resultsDir+"/trees/")
	#callFunction("iqtree -redo -s "+path_resultsDir+"/all.aln -pre "+path_tDir+"/all.aln.tree -m TEST -quiet")
	# Get original fasta alignment score. If there are any empty sequences, take them out.
	oseqs = readSeqs(path_wDir+"/or.aln")
	oscore = alignmentScore(oseqs, omitEmpty = True)
	# Get new alignment score. If there are any empty sequences, take them out.
	nseqs = readSeqs(path_resultsDir+"/all.aln")
	nscore = alignmentScore(nseqs, omitEmpty = True)
	# Write global stats to file:
	path_statsGlobal = path_sDir + "/global.stats"
	with open(path_statsGlobal, "w") as f:
		f.write("oldscore\tnewscore\n")
		f.write(str(oscore)+"\t"+str(nscore))
	statsSh=tempfile.mktemp()+".sh"
	statsScript(statsSh)
	# Write detailed stats to file...
	for sequence in dict_seqInfo:
		ogtf = dict_seqInfo[sequence]["gtf"]
		ngtf = dict_seqInfo[sequence]["res_gtfpath"]
		#sprint("comparing " + ogtf + " and " + ngtf)
		sprint("Gene coordinates written to " + ngtf)
		callFunction("bash "+statsSh+" " + ogtf + " " + ngtf +" > " + path_sDir +"/"+dict_seqInfo[sequence]["gtfId"]+".stats")

def statsScript(tmpf):
	with open(tmpf,"w") as f:
		f.write('#!/bin/bash\n\n# Input: the unfixed gtf, and the fixed gtf.\n\n# Possible events:\n#\texon removal - in-frame, simple\n#\texon addition - in-frame, simple\n#\tintron removal - in-frame, simple\n#\tintron addition - in-frame, simple\n#\tmoved start codon - in-frame, simple\n#\tmoved start codon - in-frame, intron-containing\n#\t\tmoved intron start - in-frame, simple\n#\t\tmoved intron end - in-frame, simple\n#\tmoved start codon - frame change,simple\n#\tmoved start codon - frame change, intron-containing.\n#\t\tmoved intron start - frame change, simple\n#\t\tmoved intron end - frame change, simple\n#\texon remains exactly the same\n#\texon addition - frame change, simple\n#\texon removal - frame change, simple\n#\tintron addition - frame change, simple\n#\tintron removal - frame change, simple\n#\t\tcomplex events (all other events).\n#\n# Output is in the format\n#\tevent type...inframe?...exoncontaining...changedbases\n\n#echo "====================="\n#echo $1\n\ngene=$1\nfixed=$2\n\not=`mktemp -d `\ngene_cds="$ot/geneCds"\nfixed_cds="$ot/fixedCds"\n\nawk \'$3 == "CDS"\' $gene | awk -F "\\t" \'BEGIN{OFS="\\t"}{$9="name_here"; print}\' | sort -k4,4n > $gene_cds\nawk \'$3 == "CDS"\' $fixed | awk -F "\\t" \'BEGIN{OFS="\\t"}{$9="name_here"; print}\' | sort -k4,4n > $fixed_cds\n\n# If the fixed gene is empty, just class this as gene removal and continue\n# Get length of original gene.\nif [[ ! -s $fixed_cds ]]; then \n\tremoved=`cat $gene_cds | awk \'BEGIN{a=0} {a=a+$5-$4+1} END{print 0-a}\'`\n\techo -e "removed_gene\\tinframe\\tcomplex\\t-\\t$removed"\n\texit\nfi\n\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $gene_cds > $gene_cds.adj\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_cds > $fixed_cds.adj\n\n# Get the total amount added and subtracted from the original.\n\nremoved=`bedtools subtract -a $gene_cds -b $fixed_cds.adj | awk \'BEGIN{a=0} {a=a+$5-$4+1} END{print 0-a}\'`\nadded=`bedtools subtract -b $gene_cds.adj -a $fixed_cds | awk \'BEGIN{a=0} {a=a+$5-$4+1} END{print a}\'`\necho -e "total_added\\t-\\t-\\t$added"\necho -e "total_removed\\t-\\t-\\t$removed"\n\nstrand=`head -n1 $gene_cds | cut -f7`\n\nif [[ $strand == "-" ]]; then\n\tmaxval=`cat $gene_cds $fixed_cds | cut -f4,5 | sed -r "s/\\t/\\n/g" | sort -n | tail -n1`\n\tawk -v m="$maxval" \'BEGIN{OFS="\\t"} {b=m+1-$4; a=m+1-$5; $4=a; $5=b; print}\' $gene_cds | sort -k4,4n > $gene_cds.tmp\n\tawk -v m="$maxval" \'BEGIN{OFS="\\t"} {b=m+1-$4; a=m+1-$5; $4=a; $5=b; print}\' $fixed_cds | sort -k4,4n > $fixed_cds.tmp\n\tmv $fixed_cds.tmp $fixed_cds\n\tmv $gene_cds.tmp $gene_cds\nfi\n\n# Do the start codon stuff first\n\nstarto=`head -n1 $gene_cds | cut -f 4`\nstartx=`head -n1 $fixed_cds | cut -f 4`\n\n# Check for exons that have remained exactly the same\nbedtools intersect -b $gene_cds  -a $fixed_cds -f 1 -F 1 | sed -r "s/.*/nochange\\tinframe\\tsimple\\t0/g"\nbedtools intersect -b $gene_cds  -a $fixed_cds -f 1 -F 1  > $ot/identical\n\nif [[ $starto -ne $startx ]]; then\n\tstart_chunk="$ot/start_chunk"\n\t# For bedtools and its annoying subtraction problem\n\tif [[ $starto -ge $startx ]]; then\n\t\thead -n1 $gene_cds | awk -v a=$starto -v b=$startx \'BEGIN{OFS="\\t"} {$4 = b; $5 = a -1; print}\' > $start_chunk\n\t\tbedtools intersect -b $start_chunk -a $fixed_cds > $start_chunk.is\n\t\tawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $start_chunk > $start_chunk.adj\n\t\tbedtools subtract -a $fixed_cds -b $start_chunk.adj > $fixed_cds.headless\n\t\tbedtools subtract -a $gene_cds -b $start_chunk.adj > $gene_cds.headless\n\t\tchange=`cat $start_chunk.is | awk \'BEGIN{a=0} {a=a+$5-$4+1} END{print a}\'`\n\t\t# Check frame\n\t\tif [[ $diff_frame -eq 0 ]]; then\n\t\t\tmessage1="moved_start\\tinframe"\n\t\telse\n\t\t\tmessage1="moved_start\\tframeshift"\n\t\tfi\n\t\t# Check whether any introns have been introduced by the new start codon.\n\t\tif cmp -s "$start_chunk.is" "$start_chunk"; then\n\t\t\tmessage2="simple"\n\t\telse\n\t\t\tmessage2="introncontaining"\n\t\tfi\n\t\techo -e "$message1\\t$message2\\t$change"\n\telse\n\t\thead -n1 $gene_cds | awk -v a=$starto -v b=$startx \'BEGIN{OFS="\\t"}{$4 = a; $5 = b -1; print}\' > $start_chunk\n\t\tbedtools intersect -b $start_chunk -a $gene_cds > $start_chunk.is\n\t\tawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $start_chunk > $start_chunk.adj\n\t\tbedtools subtract -a $fixed_cds -b $start_chunk.adj > $fixed_cds.headless\n\t\tbedtools subtract -a $gene_cds -b $start_chunk.adj > $gene_cds.headless\n\t\tchange=`cat $start_chunk.is | awk \'BEGIN{a=0} {a=a+$5-$4+1} END{print a}\'`\n\t\tif [[ $diff_frame -eq 0 ]]; then\n\t\t\tmessage1="moved_start\\tinframe"\n\t\telse\n\t\t\tmessage1="moved_start\\tframeshift"\n\t\tfi\n\t\t# Check whether any introns have been introduced by the new start codon.\n\t\tif cmp -s "$start_chunk.is" "$start_chunk"; then\n\t\t\tmessage2="simple"\n\t\telse\n\t\t\tmessage2="introncontaining"\n\t\tfi\n\t\techo -e "$message1\\t$message2\\t-$change"\n\tfi\n\tgene_cds=$gene_cds.headless\n\tfixed_cds=$fixed_cds.headless\nfi\n\n# Check for simple exon removal events\nloneexons_o="$ot/lone_o"\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_cds > $fixed_cds.adj\nbedtools subtract -A -a $gene_cds -b $fixed_cds.adj > $loneexons_o\nawk \'{a=$5 - $4+ 1; if(a % 3 == 0) {print "removed_exon\\tinframe\\tsimple\\t-"a} else {print "removed_exon\\tframeshift\\tsimple\\t-"a}}\' $loneexons_o\n\n# Check for simple exon addition  events\nloneexons_x="$ot/lone_x"\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $gene_cds > $gene_cds.adj\nbedtools subtract -A -b $gene_cds.adj -a $fixed_cds > $loneexons_x\nawk \'{a=$5 - $4+ 1; if(a % 3 == 0) {print "added_exon\\tinframe\\tsimple\\t"a} else {print "added_exon\\tframeshift\\tsimple\\t"a}}\' $loneexons_x\n\n# Adjust our copy of the x\'d file once we\'ve acknowledged the changes.\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $loneexons_x > $loneexons_x.adj\nbedtools subtract -A -a $fixed_cds -b $loneexons_x.adj > $fixed_cds.tmp\ncat $fixed_cds.tmp $loneexons_o > $fixed_cds\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\'  $ot/identical >  $ot/identical.adj\nbedtools subtract -A -a  $fixed_cds -b $ot/identical.adj > $fixed_cds.tmp\nbedtools subtract -A -a  $gene_cds -b $ot/identical.adj > $gene_cds.tmp\nmv $fixed_cds.tmp $fixed_cds\nmv $gene_cds.tmp $gene_cds\n\n#For intron checking -- invert.\n#cat $fixed_cds\n#cat $gene_cds\n\n# Check for simple intron removal events\nbase_gene="$ot/base_gene"\nbase_fixed="$ot/base_fixed"\n\nstarto=`sort -k4,4n $gene_cds | head -n1 | cut -f 4`\nendo=`sort -k4,4n $gene_cds | tail -n1 | cut -f 5`\nstartx=`sort -k4,4n $fixed_cds | head -n1 | cut -f 4`\nendx=`sort -k4,4n $fixed_cds | tail -n1 | cut -f 5`\n\nhead -n1 $gene_cds | awk -v a=$starto -v b=$endo \'BEGIN{OFS="\\t"} {$4 = a; $5 = b; print}\' > $base_gene\nhead -n1 $fixed_cds | awk -v a=$startx -v b=$endx \'BEGIN{OFS="\\t"} {$4 = a; $5 = b; print}\' > $base_fixed\n\ngene_introns="$ot/gene_introns"\nfixed_introns="$ot/fixed_introns"\n\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $gene_cds > $gene_cds.adj\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_cds > $fixed_cds.adj\n\nbedtools subtract -a $base_gene -b $gene_cds.adj > $gene_introns\nbedtools subtract -a $base_fixed -b $fixed_cds.adj > $fixed_introns\n\n# Check for simple intron addition events\nloneintrons_o="$ot/lonei_o"\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_introns > $fixed_introns.adj\nbedtools subtract -A -a $gene_introns -b $fixed_introns.adj > $loneintrons_o\nawk \'{a=$5 - $4 + 1; if(a % 3 == 0) {print "removed_intron\\tinframe\\tsimple\\t"a} else {print "removed_intron\\tframeshift\\tsimple\\t"a}}\' $loneintrons_o\n\n# Check for simple intron removal events\nloneintrons_x="$ot/lonei_x"\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $gene_introns > $gene_introns.adj\nbedtools subtract -A -b $gene_introns.adj -a $fixed_introns > $loneintrons_x\nawk \'{a=$5 - $4 + 1; if(a % 3 == 0) {print "added_intron\\tinframe\\tsimple\\t-"a} else {print "added_intron\\tframeshift\\tsimple\\t-"a}}\' $loneintrons_x\n\n# Check for complex events and adjust our copy of the x\'d file once we\'ve acknowledged the changes.\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $loneintrons_x > $loneintrons_x.adj\nbedtools subtract -A -a $fixed_introns -b $loneintrons_x.adj > $fixed_introns.tmp\ncat $fixed_introns.tmp $loneintrons_o > $fixed_introns\n\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_introns > $fixed_introns.adj\nbedtools subtract -a  $base_fixed -b $fixed_introns.adj > $fixed_cds\n\nbedtools intersect -a $gene_cds -b $fixed_cds -wa -wb > $ot/cds_intersection\ncut -f1-9 $ot/cds_intersection | sort | uniq -u > $gene_cds.good\ncut -f1-9 $ot/cds_intersection | sort | uniq -d > $gene_cds.junk\ncut -f10-18 $ot/cds_intersection | sort | uniq -u > $fixed_cds.good\ncut -f10-18 $ot/cds_intersection | sort | uniq -d > $fixed_cds.junk\n\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_cds.good > $fixed_cds.good.adj\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $gene_cds.good > $gene_cds.good.adj\nbedtools subtract -b $gene_cds.good.adj -a $fixed_cds.good | awk \'{a=$5 - $4 +1; if(a % 3 == 0) {print "exon_extension\\tinframe\\tsimple\\t"a} else {print "exon_extension\\tframeshift\\tsimple\\t"a}}\'\nbedtools subtract -a $gene_cds.good -b $fixed_cds.good.adj | awk \'{a=$5 - $4 +1; if(a % 3 == 0) {print "exon_contraction\\tinframe\\tsimple\\t-"a} else {print "exon_contraction\\tframeshift\\tsimple\\t-"a}}\'\n\n# Junk\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $gene_cds.junk  > $gene_cds.junk.adj\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_cds.junk > $fixed_cds.junk.adj\njunk_removed=`bedtools subtract -a $gene_cds.junk.adj -b $fixed_cds.junk | awk \'BEGIN{a=0} {a=a+$5-$4 +1} END{print 0-a}\'`\njunk_added=`bedtools subtract -b $gene_cds.junk -a $fixed_cds.junk.adj | awk \'BEGIN{a=0} {a=a+$5-$4 +1} END{print a}\'`\necho -e "other_added\\t-\\t-\\t$junk_added"\necho -e "other_removed\\t-\\t-\\t$junk_removed"\n\n#cat $gene_cds\n#cat $fixed_cds\n\n\nrm -r $ot\n')
		
########################################
# Choice function
########################################

def mostCoherentSet(sequences, path_wDir):
	# We take consensus slices to have something
	# to compare things to.
	sequences = list(sequences)
	sbm = groupSequencesByGeneregion(sequences)
	css, order = getBestSlices(sbm)
	writeSlices(css, order, path_wDir + "/slices.aln")
	# Remember where all the sequences are.
	generegionlookup, allSequences = getGeneregionLookup(order, sbm)
	checkpoints = getCheckpoints(generegionlookup)
	subsets = getSubsetStrings(css, allSequences, generegionlookup)
        winners = randomBestSubsets(subsets, checkpoints)
  	# If there are multiple winning sets, pick the best ones.
        options = {}; bestScore = -1000; bestId = 0
        for winner in winners:
                winningseqs = {}
                for i, w in enumerate(bin(winner).split("b")[1].zfill(len(sequences))):
                        if int(w) == 1:
                                k = order[generegionlookup[i]]
                                winningseqs[k] = winningseqs.get(k,[]) + [allSequences[i]]
                if len(winners) == 1 and all([len(winningseqs[k]) == 1 for k in winningseqs]):
                        return dict((k, winningseqs[k][0]) for k in winningseqs)
                else:
                        possibilities = itertools.product(*winningseqs.values())
                        for p in possibilities:
                                score = alignmentScore(p)
                                optionId = len(options)
                                options[optionId] = {"seqs": dict((order[i], p[i]) for i in range(len(order))), "score": score}
                                if score > bestScore:
                                        bestId = optionId
                                        bestScore = score
        return options[bestId]["seqs"]

def writeSlices(css, order, path_out):
	sequences = dict((k, string.join([j[0][i] for j in css],"")) for i,k in enumerate(order))
	writeSeqs([makeSeq(k, sequences[k]) for k in sequences], path_out)

def randomBestSubsets2(subsets, checkpoints):
        full = "1"*len(subsets[0])
        ss = [int(s,2) for s in subsets if not s == full]
	if not ss: return [int(full,2)]
	checksB = [int(s,2) for s in checkpoints]
        results = []
	sealed = []
        for i in range(1000):
                # Pick a random starting point. This will be naturally
                # weighted by the most abundant entries.
                st = [s for s in ss if not s in sealed]
		if not st: break
                while True:
			sp = [s for s in st if not s in sealed]
			if not sp:
				break
                        r = random.choice(sp)
                        st = [binAnd(r, i) for i in st if binCompat(i, r, checksB)]
			sp = [s for s in st if not s in sealed]
			if not sp:
				sealed += [r]
				break
                        sq = [s for s in st if not s == r]
                        if sq:
                                st = sq[:]
                        else:
                                results += [r]
				sealed += [r]
                                break
	ress = set(results)
	maxsup = 0; winners = []
	for i in ress:
		sup = support(i, ss)
		if sup > maxsup:
			winners = [i]
			maxsup = sup
		elif sup == maxsup:
			winners += [i]
        return winners

def randomBestSubsets(subsets, checkpoints):
        full = "1"*len(subsets[0])
        ss = [int(s,2) for s in subsets if not s == full]
	if not ss: return [int(full,2)]
	checksB = [int(s,2) for s in checkpoints]
        results = []
        for i in range(1000):
                # Pick a random starting point. This will be naturally
                # weighted by the most abundant entries.
                st = ss[:]
                while True:
                        r = random.choice(st)
                        st = [binAnd(r, i) for i in st if binCompat(i, r, checksB)]
                        sq = [s for s in st if not s == r]
                        if sq:
                                st = sq[:]
                        else:
                                results += st
                                break
	ress = set(results)
	maxsup = 0; winners = []
	for i in ress:
		sup = support(i, ss)
		if sup > maxsup:
			winners = [i]
			maxsup = sup
		elif sup == maxsup:
			winners += [i]
        return winners

def getGeneregionLookup(order, sbm):
	allSequences = []
	generegionlookup = {}
	for i, k in enumerate(order):
		for sequence in sbm[k]:
			allSequences.append(sequence)
			generegionlookup[len(allSequences)-1] = i
	return generegionlookup, allSequences

def getSubsetStrings(css, allSequences, generegionlookup):
	subsets = []
	for posi, pos in enumerate(css):
		for option in pos:
			posstring = ""
			for i, sequence in enumerate(allSequences):
				aa = sequence[posi]
				cs_aa = option[generegionlookup[i]]
				posstring += str(int(aa == cs_aa))
			subsets += [posstring]
	return subsets

def getCheckpoints(generegionlookup):
	checkpoints = dict((k,["0"]*len(generegionlookup)) for k in generegionlookup.values())
	for k in generegionlookup:
		checkpoints[generegionlookup[k]][k] = "1"
	checkpoints = [string.join(a, "") for a in checkpoints.values()]
	return checkpoints

def sliceAlignment(sbm):
	l = len(sbm[sbm.keys()[0]][0])
	dicty = dict((k, {}) for k in sbm)
	for i in range(0,l):
		for d in sbm:
			dicty[d][i] = [k[i] for k in sbm[d]]
	return dicty

def groupSequencesByGeneregion(sequences):
	dicty = {}
	for s in sequences:
		generegion = re.sub(r"(generegion_[0-9]*).*", r"\1", s.id)
		dicty[generegion] = dicty.get(generegion, []) + [s]
	return dicty

def getBestSlices(sbm):
	ssm = sliceAlignment(sbm)
	css = []; slen = len(sbm[sbm.keys()[0]][0])
	order = [k for k in ssm]
	for i in range(0,slen):
		sslice = [ssm[k][i] for k in order]
		c = consensus(sslice)
		if not c:
			sprint("this slice doesn't yield>..")
			sprint(sslice)
			sys.exit()
		css.append(consensus(sslice))
	return css, order

def consensus(slices):
	slc = [set(i) for i in slices]
	cp = consensusPatterns(slc)
	cpall = string.join(cp, "")
	slc = [[a for a in x if a in cpall] for x in slc]
	return grabRepresentatives(slc, cp)

def consensusPatterns(slices):
	"""Choose the best subset, one from each slice.
	   Input is a list of lists of AAs.
	"""
	slicetypes = {}
	for s in slices:
		slicetype = sorted(set(s))
		slicestring = string.join(slicetype,"")
		slicetypes[slicestring] = slicetypes.get(slicestring,0)+1
	# Need to check whether these are equivalent for my purposes
	st_choices = dict((s, [slicetypes[s]*i for i in s]) for s in slicetypes)
	# Get all the choices.
	choices = set([string.join(sorted(string.join(c,"")),"") for c in itertools.product(*st_choices.values())])
	scores = dict((c, colScore(c)) for c in choices)
	maxscore = max(scores.values())
	maxscorers = [s for s in scores if scores[s] == maxscore]
	return maxscorers

def getOrAln(dict_generegions, dict_seqInfo, path_orFa, path_orAln):
	writeSeqs(gatherSequences(dict_generegions, dict_seqInfo), path_orFa)
	align(path_orFa, path_orAln)
	
def gatherSequences(dict_generegions, dict_seqInfo):
	for generegion, dg in dict_generegions.items():
		for sequence in dg["sequences"]:
			yield makeSeq("original_"+generegion+"."+sequence, dict_seqInfo[sequence]["sequence_aa"])

def grabRepresentatives(slc, cp):
	reps = Counter(string.join([string.join(s,"") for s in slc],""))
	ret =[]
	for c in cp:
		poss = Counter(c)
		essential = []
		for s in slc:
			res = s
			for x in s:
				if reps[x] == poss[x]: res = [x]
			essential.append(res)
		# I haven't been able to construct an example where the product of /essential is
		# not the set we want. I'm sure they must exist, but for now I'm just going to
		# keep them all (double checking that they match the required pattern), and leave
		# the special cases to their special ways.
		for i in itertools.product(*essential):
			# Note the consensus patterns must be sorted
			if string.join(sorted(i), "") == c:
				ret.append(i)
			else:
				pass
	return ret

########################################
# Part combinations
########################################

def getCombinationsFasta(partslist, combis, cds, generegion, minintron, minexon, checkStart=True):
	labels={}
	if not partslist: return labels
	for j, c in enumerate(combis):
		partstring = ""; so_cds =""
		clist = []; parts = []
		for i in range(0, len(c)/2):
			part  = partslist[i]
			gtf   = part["gtfline"]
			nb    = c[2*i]; ne = c[2*i + 1]
			partc = copy.deepcopy(part)
			partc["gtfline"][3] = nb; partc["gtfline"][4] = ne
			# Have to extract the frame carefully here
			#partc["gtfline"][7] = [a for a in partc["optionsLeft"] if c[2*i] in partc["optionsLeft"][a]][0]
			partc["gtfline"][7] = (3-len(so_cds))%3
			clist += [nb, ne]
			parts.append(partc)
			# Add to the cumulative gene
			so_cds += cds[int(nb)-1: int(ne)]
			partstring += "b"+str(nb)+"e"+str(ne)
		# The beginning might not be a start codon, so the frame isn't necessarily zero.
		# Use the frame of the gtf part to figure out the frame of the variant.
		firstpos   = clist[0]
		firstpart  = partslist[0]
		localframe = [a for a in firstpart["optionsLeft"] if firstpos in firstpart["optionsLeft"][a]][0]
		aa = translatePart(so_cds[localframe:])
		# Do a bunch of checks to make extra sre the gene is well-formed
		# Check that:
		# The gtf parts don't overlap
		# There are no stop codons that aren't at the end of the gene.
		# The thing has a viable start codon (if the first exon becomes very short the start codon can get destroyed).
		startCheck  = (not checkStart) or (isStartCodon(so_cds[0:3]) and localframe == 0)
		containsN   = "n" in so_cds[localframe:].lower()
		noStopCheck = not "*" in aa[:len(aa)-2]
		wellOrdered = all(clist[p] <= clist[p+1] for p in xrange(len(clist)-1))
		goodIntrons = all(clist[(2*p)+1] + minintron <  clist[2*(p+1)] for p in xrange(len(clist)/2 -1))
		# If the gene has made the cut, then happy days: continue.
		if wellOrdered and goodIntrons and noStopCheck and startCheck and not containsN:
			label = generegion + ".superopt_"+str(j)+"."+partstring	
			labels[label] = {"parts": parts, "aa": aa, "seq": makeSeq(label, aa), "generegion": generegion}
	return labels

def getPartCombinations(partslist, cds, checkInitTerm = False, reject = []):
	boundaries = flattenLol([[part["leftFlat"], part["rightFlat"]] for part in partslist])
	boundaries = [[a for a in b if not a in reject] for b in boundaries]
	combis = list(itertools.product(*boundaries))
	combiscopy=combis[:]
	for c in combis:
		for i in range(0, len(c)/2 - 1):
			donor    = partslist[i]
			acceptor = partslist[i + 1]			
			framesd  = [a for a in donor["optionsRight"] if c[2*i + 1] in donor["optionsRight"][a]]
			framesr  = [a for a in acceptor["optionsLeft"] if c[2*i + 2] in acceptor["optionsLeft"][a]]
			# Sometimes we may try to fuse a non-donor or non-acceptor (i.e. an initial or terminal
			# exon) to something else. We can't do this.
			if not (isDonorSite(c[2*i + 1], cds) and isAcceptorSite(c[2*i + 2], cds)):
				if c in combiscopy: combiscopy.remove(c)
			if not any((framed + framer) % 3 == 0 for (framed, framer) in itertools.product(framesd, framesr)):
				if c in combiscopy: combiscopy.remove(c)
	return combiscopy


def fetchByLabel(lbls, wnrs):
	return dict((gr, [] if wnrs[gr].id == gr + ".blank" else lbls[gr][wnrs[gr].id]["parts"]) for gr in wnrs)


########################################
# Basic fixing
########################################

def fixIt(adj, parts, d_gr, path_wDir, minintron, minexon, path_winnersAln):
	"""This function repeatedly adds in small chunks of gene and picks the best combination at
	   each stage. Each stage is wiggled recursively to aim for optimal results.
	"""
	# Run through each set of adjacent gene parts, sequentially adding them in.
	res = idict(d_gr, [])
	for a, ar in enumerate(adj):
		for r in ar:
			generegion, partid = parseId(r)
			if not generegion in res: res[generegion] = []
			latestpart = parts[generegion][partid]
			alreadyThere = [x for x in res[generegion] if int(x["id"]) == int(partid)]
			if alreadyThere:
				alreadyThere[0]["status"] = "fresh"
			if not alreadyThere:
				latestpart["status"] = "fresh"
				# Make sure the previous part is waiting.
				# Can't be terminal
				if res[generegion]:
					res[generegion][-1]["status"] = "waitingA"
					res[generegion][-1]["terminal"] = False
				res[generegion].append(latestpart)
			firstPartId = min([int(x["id"]) for x in res[generegion]])
			for x in res[generegion]:
				if x["id"] == firstPartId:
					x["initial"] = True
		path_fix = makeIfAbsent(path_wDir +"/fix/fix_" + str(a))	
		sprint("Fixing "+str(a))
		res = incrementalFixRecursive(res, path_fix, minintron, minexon, path_winnersAln, d_gr)
	sprint("finishing genes.")
	path_final  = makeIfAbsent(path_wDir +"/fix/fix_final")
	res = incrementalFixRecursive(res, path_final, minintron, minexon, path_winnersAln, d_gr, tTerminal=True)
	return res

def incrementalFix(parts, p_lDir, mi, mx, path_winnersAln, d_gr, refine=False, tTerminal = False):
	"""Takes a set of parts and wiggles the ends around.
	"""
	path_fasta = p_lDir + "/options.all.fasta"
	# Wiggle each generegion separately and then consider together at the end.
	labels, options = preparePartSets(parts, path_winnersAln, tTerminal, mi, mx, p_lDir, d_gr)
	# Finally, compare the options to choose the winners.
	return processLabels(labels, options, p_lDir, path_winnersAln, False, refine)

def incrementalFixRecursive(parts, path_fDir, mi, mx, path_winnersAln, d_gr, tTerminal=False):
	iteration = 0; iterate = True
	inparts = copy.deepcopy(parts)
	results = idict(parts, [])
	stringsets = []; resses = {}; prevAlns = {}
	while iterate:
		path_iDir = makeIfAbsent(path_fDir + "/iteration_" + str(iteration))
		res, prevAln = incrementalFix(inparts, path_iDir, mi, mx, path_winnersAln, d_gr, tTerminal=tTerminal)
		partstrings = dict((a,getPartString(a, res[a])) for a in res)
		for i,e in enumerate(stringsets):
			if all([partstrings[k] == e[k] for k in partstrings]):
				iterate = False
				choiceIndices = range(0, iteration)
		resses[iteration]   = res
		prevAlns[iteration] = prevAln+".re"
		align(prevAln, prevAln + ".re")
		stringsets += [partstrings]
		iteration  +=1
		inparts     = copy.deepcopy(res)
	# Pick whichever res has the best alignment score.
	# If there's only one, just return that.
	winner  = chooseIteration(prevAlns, choiceIndices)
	res     = refineStatuses(resses[winner])
	return res

def statusCheck(data, list_good):
	# Can accept both strings and part lists...
	if type(data) == dict and "status" in data: status = data["status"]
	elif type(data) == str: status = data
	else: raise Exception
	return any(status.startswith(a) for a in list_good)

def killBadParts(res, partsList):
	parts = copy.deepcopy(res)
	for part in parts:
		if not (part["leftFlat"] and part["rightFlat"]):
			partsList = [a for a in partsList if int(a["id"]) != int(part["id"])]
			prevIds = [int(a["id"]) for a in partsList if int(a["id"]) < int(part["id"])]
			if prevIds:
				for opart in partsList:
					if opart["id"] == max(prevIds):
						opart["status"] = "waitingC"
						opart["gtfline"][4] = opart["ogtfline"][4]
	return parts, partsList

def scanLeft(i, partsList, gtf, part, status, cds):
	# If the status is fresh, lock the RHS.
	if statusCheck(status, ["fresh"]): part["optionsRight"][getEndFrame(gtf)] = [int(gtf[4])]
	# If the part is at the beginning, look for a start codon.
	if i == 0: part["optionsLeft"][0] += findStarts(getProtoexon(gtf), cds, True, True, True, False)
	# Otherwise, look for an acceptor
	else: part["optionsLeft"] = findSpliceSites(getProtoexon(gtf), cds, left=True, right=False, happyWithCentral=False)[0]

def scanRight(i, partsList, gtf, part, status, cds, tTerminal):
	# Try to find donors, too.
	if statusCheck(status, ["waiting"]) and not (tTerminal and i == len(partsList) -1):
		part["optionsRight"] = findSpliceSites(getProtoexon(gtf), cds, left=False, right=True, happyWithCentral=False)[1]
	# If this is the last part, also try adding in a terminal option.
	if statusCheck(status, ["waiting"]) and i == len(partsList) -1:
		stops = findStops(getProtoexon(gtf), cds)
		if tTerminal: part["optionsRight"] = idict([0,1,2], [])
		if stops: part["optionsRight"][0] += [stops[0]]
		# Allow for microstops. We can't tell at this point but we'll double check later.
		if gtfLength(gtf) <= 3: part["optionsRight"][0] += [int(gtf[4])+1]

def prepareParts(partsListIn, cds, tTerminal=False):
	partsList = copy.deepcopy(partsListIn)
	finished  = False
	while not finished:
		res = []
		for i, part in enumerate(partsList):
			part   = initialiseOptions(part)
			gtf    = part["gtfline"]
			status = part["status"]
			if not "ogtfline" in part.keys(): part["ogtfline"] = gtf[:]
			if statusCheck(status, ["fresh", "double"]):   scanLeft(i, partsList, gtf, part, status, cds)
			if statusCheck(status, ["waiting", "double"]): scanRight(i, partsList, gtf, part, status, cds, tTerminal)
			# flatten options and add to result.
			res.append(flattenOptions(part))
		# Some parts will prevent the gene from being constructed properly. Check for these.
		# If we can't find boundaries for a part, we have to kill it.
		if all([part["leftFlat"] and part["rightFlat"] for part in res]): finished = True
		else: res, partsList = killBadParts(res, partsList)
	return res

def getReworkedCombinations(plist, cds, generegion, minintron, minexon, tTerminal = False):
	plist = prepareParts(plist, cds, tTerminal)
	return getCombinationsFasta(plist, getPartCombinations(plist, cds), cds, generegion, minintron, minexon) if plist else {}

def noFix(parts):
	pc = parts[:]
	for p in pc: p["status"] = "done"
	return pc

def tinypart(plist):
	# Check if any of the most recent parts are too small to be useful.
	if len(plist) < 2: return False
	lastpartgtf = plist[-2]["gtfline"]
	return (int(lastpartgtf[4]) -int(lastpartgtf[3]) < 40)

def checkIntegrity(plist, cds, gr, mi, mx, labels_l, tryfix, tTerminal):
	if tryfix and ((not labels_l) or tinypart(plist)):
		plist1, plist2 = reworkEnds(plist)
		labels_l1 = getReworkedCombinations(plist1, cds, gr, mi, mx, tTerminal = tTerminal)
		labels_l2 = getReworkedCombinations(plist2, cds, gr, mi, mx, tTerminal = tTerminal)
		for l in labels_l1: labels_l[l] = labels_l1[l]
		for l in labels_l2: labels_l[l] = labels_l2[l]
	return labels_l

def checkStops(plist, cds, gr, mi, mx, labels_l):
	plist1 = copy.deepcopy(plist[:-1])
	plist1[-1]["terminal"] = True
	plist1[-1]["status"] = "waitingQ"
	plist1   = prepareParts(plist1, cds, tTerminal = True)
	labels_t = getReworkedCombinations(plist1, cds, gr, mi, mx, tTerminal = True)
	for l in labels_t: labels_l[l] = labels_t[l]
	return labels_l

def preparePartSets(parts, path_winnersAln, tTerminal, minintron, minexon, path_fDir, d_gr):
	labels = idict(parts, {}); options = idict(parts, [])
	write(parts, path_fDir +"/prevparts.tpy")
	for gr in parts:
		if parts[gr]:
			cds = d_gr[gr]["cdsBase_data"]
			# Prepare the initial round of parts
			tryfix = not d_gr[gr]["keep"]
			plist = prepareParts(parts[gr] if tryfix else noFix(parts[gr]), cds, tTerminal=tTerminal)
			write(plist, path_fDir+"/plist.tpy")
			if plist:
				writeSeqs([makeSeq(str(p["id"]), p["part"]) for p in parts[gr]], path_fDir +"/parts." + gr+".fasta")
				# Now get some part combinations
				labels_l = getCombinationsFasta(plist, getPartCombinations(plist, cds), cds, gr, minintron, minexon)
				# Check if any of the combis are empty. If so, try and rearrange the most recent bits.
				# At best we should be able to get the previous lot back.
				labels_l = checkIntegrity(plist, cds, gr, minintron, minexon, labels_l, tryfix, tTerminal = tTerminal)
				# If this is the terminal step, double check that removing the last part of
				# each gene does not improve things.
				if tTerminal and len(plist) > 1: labels_l = checkStops(plist, cds, gr, minintron, minexon, labels_l)
				# Get the options
				labels[gr], options[gr] = labels_l, [l["seq"] for l in labels_l.values()]
	return labels, options

def refineStatuses(parts):
	partsnew = idict(parts, [])
	for generegion in parts:
		partsDict = dict((x["id"], x) for x in parts[generegion])
		partIds   = [int(x["id"]) for x in parts[generegion]]
		for i, partId in enumerate(partIds):
			part = partsDict[partId]
			part["status"] = "waitingD" if statusCheck(part, ["fresh", "double"]) or (not part["terminal"] and int(partId) == max(partIds)) else "done"
			partsnew[generegion].append(part)
	return partsnew

def setStatusById(plist, pid, status):
	for p in plist:
		if p["id"] == pid: p["status"] = status

def setTerminalById(plist, pid):	
	for p in plist:
		if p["id"] == pid: p["terminal"] = True

def reworkEnds(partslist):
	# Release two variants on the partslist
	# 1. Removing the penultimate part
	# 2. Removing the last part.
	p1 = copy.deepcopy(partslist)
	p2 = copy.deepcopy(partslist)
	freshIds = [a["id"] for a in partslist if statusCheck(a, ["fresh"])]
	# We can find ourselves without any freshids, if the end reworking is taking place after a terminal part has been removed.
	lastId = int(freshIds[0]) if freshIds else max([int(a["id"]) for a in partslist])
	prevIds = [int(a["id"]) for a in partslist if int(a["id"]) < lastId]
	# If there is no previous id, simply delete this gene part and return empty.
	if not prevIds: return [],[]
	# If there is only one previous id, try both omitting the fresh id and omitting the original.
	# In this case we have to set the fresh id to initial.
	prevId = max(prevIds)
	if len(prevIds) == 1:
		p1 = [v for v in p1 if not int(v["id"]) == prevId]
		setStatusById(p1, prevId, "doubleZ")
		p2 = [v for v in p1 if not int(v["id"]) == lastId]
	# If there are multiple previous ids, we must attempt to unravel the pen-penultimate one.
	elif len(prevIds) > 1:
		p1 = [v for v in p1 if not int(v["id"]) == prevId]
		prevprevId = max([int(a["id"]) for a in partslist if int(a["id"]) < prevId])
		setStatusById(p1, prevprevId, "waitingE")
		setStatusById(p1, lastId, "doubleP")
		p2 = [v for v in p1 if not int(v["id"]) == lastId]
	oLastPart = [v for v in partslist if int(v["id"]) == lastId][0]
	if oLastPart["terminal"]:
		if p1: setTerminalById(p1, max([int(a["id"]) for a in p1]))
		if p2: setTerminalById(p2, max([int(a["id"]) for a in p2]))
	return p1, p2

def initialiseOptions(part):
	for bag in ["optionsLeft", "optionsRight"]:
		part[bag] = idict([0,1,2], [])
	gtf = part["gtfline"]
	part["optionsLeft"][int(gtf[7])] = [int(gtf[3])]
	part["optionsRight"][getEndFrame(gtf)] = [int(gtf[4])]
	return part

def flattenOptions(part):
	sl = part["optionsLeft"]; sr = part["optionsRight"]
	part["leftFlat"]  = sl[0] + sl[1] + sl[2]
	part["rightFlat"] = sr[0] + sr[1] + sr[2]
	return part

def chooseIteration(prevAlns, choiceIndices):
	winningScore = -float("inf")
	for i in choiceIndices:
		score = getScore(prevAlns[i])
		if score > winningScore:
			winner, winningScore = i, score
	return winner


########################################
# Filtering
########################################

def filterChanges(adj, cParts, cCoords, path_fltrs, dict_cdsBases, path_aln, minintron, minexon, dict_generegions, cMcount):
	# Sort the genes into piles based on gene region.
	compare = {}
	for gene in cParts:
		generegion = getGr(gene)
		compare[generegion] = compare.get(generegion, []) + [gene]
	# For each gene region, find introns and exons that only exist in the 
	# new gene. These come out as generegion coordinates.
	novEx, novIn, novFl = getNovelRegions(compare, cParts, cCoords, dict_generegions)
	# Grab the intervals
	alnseqs = readSeqs(path_aln, "fasta")
	inEx, inIn, inFl = getInspectionIntervals(novEx, novIn, novFl, cParts, cCoords, alnseqs, compare, path_aln)
	# Filter the intervals
	reject = filterIntervals(inEx, inIn, inFl, cParts, cCoords, alnseqs, compare)
	# Remove anything from the reject pile that was present in any of the original gene.
	reject = allowOriginals(reject, cParts, compare)
	return compareParts(adj, cParts, cMcount, path_fltrs, dict_cdsBases, path_aln, minintron, minexon, dict_generegions,  d_reject = reject)

def getIntervals(nov, cParts, compare, alnlen, cCoords, path_aln, grab):
	featureIntervals = {}
	for gene in nov:
		oGenes = [g for g in compare[getGr(gene)] if not g == gene]
		featureIntervals[gene] = flattenLol(grab(feature, cParts, gene, oGenes, alnlen, cCoords, path_aln) for feature in nov[gene])
	return featureIntervals

def getInspectionIntervals(novEx, novIn, novFl, cParts, cCoords, alnseqs, compare, path_aln):
	callFunction("echo \"\" > " + path_aln +".fil")
	exonIntervals = {}; intronIntervals = {}; flankIntervals = {}
	alnlen = len(alnseqs[0])
	exonIntervals   = getIntervals(novEx, cParts, compare, alnlen, cCoords, path_aln, grabExonIntervals)
	intronIntervals = getIntervals(novIn, cParts, compare, alnlen, cCoords, path_aln, grabIntronIntervals)
	flankIntervals  = getIntervals(novFl, cParts, compare, alnlen, cCoords, path_aln, grabFlankIntervals)
	return exonIntervals, intronIntervals, flankIntervals

def grabExonIntervals(exon, cParts, gene, origGenes, alnlen, cCoords, path_aln):
	# This 100% should exist and should be unique.
	rPartId = [a for a in cParts[gene] if int(cParts[gene][a]["gtfline"][3]) == exon[0] and int(cParts[gene][a]["gtfline"][4]) == exon[1]][0]
	rgtf    = cParts[gene][rPartId]["gtfline"]
	intervalProbes = []
	# Go through each of the original genes, inspecting only the junctions
	# that appear to have changed.
	# LHS first
	rPrevIds = [a for a in cParts[gene] if a < rPartId]
	rInterval = [int(cParts[gene][max(rPrevIds)]["gtfline"][4]) if rPrevIds else 0, int(cParts[gene][rPartId]["gtfline"][3])]
	for o in origGenes:
		oPartIds = [k for k in cParts[o] if gtfLinesOverlap(cParts[o][k]["gtfline"], rgtf)]
		oPrevIds = [k for k in cParts[o] if all(k < l for l in oPartIds)]
		oInterval = [int(cParts[o][max(oPrevIds)]["gtfline"][4]) if oPrevIds else 0, int(cParts[o][min(oPartIds)]["gtfline"][3])] if oPartIds else [] 
		intervals = grabAlnIntervalsEx(oInterval, rInterval, gene, o, cParts, cCoords, path_aln +".fil", alnlen, "l")
		if intervals: intervalProbes.append(grabIntervalProbe(rPrevIds, [rPartId], intervals, gene, o, flavour = "ex", direction = "right"))
	# RHS next.
	rNextIds = [a for a in cParts[gene] if a > rPartId]
	rInterval = [int(cParts[gene][rPartId]["gtfline"][4]), int(cParts[gene][min(rNextIds)]["gtfline"][3]) if rNextIds else -1 ]
	for o in origGenes:
		oPartIds = [k for k in cParts[o] if gtfLinesOverlap(cParts[o][k]["gtfline"], rgtf)]
		oNextIds = [k for k in cParts[o] if all(k > l for l in oPartIds)]
		oInterval = [int(cParts[o][max(oPartIds)]["gtfline"][4]), int(cParts[o][min(oNextIds)]["gtfline"][3]) if oNextIds else -1] if oPartIds else [] 
		intervals = grabAlnIntervalsEx(oInterval, rInterval, gene, o, cParts, cCoords, path_aln +".fil", alnlen, "r")
		if intervals: intervalProbes.append(grabIntervalProbe([rPartId], rNextIds, intervals, gene, o, flavour = "ex", direction = "left"))
	return intervalProbes

def grabFlankIntervals(flank, cParts, gene, origGenes, alnlen, cCoords, path_aln):
	# Each flank represents a single end part.
	# This should 100% exist and should be unique
	cpg = cParts[gene].keys()
	rRight  = [cpg[i] for i in cpg if int(cParts[gene][cpg[i]]["gtfline"][3]) == flank[1]+1]
	rLeft = [cpg[i] for i in cpg if int(cParts[gene][cpg[i]]["gtfline"][4]) == flank[0]-1]
	intervalProbes = []
	if rRight:
		rPartId = rRight[0];
		rgtf = cParts[gene][rPartId]["gtfline"][:]
		rgtf[3] = flank[0]
		rgtf[4] = cParts[gene][rPartId]["gtfline"][3]-1
		rInterval = [rgtf[3],rgtf[4]]
		for o in origGenes:
			oflank = int(cParts[o][min(cParts[o].keys())]["gtfline"][3])-1
			oInterval = [flank[0], oflank]
			intervals = grabAlnIntervalsFl(oInterval, rInterval, gene, o, cParts, cCoords, path_aln +".fil", alnlen, "r")
			if intervals: intervalProbes.append(grabIntervalProbe([], [rPartId], intervals, gene, o, flavour = "fl", direction = "right"))
	elif rLeft:
		rPartId = rLeft[0];
		rgtf = cParts[gene][rPartId]["gtfline"][:]
		rgtf[4] = flank[1]
		rgtf[3] = int(cParts[gene][rPartId]["gtfline"][4])+1
		[rgtf[3],rgtf[4]]
		rInterval = [rgtf[3],rgtf[4]]
		for o in origGenes:
			oflank = int(cParts[o][max(cParts[o].keys())]["gtfline"][4])+1
			oInterval = [oflank, flank[1]]
			intervals = grabAlnIntervalsFl(oInterval, rInterval, gene, o, cParts, cCoords, path_aln +".fil", alnlen, "l")
			if intervals: intervalProbes.append(grabIntervalProbe([rPartId], [], intervals, gene, o, flavour = "fl", direction = "left"))
	return intervalProbes

def grabIntronIntervals(intron, cParts, gene, origGenes, alnlen, cCoords, path_aln):
	# Rather than being a single part here, each intron represents a pair of gene parts.
	# This should 100% exist and should be unique
	cpg = cParts[gene].keys()
	rPartIds = [[cpg[i], cpg[i+1]] for i in cpg[:-1] if int(cParts[gene][cpg[i]]["gtfline"][4]) == intron[0]-1 and int(cParts[gene][cpg[i+1]]["gtfline"][3]) == intron[1]+1][0]
	# Construct an intron gtfline for this intron.
	rgtf = cParts[gene][rPartIds[0]]["gtfline"][:]
	rgtf[3] = cParts[gene][rPartIds[0]]["gtfline"][4]+1
	rgtf[4] = cParts[gene][rPartIds[1]]["gtfline"][3]-1
	rInterval = [rgtf[3],rgtf[4]]
	intervalProbes = []
	# Go through each of the original genes, inspecting only the introns
	# that appear to have changed.
	for o in origGenes:
		introngtfs = getIntrons(cParts[o])
		oIntronIds = [k for k in introngtfs if gtfLinesOverlap(introngtfs[k]["gtfline"], rgtf)]
		# If there are no intron ids here, just carry on.
		# This Ross will be sorted out anyway.
		if not oIntronIds: continue	
		oInterval = [introngtfs[min(oIntronIds)]["gtfline"][3], introngtfs[max(oIntronIds)]["gtfline"][4]]
		intervals = grabAlnIntervalsIn(oInterval, rInterval, gene, o, cParts, cCoords, path_aln +".fil", alnlen, force = len(oIntronIds) > 1)
		if intervals: intervalProbes.append(grabIntervalProbe([rPartIds[0]], [rPartIds[1]], intervals, gene, o, flavour = "in", direction = "both"))
	return intervalProbes

def getFeatures(parts, ubound):
	partIds = sorted(parts.keys(), key = lambda x: int(parts[x]["gtfline"][3]))
	exons = []; introns = []; flanks = []
        for i in range(len(partIds) -1):
                gtf1 = parts[i]["gtfline"]
		gtf2 = parts[i+1]["gtfline"]
		if i == 0: exons += [[int(gtf1[3]), int(gtf1[4])]]
		exons += [[int(gtf2[3]), int(gtf2[4])]]
		introns += [[int(gtf1[4])+1, int(gtf2[3])-1]]
	# Add the flanking regions. These will be processed similarly to introns.
	flanks += [[0,int(parts[partIds[0]]["gtfline"][3])-1]]
	flanks += [[int(parts[partIds[-1]]["gtfline"][4])+1, ubound-1]]
	return exons, introns, flanks

def getIntervalJuncs(intervalGroup, cParts, resgene):
	res = []
	# Interval group should contain either one item or two.
	for i in intervalGroup:
		if "r" in i and i["direction"] in ["right", "both"]: res += [cParts[resgene][i["r"]]["gtfline"][3]]
		if "l" in i and i["direction"] in ["left", "both"]: res += [cParts[resgene][i["l"]]["gtfline"][4]]
	return list(set(res))

def getIntervalGroup(tag1, tag2, dirn, pIds, intervals, a, ivalGroups):
	partners = [i for i in intervals if tag1 in i and pIds.index(i[tag1]) == pIds.index(a[tag2]) + dirn \
		and i["ogene"] == a["ogene"] and not i == a]
	if partners:
		if not anyperm([a, partners[0]], ivalGroups):
			return [[a, partners[0]]]
		elif not [a] in ivalGroups:
			return [[a]]
	return []
	
def getIntervalGroups(cParts, resgene, intervals):
	# They either should all be exonic or all not.
	if any([a["flavour"] == "ex" for a in intervals]):
		intervalGroups = []; partIds = cParts[resgene].keys()
		for a in intervals:
			# grab any adjacent parts, to be considered together
			if "r" in a: intervalGroups += getIntervalGroup("l", "r", -1, partIds, intervals, a, intervalGroups)
			if "l" in a: intervalGroups += getIntervalGroup("r", "l", 1, partIds, intervals, a, intervalGroups)
		return intervalGroups
	else:
		return deDup([[a] for a in intervals])

def allowOriginals(reject, cParts, compare):
	# Removes any coordinates in the reject pile that existed in the original genes.
	# We only want to exclude new stuff that we've found
	res = {}
	for generegion in compare:
		if not generegion in reject: continue
		origGenes = [a for a in compare[generegion] if not "result" in a]
		res[generegion] = flattenLol([a for a in reject[generegion].values()])
		for o in origGenes:
			coords = flattenLol([[int(i["gtfline"][3]), int(i["gtfline"][4]),] for i in cParts[o].values()])
			res[generegion] = [r for r in res[generegion] if not r in coords]
	return res

def aq(alnNc, rgene, ogene, compare, chop):
        for grpair in itertools.combinations([a for a in compare if not a == getGr(rgene)], 2):
		region = chop if not alnNc or len(alnNc[0]) == 0 else alnNc
                oseqs = [a for a in region if (not "result" in a.id and getGr(a.id) in grpair) or a.id == ogene]
                rseqs = [a for a in region if ("result" in a.id and getGr(a.id) in grpair) or a.id == rgene]
                if not rseqs or not len(rseqs[0]): return False
                nonzerocols = [i for i in range(len(rseqs[0])) if not all(k[i] == "-" for k in rseqs)]
		#Judge the changed and unchanged regions separately.
		og = [a for a in alnNc if  a.id == ogene][0]
		rg = [a for a in alnNc if  a.id == rgene][0]
		changed      = [i for i in range(len(og)) if og[i] != rg[i]]
		unchanged    = [i for i in range(len(og)) if og[i] == rg[i]]
		changed_nz   = [i for i in changed if not all(k[i] == "-" for k in rseqs)]
		unchanged_nz = [i for i in unchanged if not all(k[i] == "-" for k in rseqs)]
		# The flanking alignment needs to be good as a start.
		if unchanged_nz:
			if not alignmentScore(chopAlignment(alnNc, unchanged_nz))/len(unchanged_nz) > 3: continue
		# If the flanking alignment is okay, then we just need to check that the insides are okay also.
		# This is quite a strict criterion.
		if not changed_nz: return False
		if changed_nz:
			if alignmentScore(chopAlignment(alnNc, changed_nz))/len(changed_nz) > 3: return False
	return True

def minimalGains(alnNc, rgene, ogene, compare):
	if not alnNc: return True
	for grpair in itertools.combinations([a for a in compare if not a == getGr(rgene)], 2):	
		# Approach is as follows. Take clumps of three genes. Get the nonconstant portions
		# of their alignment in the region of interest.
		ogenes = [a for a in alnNc if (not "result" in a.id and getGr(a.id) in grpair) or a.id == ogene]
		rgenes = [a for a in alnNc if ("result" in a.id and getGr(a.id) in grpair) or a.id == rgene]
		og = [a for a in alnNc if  a.id == ogene][0]
		rg = [a for a in alnNc if  a.id == rgene][0]
		# Grab the scores
		oscore = alignmentScore([a for a in alnNc if (not "result" in a.id and getGr(a.id) in grpair) or a.id == ogene])
		rscore = alignmentScore([a for a in alnNc if ("result" in a.id and getGr(a.id) in grpair) or a.id == rgene])
		# If contiguous in both cases and less than or equal to 3 aa change, don't allow.
		# If contiguous and longer, require a total increase of 4 or more
		if not any((og[i] == "-" and rg[i] == "-") for i in range(len(rg))):
			if len(alnNc[0]) < 4: return True
			if rscore > oscore + 4: return False
		# If noncontiguous, require a total increase of 2 or more
		else:
			if rscore > oscore + 2: return False
	return True

def parallels(alnNc, rgene, ogene, compare):
	if not alnNc or len(alnNc[0]) == 0 : return False
	for grpair in itertools.combinations([a for a in compare if not a == getGr(rgene)], 2):
		oseqs = [a for a in alnNc if (not "result" in a.id and getGr(a.id) in grpair) or a.id == ogene]
		rseqs = [a for a in alnNc if ("result" in a.id and getGr(a.id) in grpair) or a.id == rgene]
		# Check to see if there are any columns that are all empty in the original and not in the new.
		if not oseqs: continue
		for i in range(len(oseqs[0])):	
			if all(a[i] == "-" for a in oseqs) and any(a[i] != "-" for a in rseqs): return True
	return False

def getNovelRegions(compare, cParts, cCoords, dict_generegions):
	# Extract all the exons and introns of the result gene as regions.
	# We look in the original gene set to see if there is anything that is present
	# in the original bunch that is not in the new gene.
	novEx = {}; novIn = {}; novFl = {}
	for generegion in compare:
		results = [a for a in compare[generegion] if "result" in a]
		if not results: continue
		# There should be at most one result.
		result    = results[0]
		# Grab the introns and exons from the old and new.
		ubound = dict_generegions[generegion]["baselen"]
		rExons, rIntrons, rFlanks = getFeatures(cParts[result], ubound)
		origs     = [a for a in compare[generegion] if not a == result]
		oFeatures = [getFeatures(cParts[o], ubound) for o in origs]
		oExons    = deDup(flattenLol(([o[0] for o in oFeatures])))
		oIntrons  = deDup(flattenLol(([o[1] for o in oFeatures])))
		oFlanks   = deDup(flattenLol(([o[2] for o in oFeatures])))
		# Check which of these are novel.
		novEx[result] = [e for e in rExons if not e in oExons]
		novIn[result] = [i for i in rIntrons if not i in oIntrons]
		novFl[result] = [i for i in rFlanks if not i in oFlanks]
	return novEx, novIn, novFl

def filterIntervals(inEx, inIn, inFl, cParts, cCoords, alnseqs, compare):
	# Filter the intervals based on various criteria.
	reject = {}
	dontreject = {}
	for inFeature in [inEx, inIn, inFl]:
		for resgene in inFeature:	
			# If any the features are exonic, we try to group them into
			# pairs that we can compare together. Otherwise consider them individually.
			intervals = inFeature[resgene]
			if not intervals: continue
			intervalGroups = getIntervalGroups(cParts, resgene, intervals)
			# Now consider each of the interval groups. There should be
			# at most two items in each interval group. When there are two 
			# items in an interval group, there will be two parts to the junction, so
			# consider them both.
			for i in intervalGroups:
				ival = list(set(flattenLol([interv["interval"] for interv in i])))
				rjuncs = getIntervalJuncs(i, cParts, resgene)
				chop = chopAlignment(alnseqs, ival)
				generegion = re.sub(r"(generegion_[0-9]*).*", r"\1", resgene)
				# Make containers for the results of the different tests.
				if not generegion in reject: reject[generegion] = {}
				if not generegion in dontreject: dontreject[generegion] = {}
				ogene = i[0]["ogene"]
				###############################
				# Alignment improvement check.
				###############################
				# Grab the alignment scores for the original alignment and the new
				seqo = [a for a in alnseqs if a.id == ogene][0]
				seqr = [a for a in alnseqs if a.id == resgene][0]
				nonconstantregion  = [a for a in ival if 0 <= a < len(seqo) and seqo[a] != seqr[a]]
				nonconstantslopped = flattenLol([range(a-10,a+10) for a in nonconstantregion])
				alnNcSlop = chopAlignment(alnseqs, sorted(set(nonconstantslopped)))
				oldscore = alignmentScore([a for a in alnNcSlop if not "result" in a.id and grPick(a.id, ogene)])
				newscore = alignmentScore([a for a in alnNcSlop if "result" in a.id])
				# Check whether the new score is better than the old score.
				# Form specific rejections groups for specific tests.	
				if not "alnimprove" in reject[generegion]: reject[generegion]["alnimprove"] = []
				if not "alnimprove" in dontreject[generegion]: dontreject[generegion]["alnimprove"] = []
				if not oldscore <= newscore:
					reject[generegion]["alnimprove"] += rjuncs
				else:
					dontreject[generegion]["alnimprove"] += rjuncs
				################################
				# Alignment quality check.
				################################
				# Check the alignment quality of the region of interest. If the
				# alignment isn't very good to start with, or if the resulting alignment
				# isn't very good, reject.
				if not "aq" in reject[generegion]: reject[generegion]["aq"] = []
				if not "aq" in dontreject[generegion]: dontreject[generegion]["aq"] = []
				if aq(alnNcSlop, resgene, ogene, compare, chop):
					reject[generegion]["aq"] += rjuncs
				else:
					dontreject[generegion]["aq"] += rjuncs
                                ###############################
                                # Minimal gains check.
                                ###############################
				# Minimal gains: nonconstant portion must exhibit an improvement of >5
				# in one of the triplet scores.
				if not "minimalgains" in reject[generegion]: reject[generegion]["minimalgains"] = []
				if not "minimalgains" in dontreject[generegion]: dontreject[generegion]["minimalgains"] = []
				chopo = [a for a in chop if a.id == ogene][0]
				chopr = [a for a in chop if a.id == resgene][0]
				nonconstantregion = [a for a in range(len(chopo)) if chopo[a] != chopr[a]]
				alnNc = chopAlignment(chop, nonconstantregion)
				if minimalGains(alnNc, resgene, ogene, compare):
					reject[generegion]["minimalgains"] += rjuncs
				else:
					dontreject[generegion]["minimalgains"] += rjuncs
				################################
				# Parallel junctions check
				################################
				# Check for regions where a triplet was blank before, but
				# has become at least partially not blank, including GOI.
				if not "parallels" in reject[generegion]: reject[generegion]["parallels"] = []
				if not "parallels" in dontreject[generegion]: dontreject[generegion]["parallels"] = []
				if parallels(alnNc, resgene, ogene, compare):
					reject[generegion]["parallels"] += rjuncs
				else:
					dontreject[generegion]["parallels"] += rjuncs
#qq				if generegion =="generegion_3":
#qq					print reject[generegion]
#qq	                                print dontreject[generegion]
	for generegion in reject:
		for test in reject[generegion]:
			reject[generegion][test] = list(set([a for a in reject[generegion][test] if not a in dontreject[generegion][test]]))
	return reject

def grPick(agene, ogene):
	return agene == ogene or not getGr(agene) == getGr(ogene)

def getIntrons(pdict):
	res = {}
	keys = sorted(pdict.keys())
	for i in range(len(pdict)-1):
		gtf    = pdict[i]["gtfline"][:]
		gtf[3] = int(pdict[i]["gtfline"][4])+1
		gtf[4] = int(pdict[i+1]["gtfline"][3])-1
		res[i] = {"gtfline": gtf, "ids":[keys[i], keys[i+1]]}
	return res

def grabIntervalProbe(leftright, rightleft, intervals, gene, o, flavour = "", direction = ""):
	res = {"interval": intervals, "gene": gene, "ogene": o}
	if flavour: res["flavour"] = flavour
	if leftright: res["l"] = max(leftright)
	if rightleft: res["r"] = min(rightleft)
	res["direction"] = direction
	return res

def grabAlnIntervalsIn(oInterval, rInterval, gene, o, cParts, cCoords, path_fil, alnlen, force = False):
	iLocs = []
	if not oInterval == rInterval or force:
		cpg = cParts[gene].keys()
		opg = cParts[o].keys()
	        rPartIds = [[cpg[i], cpg[i+1]] for i in cpg[:-1] if int(cParts[gene][cpg[i]]["gtfline"][4]) == rInterval[0]-1 \
					and int(cParts[gene][cpg[i+1]]["gtfline"][3]) == rInterval[1]+1][0]
		riLoc = cCoords[gene][rPartIds[0]][-15:-1] +  cCoords[gene][rPartIds[1]][0:15]
		if oInterval:
			oPartLid = [i for i in opg[:-1] if int(cParts[o][i]["gtfline"][4]) == oInterval[0]-1][0]
			oPartRid = [i for i in opg[1:] if int(cParts[o][i]["gtfline"][3]) == oInterval[1]+1][0]
			oiLoc = cCoords[o][oPartLid][-15:-1] + cCoords[o][oPartRid][0:15]
			iLocs = interpolate(riLoc + oiLoc)
		else:
			iLocs=interpolate(riLoc)
		iLocs = cdsToAaLocs(iLocs)
		filstring = string.join(["Y" if i in iLocs else "-" for i in range(alnlen)], "")
		callFunction("echo \">"+gene+"."+o+"\n"+filstring+"\">> " + path_fil)
	return iLocs

def grabAlnIntervalsFl(oInterval, rInterval, gene, o, cParts, cCoords, path_fil, alnlen, dirn, force = False):
	iLocs = []
	ik1, ik2  = [0,15] if dirn == "r" else [-15,-1]
	if not oInterval == rInterval or force:
		cpg = cParts[gene].keys()
		opg = cParts[o].keys()
		if dirn == "l":
			rPartId = [i for i in cpg if int(cParts[gene][i]["gtfline"][4]) == rInterval[0]-1][0]
			oPartId = [i for i in opg if int(cParts[o][i]["gtfline"][4]) == oInterval[0]-1][0]
			riLoc = cCoords[gene][rPartId][ik1:ik2]
			oiLoc = cCoords[o][oPartId][ik1:ik2]
			iLocs = interpolate(riLoc + oiLoc)
		else:
			rPartId = [i for i in cpg if int(cParts[gene][i]["gtfline"][3]) == rInterval[-1]+1][0]
			oPartId = [i for i in opg if int(cParts[o][i]["gtfline"][3]) == oInterval[-1]+1][0]
			riLoc = cCoords[gene][rPartId][ik1:ik2]
			oiLoc = cCoords[o][oPartId][ik1:ik2]
			iLocs = interpolate(riLoc + oiLoc)
		iLocs = cdsToAaLocs(iLocs)
		filstring = string.join(["T" if i in iLocs else "-" for i in range(alnlen)], "")
		callFunction("echo \">"+gene+"."+o+"\n"+filstring+"\">> " + path_fil)
	return iLocs
	
def grabAlnIntervalsEx(oInterval, rInterval, gene, o, cParts, cCoords, path_fil, alnlen, direction, force = False):
	switchKey = 1 if direction == "l" else 0
	gtfKey    = 3 if direction == "l" else 4
	ik1, ik2  = [0,15] if direction == "l" else [-15,-1]
	cipher    = "D" if direction == "l" else "R"
	iLocs     = []
	if not oInterval == rInterval or force:
		rpart = [a for a in cParts[gene] if int(cParts[gene][a]["gtfline"][gtfKey]) == rInterval[switchKey]][0]
		riLoc = cCoords[gene][rpart][ik1:ik2]
		oiLoc = []
		if oInterval:
			opart = [a for a in cParts[o] if int(cParts[o][a]["gtfline"][gtfKey]) == oInterval[switchKey]][0]
			oiLoc = cCoords[o][opart][ik1:ik2]
			iLocs = interpolate(oiLoc + riLoc)
		else:
			iLocs = riLoc
		iLocs = cdsToAaLocs(iLocs)
		filstring = string.join([cipher if i in iLocs else "-" for i in range(alnlen)], "")
		callFunction("echo \">"+gene+"."+o+"\n"+filstring+"\">> " + path_fil)
	return iLocs

	
########################################
# All the other functions
########################################

def prepareWinnerAlignment(path_wDir, winners, dict_generegions, orFa = ""):
	path_winnersOut = blankFile(path_wDir + "/all.options.winners")#qe
	path_winnersAln = path_winnersOut + ".aln"#qe
	for i in winners:
		generegion, option = parseId(i)
		path_fasta = dict_generegions[generegion]["options"][int(option)]["fasta"]
		callFunction("cat " + path_fasta + " | sed -r \"s/\*/X/g\" >> " + path_winnersOut)
	if orFa: 
		callFunction("cat " + orFa + " >> " + path_winnersOut)
	align(path_winnersOut, path_winnersAln)
	return path_winnersAln

def boundAlignManual(partsInfo, cds, path_stub, slookup_rev):
	"""Align a set of boundaries from the same region.
	"""
	allGtf = flattenLol([[q["gtfline"] for q in partsInfo[p]["parts"]] for p in partsInfo])
	writeCsv(allGtf, path_stub + ".all.gtf")
	mergeFriendly(path_stub + ".all.gtf", path_stub + ".regions.gtf", sort = True)
	regions = dict(enumerate(readCsv(path_stub + ".regions.gtf")))
	# Get a cds alignment, to convert later to an aa alignment.
	# Do this separately for each part in the merge, so that non-compatible parts
	# are not placed on top of each other.
	regionsblank = dict((r, "-"*(int(regions[r][4]) - int(regions[r][3]) + 1)) for r in regions)
	partsout = []; partstrings={}; aaout = []
	for p in partsInfo:
		# Each part will be contained in exactly one of the regions.
		regionlines=copy.deepcopy(regionsblank)
		regionsstring = ""
		for part in partsInfo[p]["parts"]:
			regionposs = [r for r in regions if overlapInFrame(regions[r], part["gtfline"])]
			if regionposs:
				region = regionposs[0]
				line = "-"*(int(part["gtfline"][3]) -int(regions[region][3])) \
						+ cds[int(part["gtfline"][3])-1: int(part["gtfline"][4])] \
						+ "-"*(-int(part["gtfline"][4]) +int(regions[region][4]))
				regionlines[region] = addRegionLine(regionlines[region], line)
		for region in sorted(regions.keys()):
			regionsstring += regionlines[region]
		partstrings[p] = regionsstring
		partsout += [[">"+slookup_rev[p]], [regionsstring]]
	# The method assumes that each gtf is properly written (and starts in-frame), as it should be.
	# Might need to keep an eye on this though...
	aastrings={}
	for p in partstrings:
		aalength = len(partstrings)/3
		aaRaw = partsInfo[p]["aa"]
		aastring = ""
		nucpos = 0
		for i, e in enumerate(partstrings[p]):
			if i % 3 == 0:
				if e == "-":
					aastring += "-"
				else:
					# Ends aren't always full.
					if nucpos / 3 < len(aaRaw):
						aastring += aaRaw[nucpos/3]
			if e != "-": nucpos = nucpos + 1
		aastrings[p] = aastring
	aalength=max([len(aastrings[aa]) for aa in aastrings])
	# There can sometimes be slight length imbalances, due to hanging codons. Fill them in with blanks.
	for aa in aastrings:
		aastring = aastrings[aa]
		if not len(aastring) == aalength:
			aastring = aastring + "-"*(aalength-len(aastring))
		aaout += [[">"+slookup_rev[aa]], [aastring]]
	writeCsv(partsout, path_stub+".regions")
	writeCsv(aaout, path_stub+".regions.aa")
	return path_stub+".regions.aa"

def addRegionLine(existing, new):
	# Add a region over the existing one, such that adding "-" counts for nothing.
	if len(existing) != len(new): raise Exception()
	return string.join([existing[i] if new[i] == "-" else new[i] for i in range(len(existing))].join(""))

def prepareChooseWinners(dict_options, path_wDir):
	path_allFasta = path_wDir + "/all.options.fasta"
	path_allAln   = path_allFasta + ".aln"
	# Sort out the sequences for processing
	seqlookup={}; seqlookup_rev = {}; towrite={}
	opts_copy = copy.deepcopy(dict_options)
	for generegion in opts_copy:
		seqlookup[generegion]={}; towrite[generegion] = []; seqlookup_rev[generegion]={}
		for option in opts_copy[generegion]:
			seqid   = option.id
			shortid = generegion + ".option_" + str(len(seqlookup[generegion]))
			seqlookup[generegion][shortid]   = seqid
			seqlookup_rev[generegion][seqid] = shortid
			option.id = shortid
			towrite[generegion].append(option)
	return path_allFasta, path_allAln, seqlookup, seqlookup_rev, towrite

def chooseWinnersVanilla(dict_options, path_wDir, tryBlank = {}):
	path_allFasta, path_allAln, seqlookup, seqlookup_rev, towrite = prepareChooseWinners(dict_options, path_wDir)
	if not tryBlank and all(len(a) <= 1 for a in dict_options.values()):
		return dict((k, dict_options[k][0]) for k in dict_options if dict_options[k])
	with open(path_allFasta, "w") as f:
		for generegion in towrite: writeSeqs(towrite[generegion], f)
	align(path_allFasta, path_allAln, safety=False)
	return chooseThem(path_allAln, tryBlank, seqlookup, seqlookup_rev, path_wDir)

def chooseWinnersRef(dict_options, path_wDir, tryBlank = {}, path_refAln = "", doubleCheck = False, spanningBand = False, orAln = False, parallelCheck = False):
	path_allFasta, path_allAln, seqlookup, seqlookup_rev, towrite = prepareChooseWinners(dict_options, path_wDir)
	if not tryBlank and all(len(a) <= 1 for a in dict_options.values()):
		return dict((k, dict_options[k][0]) for k in dict_options if dict_options[k])
	# initial microexons can cause zero-length aa strings to be formed, which can't be aligned.
	nonempty = [a for a in flattenLol(towrite.values()) if len(a) !=0]
	empty    = [a for a in flattenLol(towrite.values()) if len(a) ==0]
	writeSeqs(nonempty, path_allFasta)
	# Store the realigned reference so that we can compare the original.
	path_refOut = path_allAln + ".ref"
	alignRef(path_allFasta, path_allAln, path_refAln, path_refOut)
	# Add blank sequences to the alignment file.
	a = readSeqs(path_allAln)
	writeSeqs(a + [makeSeq(t.id, "-"*len(a[0])) for t in empty], path_allAln)
	# If doubleCheck is on, we extract the nonconstant portion from the alignment and
	# align that using l-ins-i (it'll be a relatively narrow column so won't eat up
	# too much compute power). We then compare the alignment score for the two and 
	# go with which ever one is better.
	if doubleCheck: path_allAln, path_refOut = realignBand(path_allAln, path_refOut, spanningBand = spanningBand)
	# If this option is turned on, check for "parallel junction expansion": that is,
	# events that force gaps to be created in the original alignment. Simply reject any such events.
#	if orAln and parallelCheck: path_allAln = rejectParallels(path_allAln, path_refOut)
	return chooseThem(path_allAln, tryBlank, seqlookup, seqlookup_rev, path_wDir)

def rejectParallels(path_allAln, path_refAln):
	# Reject anything which has something when the original alignment has all gaps.
	# The ref and all alignments should be the same length.
	resSeqs = readSeqs(path_allAln)
	refSeqs = [a for a in readSeqs(path_refAln) if "orig" in a.id]
	sbm = groupSequencesByGeneregion(resSeqs)
	path_noparallels = path_allAln + ".nopara"
	toKeep = getNc(sbm)
	if not resSeqs or not refSeqs:
		writeSeqs([], path_noparallels)
		return path_noparallels
	if not len(resSeqs[0]) == len(refSeqs[0]): raise
	reject = []
	resSeqs = chopAlignment(resSeqs, toKeep)
	refSeqs = chopAlignment(refSeqs, toKeep)
	for i in range(len(resSeqs[0])):
		if all(a[i]=="-" for a in refSeqs):
			for a in resSeqs:
				if not a[i] == "-": reject.append(a.id)
	writeSeqs([a for a in resSeqs if not a.id in reject], path_noparallels)
	writeSeqs(resSeqs + refSeqs, path_noparallels+".check")
	return path_noparallels	

def realignBand(path_inAln, path_refAln, thebuffer = 20, spanningBand = False):
	# Find the window in the alignment where the interesting stuff is happening.
	# Then realign this with l-ins-i, and double check the results. Since the
	# aligner may have made a mistake (particularly if the alignment is messy),
	# there might be more than one such region.
	path_column    = path_inAln + ".column.fasta"
	path_columnref = path_inAln + ".column.fasta.ref"
	# Don't go any further if there are no input sequences
	inseqs = readSeqs(path_inAln)
	if not inseqs: return path_inAln, path_refAln
	# Remove any all-blanks...
	seqlen           = len(inseqs[0])
	removecols       = [i for i in range(seqlen) if all(s[i] == "-" for s in inseqs)]
	inseqs_nonblank  = chopAlignment(inseqs, removecols, True)
	refseqs_nonblank = chopAlignment(readSeqs(path_refAln), removecols, True)
	# Get to work. Want to chop out the windows where everything is fine: these 
	# are not interesting to us.
	sbm = groupSequencesByGeneregion(inseqs_nonblank)
	if len(sbm.keys()) < 2: return path_inAln, path_refAln
	# The squeezemap is a binary representation of whether a column should be considered.
	# Grab the coordinates of the nonuniform bits, with a buffer added.
	toKeep = getNc(sbm, thebuffer)
	# There's a choice between taking the span of all of the keep regions and
	# Considering them separately. The spanning region is more accurate but 
	# Could be costly if the regions are distant
	if not toKeep: return path_inAln, path_refAln
	toKeep = range(min(toKeep), max(toKeep)+1) if spanningBand else sorted(set(toKeep))
	sbmalt = dict((k, chopAlignment(sbm[k], toKeep)) for k in sbm)
	refseqs_keep = chopAlignment(refseqs_nonblank, toKeep)
	# Write out the newly found columns.
	for s in refseqs_keep: s.id = "reference."+s.id
	writeSeqs(flattenLol(sbmalt.values()), path_column)
	writeSeqs(refseqs_keep, path_columnref)
	# Retain the lateral context of the band that we've taken out.
	# This involves just adding in the bits we ignored.
	squeezeFa = grabSqueezeMap(inseqs, removecols, toKeep, sbm)
	# Align and compare with the original alignment.
	path_columnAll = path_column + ".all"
	path_columnAln = path_column + ".aln"
	writeSeqs(flattenLol(sbmalt.values()) + refseqs_keep, path_columnAll)
	writeSeqs(squeezeFa, path_column+".squeeze.fa")
	# Align the bastard.
	align(path_columnAll, path_columnAln)
	oldscore = alignmentScore([a for a in readSeqs(path_column) if not "ref" in a.id])
	newscore = alignmentScore([a for a in readSeqs(path_columnAln) if not "ref" in a.id])
	# If the realigned sequences are better, use those. If not use the original ref alignment.
	path_sequences = path_column if oldscore > newscore else path_columnAln
	writeSeqs([a for a in readSeqs(path_sequences) if not "ref" in a.id], path_column+".bestaln")
	writeSeqs([a for a in readSeqs(path_sequences) if "ref" in a.id], path_column+".bestaln.ref")
	return path_column+".bestaln", path_column+".bestaln.ref"

def getNc(sbm, thebuffer = 0):
	seqlen     = len(sbm[sbm.keys()[0]][0])
	squeezeMap = [0 if all(len(set([a[i] for a in sbm[k]])) == 1 for k in sbm) else 1 for i in range(seqlen)]
	toKeep     = flattenLol(range(i-thebuffer, i+thebuffer + 1) for (i,e) in enumerate(squeezeMap) if e == 1)
	return [i for i in toKeep if i >= 0 and i < seqlen]

def grabSqueezeMap(oseqs, removecols, toKeep, sbm):
	squeezeFa = []
	sbmo = groupSequencesByGeneregion(oseqs)
	for k in sbm:
		if not sbmo[k]: continue
		# The sequences by definition are identical when chopped, so just take the first.
		s = sbmo[k][0]
		s_adj  = [s[i] for i in range(len(s)) if not i in removecols]
		s_adjt = [s[i] for i in range(len(s)) if i in removecols]
		ts = string.join([("@" if i in toKeep else s_adj[i]) for i in range(len(s_adj))] + s_adjt, "")
		squeezeFa.append(replaceSeq(s, Seq(ts)))
	return squeezeFa

def chooseWinnersSeeded(dict_options, path_wDir, seed = "", seedInfo = {}, cdsBases = {}, tryBlank = {}):
	path_allFasta, path_allAln, seqlookup, seqlookup_rev, towrite = prepareChooseWinners(dict_options, path_wDir)
	if not tryBlank and all(len(a) <= 1 for a in dict_options.values()):
		return dict((k, dict_options[k][0]) for k in dict_options if dict_options[k])
	if seed == "auto":
		seeds = []
		for generegion in towrite:
			path_dFa = path_allFasta + "."+generegion+".fa"
			writeSeqs(towrite[generegion], path_dFa)
			seeds.append(path_dFa)
		alignSeeded(seeds, path_allAln, prealigned = False, safety = False)
	elif seed == "manual":
		# Manually align the seed sequences. This will only work for sequences from 
		# the same gene region.
		seeds = []
		for generegion in seedInfo:
			if seedInfo[generegion]:
				seeds.append(boundAlignManual(seedInfo[generegion], cdsBases[generegion], path_allFasta + "."+generegion, seqlookup_rev[generegion]))
		alignSeeded(seeds, path_allAln, prealigned = True, safety = True)
	return chooseThem(path_allAln, tryBlank, seqlookup, seqlookup_rev, path_wDir)

def chooseThem(path_allAln, tryBlank, seqlookup, seqlookup_rev, path_wDir):
	# Gather the sequences and pick the best set.
	sequences = [a for a in readSeqs(path_allAln) if not (a.id).startswith("prev.")]
	if not sequences: return {}
	# The propercounts option ensures we fill in absent options with blanks.
	seqlen = len(sequences[0].seq)
	for k in [a for a in tryBlank if tryBlank[a]]:
		sequences.append(makeSeq(k + ".blank", "-"*seqlen))
		if not k in seqlookup:
			seqlookup[k] = {}; seqlookup_rev[k] = {}
		seqlookup[k][k + ".blank"] = k + ".blank"
		seqlookup_rev[k][k + ".blank"] = k + ".blank"
	writeSeqs(sequences, path_wDir+"/options.aln")
	sbm = groupSequencesByGeneregion(sequences)
        if len(sbm.keys()) == 1: return {sbm.keys()[0]: makeSeq(sbm.keys()[0]+".blank", "")}
	winner = mostCoherentSet(sequences, path_wDir); res = {}
	for k, w in winner.items():
		w.id   = seqlookup[k][winner[k].id]
		res[k] = w
	return res

def getWinnerAlnCoords(winners, path_winnersAln, dict_generegions, path_out):
	aligned = [a for a in readSeqs(path_winnersAln) if not "orig" in a.id]
	orders  = {}; 	gtfs    = {}
	for w in winners:
		generegion, option = parseId(w)
		winnerpick = dict_generegions[generegion]["options"][int(option)]
		dict_generegions[generegion]["winner"] = winnerpick
		wparts     = winnerpick["parts"]
		partIds    = sorted(wparts.keys())
		dict_generegions[generegion]["winner"] = wparts
		gtfd  = dict((a, wparts[a]["gtfline"]) for a in partIds)
		# By design this gtf should be clean and safe.
		order = gtfd.keys()
		orders[w] = order
		gtfs[w]   = [gtfd[o] for o in order]
		for key in partIds:
			parth = wparts[key]
			parth["initial"] = (int(key) == 0)
			parth["terminal"] = (int(key) == len(partIds) - 1)
			parth["id"] = int(key)
	locsAa, locs, goodGtfs = flattenAlignment(aligned, gtfs, preSorted = True, path_out = path_out)
	dict_partCoords = {}
	for w in winners:
		generegion, option = parseId(w)
		dict_partCoords[generegion] = {}
		for i,o in enumerate(orders[w]):
			l = locs[w][i]
			dict_partCoords[generegion][o] = [min(l), max(l)]
	return dict_partCoords

def prepareOptions(dict_generegions):
	options = idict(dict_generegions.keys(), [])
	for generegion in dict_generegions:
		dict_options = dict_generegions[generegion]["options"]
		options[generegion] = flattenLol([readSeqs(option["fasta"]) for option in dict_options.values()])
	winners = chooseWinnersVanilla(options, path_wDir).values()
	return [a.id for a in winners]

def prepareCompare(res, dict_generegions, dict_seqInfo, path_dbl):
	# Extract all the information from a res necessary to compare
	# Its parts with the original. We're going to do the comparison
	# A couple of times so it's good to have this separate.
	comparatorGtf         = {}
	comparatorAa	      = {}
	comparatorParts       = {}
	comparatorGeneregions = {}
	comparatorSequences   = {}
	comparatorMcount      = {}
	# Work through all the gene regions and add in the genes.
	for generegion in dict_generegions:
		cds = dict_generegions[generegion]["cdsBase_data"]
		comparatorMcount[generegion] = 1
		# Do all the original sequences first.
		for sequence in dict_generegions[generegion]["sequences"]:
			seqid = generegion+"-original-"+sequence
			# Our result shouldn't contain non-viable genes. HOWEVER, the original input
			# could well do. We must protect ourselves against this.
			# In the case that there is an invalid gene input to begin with and no other result
			# is acquired, output for said gene will simply be blank.
			local_cds = dict_seqInfo[sequence]["sequence_cds"]
			aa        = dict_seqInfo[sequence]["sequence_aa"]
			if containsStop(local_cds[:-3]): continue
			comparatorGtf[seqid]	   = readCsv(dict_seqInfo[sequence]["gtfRel"])
			# Make sure it ends properly!
			comparatorAa[seqid] 	     = re.sub(r"([A-Z])$", r"\1*", aa)
			comparatorGeneregions[seqid] = generegion
			comparatorSequences[seqid]   = sequence
			comparatorMcount[generegion] += 1
		# Finally add in the result. Careful though, it might not exist!
		if generegion in res:
			newid = generegion+"-result"
			comparatorGtf[newid] = [p["gtfline"] for p in res[generegion]]
			comparatorAa[newid]  = translateGtfLive(comparatorGtf[newid], cds)
			comparatorGeneregions[newid] = generegion
			comparatorSequences[newid]   = "result"
	# Grab the parts
	comparatorParts = grabComparatorParts(comparatorGtf, comparatorGeneregions, comparatorAa, path_dbl)
	sequences       = [makeSeq(c, comparatorAa[c]) for c in comparatorAa]
	# Make an alignment
	path_sequences  = path_dbl+"/comparator.fa"
	path_aln        = path_dbl+"/comparator.aln"
	writeSeqs(sequences, path_sequences)
	align(path_sequences, path_aln)
	# Grab the local coords for the alignment.	
	alignedseqs = dict((s.id, s) for s in readSeqs(path_aln))
	locsAa, locs, goodGtfs = flattenAlignment(alignedseqs.values(), comparatorGtf, path_aln + ".flat")
	return locs, comparatorParts, comparatorMcount, path_aln

def reviveCasualties(res, d_gr, d_si, path_wDir):
	# Designed to revive any genes that have been killed.
	casualties = [k for k in d_gr if not k in res]
	safe = [k for k in d_gr if k in res]
	if casualties:
		options = dict((k, [res[k]]) for k in safe)
		for k in casualties: options[k] = [makeMockPartSets(d_si[s]) for s in d_gr[k]["sequences"]]
		return scrapeBucket(options, path_wDir, d_gr)
	else:
		return res

def scrapeBucket(options, path_wDir, d_gr):
	# If there is only one sequences per generegion, return this list. Otherwise choose between.
	if all(len(s) == 1 for s in options.values()):
		return dict((k, options[k][0]) for k in options)
	else:
		casDir = makeIfAbsent(path_wDir+"/casualties")
		labels = idict(options, {})
		aas    = idict(options, [])
		for k in options:
			for l, s in enumerate(options[k]):
				labels[k][k + "."+ str(l)] = s
				aas[k].append(makeSeq(k + "."+ str(l), translatePart(getCdsLive([a["gtfline"] for a in s], d_gr[k]["cdsBase_data"]))))
		winners = chooseWinnersVanilla(aas, casDir)
		resn = {}
		for k in winners:
			resn[k] = labels[k][winners[k].id]
		return resn

def checkAqImprovement(resO, d_gr, d_si, path_wDir, orScore):
	seqs = {}
	aas = [makeSeq(k + "."+ str(random.random()), translatePart(getCdsLive([a["gtfline"] for a in resO[k]], d_gr[k]["cdsBase_data"]))) for k in resO]
	writeSeqs(aas, path_wDir+"/res.check.fa")
	align(path_wDir+"/res.check.fa", path_wDir+"/res.check.aln")
	alnseqs = readSeqs(path_wDir+"/res.check.aln")
	if alignmentScore(alnseqs) < orScore:
		for k in d_gr:
			seqs[k] = [makeMockPartSets(d_si[s]) for s in d_gr[k]["sequences"]]
		return scrapeBucket(seqs, path_wDir, d_gr)
	else: return resO

def makeMockPartSets(s):
	# All we really need is the gtflines.
	gtf = readCsv(s["gtfRel"])
	gtf = safeGtf(gtf)
	partslist = [{"gtfline": line} for line in gtf]
	return partslist

def compareWithOriginal(res, dict_generegions, dict_seqInfo, path_wDir):
	# Align the original and the new sequences.
	# Group them into overlap regions.
	# Choose winners for each as usual.
	sprint("Double-checking against original sequences...")
	path_dbl      = makeIfAbsent(path_wDir +"/doublecheck")
	dict_cdsBases = dict((k, dict_generegions[k]["cdsBase_data"]) for k in dict_generegions)
	# Prepare containers for all the information we're about to gather.
	# Grab the dna sequences for each generegion.
	cCoords, cParts, cMcount, path_aln = prepareCompare(res, dict_generegions, dict_seqInfo, path_dbl)
	groups, adj = groupParts(cCoords, ranges=True)
	# Do a first round of refinements by making a nice patchwork of old and new bits
	# based on what are the best choices.
	resRefined  = compareParts(adj, cParts, cMcount, path_dbl, dict_cdsBases, path_aln, minintron, minexon, dict_generegions)
	# The second round of refinements applies various filters to make sure nothing 
	# too exciting is going on. (Aq, Mg, Pl filters.)
	sprint("Applying filters...")
	path_fltrs  = makeIfAbsent(path_wDir+"/filters")
	cCoords, cParts, cMcount, path_aln = prepareCompare(resRefined, dict_generegions, dict_seqInfo, path_fltrs)
	groups, adj = groupParts(cCoords, ranges=True)
	resFiltered = filterChanges(adj, cParts, cCoords, path_fltrs, dict_cdsBases, path_aln, minintron, minexon, dict_generegions, cMcount)
	# We then do one more run of compare parts to double check that the filtered parts
	# work together to form viable genes. Shouldn't be too much to check here.
	return resFiltered

def grabComparatorParts(comparatorGtf, comparatorGeneregions, comparatorAa, path_dbl):
	comparatorParts = {}
	for newid in comparatorGtf:
		comparatorParts[newid] = {}
		fullAa = comparatorAa[newid]
		progress = 0
		gtflength = len(comparatorGtf[newid])
		for i,line in enumerate(sorted(comparatorGtf[newid], key = lambda x: int(x[3]))):
			part = {}
			cdslength 	   = int(line[4]) - int(line[3])-int(line[7]) + 1
			partlength	   = int(math.ceil(cdslength / 3.0))
			part["gtfline"]    = line
			part["part"]	   = fullAa[progress: progress + partlength]
			part["partId"]	   = newid + "." + str(i)
			part["initial"]	   = (i == 0)
			part["terminal"]   = (i == gtflength)
			part["partlength"] = partlength
			part["generegion"]   = comparatorGeneregions[newid]
			progress += partlength
			comparatorParts[newid][i] = part
	sequences = [makeSeq(c, comparatorAa[c]) for c in comparatorAa]
	path_sequences = path_dbl+"/comparator.fa"
	return comparatorParts

def initialiseEndOptions(part):
	part["optionsLeft"] = idict([0,1,2], [])
	part["optionsRight"] = idict([0,1,2], [])
	part["leftFlat"] = []
	part["rightFlat"] = []
	part["initialSites"] = []
	part["terminalSites"] = []
	part["donorSites"] = []
	part["acceptorSites"] = []
	return part

def initialiseProbe(parts, probeId):
	probe = copy.deepcopy(parts[0])
	probe["gtfline"][3] = min([p["gtfline"][3] for p in parts])
	probe["gtfline"][4] = max([p["gtfline"][4] for p in parts])
	probe["status"] = "fresh"
	probe["evidence"] = 0
	probe["id"] = probeId
	probe = initialiseEndOptions(probe)
	return probe

def addPartToProbe(probe, part, cds):
	gtf	 = part["gtfline"]
	frame 	 = int(gtf[7])
	endframe = getEndFrame(gtf)
	probe["ogtfline"] = part["gtfline"][:]
	probe["optionsLeft"][frame]	= l = list(set(probe["optionsLeft"][frame] + [int(gtf[3])]))
	probe["optionsRight"][endframe] = r = list(set(probe["optionsRight"][endframe] + [int(gtf[4])]))
	probe["evidence"] += 1
	probe["constituentIds"] = probe.get("constituentIds", []) + [part["partId"]]
	# Find out which bits are terminal, initial, donor, acceptor.
	if isInitialPart(gtf, cds):  probe["initialSites"]  += [int(gtf[3])]
	if isTerminalPart(gtf, cds): probe["terminalSites"] += [int(gtf[4])]
	if isAcceptorPart(gtf, cds): probe["acceptorSites"] += [int(gtf[3])]
	if isDonorPart(gtf, cds):    probe["donorSites"]    += [int(gtf[4])]
	# Add the flat options
	probe["leftFlat"] = list(set(l+probe["leftFlat"]))
	probe["rightFlat"] = list(set(r+probe["rightFlat"]))
	# Allow the options to be restored later if necessary
	probe["lbackupflat"] = copy.deepcopy(probe["leftFlat"])
	probe["rbackupflat"] = copy.deepcopy(probe["rightFlat"])
	probe["lbackup"] = copy.deepcopy(probe["optionsLeft"])
	probe["rbackup"] = copy.deepcopy(probe["optionsRight"])
	return probe

def createPart(q, status, startframe=-1, pid = ""):
	p   = copy.deepcopy(q)
	gtf = p["gtfline"]
	sf  = startframe if startframe != -1 else int(gtf[7])
	p["status"]          = status
	# Back up the options from before
	p["rbackupflat"]     = p["rightFlat"]
	p["lbackupflat"]     = p["leftFlat"]
	p["rbackup"]         = p["optionsRight"]
	p["lbackup"]         = p["optionsLeft"]
	# Replace the actual options with the new data.
	p["optionsLeft"]     = idict([0,1,2], [])
	p["optionsRight"]    = idict([0,1,2], [])
	p["leftFlat"]        = [int(gtf[3])]
	p["rightFlat"]       = [int(gtf[4])]
	p["optionsLeft"][sf] = [int(gtf[3])]
	p["optionsRight"][getEndFrame(gtf)] = [int(gtf[4])]
	if pid: p["id"] = pid
	return p

def prepareProbe(parts, i, cds):
	# Now we add in the options for each part based on what we have here.
	probe = initialiseProbe(parts, i)
	for p in parts: probe = addPartToProbe(probe, p, cds)
	# If there isn't sufficient representation at this point, we would like to
	# probe not having the part there at all.
	return probe

def comparePartsHeart(i,a, dict_parts, dict_mcount, path_refAln, d_reject, d_gr, res, tryempty, prevProbes, mi, mx, path_chkDir):
	# Create a set of peobes. Make sure that the only parts added in are parts that haven't already been added in.
	# If there are no parts available for a gene region, assume that neither new nor old has anything there.
	probes = {}
	for k in d_gr:
		parts = flattenLol([[copy.deepcopy(p) for p in q.values() if p["partId"] in a and p["generegion"] == k] for q in dict_parts.values()])
		# Make sure we haven't already added in these parts!
		if not k in res: res[k] = []
		if res[k]: parts = [copy.deepcopy(part) for part in parts if not any(["constituentIds" in q and part["partId"] in q["constituentIds"] for q in res[k]])]
		# If there are no parts, assume both the result and the original had nothing to offer.
		if not parts: continue
		probes[k]   = prepareProbe(parts, i, d_gr[k]["cdsBase_data"])
		tryempty[k] = probes[k]["evidence"] < dict_mcount[k]
	# Get the part combinations for option probing.
	# Need to be a little careful.
	dict_opts = {}; labels = {}
	write (d_reject, path_chkDir+"/reject.tpy")
	write(probes, path_chkDir+"/probes.tpy")
	write(prevProbes, path_chkDir+"/prevProbes.tpy")
	write(res, path_chkDir+"/plist.tpy")
	for k in res:
		plist  = copy.deepcopy(res[k]) if res[k] else []
		latest = [probes[k]] if k in probes else []
		if not plist and not latest: continue
		labels[k] = carefulPartCombinations(plist, latest, d_gr[k]["cdsBase_data"], k, mi, mx, tryempty[k] or not latest, reject=d_reject.get(k,[]))
	# Also try adding back in any probes from last time that were omitted.
	# Aims to mitigate mistakes.
	for k in res:
		allIds = [a["partId"] for a in res[k]]
		if k in prevProbes and k in res and not any([a in allIds for a in prevProbes[k]["constituentIds"]]):
			plist  = copy.deepcopy(res[k]) if res[k] else []
			latest = [prevProbes[k]] + ([] if not k in probes else [probes[k]])
			labels_l = carefulPartCombinations(plist, latest, d_gr[k]["cdsBase_data"], k, mi, mx, tryempty[k], reject=d_reject.get(k,[]), dr = path_chkDir)
			for l in labels_l: labels[k][l] = labels_l[l]
	dict_opts = dict((k,[l["seq"] for l in labels[k].values()]) for k in labels)
	return labels, dict_opts, probes, tryempty

def compareParts(adj, dict_parts, dict_mcount, path_wdir, dict_cdsBases, path_refAln, mi, mx, d_gr, d_reject = {}, refine=True):
	# Add in the parts as usual. At each stage pick which part works best in the context.
	# At each stage, if there are not k+1 options, where k is the total number of original
	# sequences for the generegion,  we assume removal of the part to be one of the options.
	# of the part. Note two identical options is still two options. Note also that we need to
	# pay attention to whether or not a part is initial, terminal, whatever.
	res      = idict(d_gr, {})
	tryempty = idict(d_gr, False)
	probes   = []
	# Now run through all the different part options in their groups, slot them in and
	# sequentially pick the best of them.
	for i, a in enumerate(adj):
		path_chkDir = makeIfAbsent(path_wdir + "/check_" + str(i))
		labels, options, probes, tryempty = comparePartsHeart(i, a, dict_parts, dict_mcount, path_refAln, \
											d_reject, d_gr, res, tryempty, probes, mi, mx, path_chkDir)
		write(labels, path_chkDir+"/labels.tpy")
		# Now align these lovely options and pick the best one.
		res, path_aln = processLabels(labels, options, path_chkDir, path_refAln, True, refine)
		sealDone(res)
	res = checkForStartsAndStops(res, dict_cdsBases, dict_parts, dict_mcount.keys(), path_wdir, mi, mx, d_reject)
	return res

def sealDone(res):
	for plist in res.values():
		for part in plist:
			if statusCheck(part["status"], ["done", "waiting"]):
				part["leftFlat"] = [int(part["gtfline"][3])]
			if statusCheck(part["status"], ["done"]):
				part["rightFlat"] = [int(part["gtfline"][4])]

def getScore(prevAln):
	# The lengths of the bands here might not be the same size. We correct for this
	# By adding in the extra bits of sequenc from the original alignment.
	# The order in which things are added back in doesn't matter, since the alignment
	# score is columnwise.
	squeezeMap = re.sub(r"bestaln", r"squeeze.fa", prevAln)
	if os.path.exists(squeezeMap):
		# There should be one for each gene region.
		filled = dict((getGr(s.id), s) for s in readSeqs(prevAln))
		for s in readSeqs(squeezeMap):
			generegion = getGr(s.id)
			if generegion in filled:
				filled[generegion].seq += Seq(str(s.seq).replace("@", ""))
		writeSeqs(filled.values(), os.path.dirname(prevAln)+"/all.filled.fa")
		return alignmentScore(filled.values(), scaled=False)
	else:
		return alignmentScore(readSeqs(prevAln), scaled=False)

def processLabels(labels, options, path_fDir, path_refAln, spanningBand, refine):
	winners  = chooseWinnersRef(options, path_fDir, path_refAln = path_refAln, doubleCheck = True, orAln = True, parallelCheck = True)
	res      = fetchByLabel(labels, winners)
	path_aln = path_fDir +"/result.aln"
	writeSeqs(winners.values(), path_aln)
	write(res, path_fDir+"/results.unrefined.typ")
	if refine: res = refineStatuses(res)
	write(res, path_fDir +"/results.tpy")
	return res, path_aln

def assignSiteTypes(part, gtf, cds):
	if isAcceptorPart(gtf, cds): part["acceptorSites"] += [int(gtf[3])]
	if isTerminalPart(gtf, cds): part["terminalSites"] += [int(gtf[4])]
	if isDonorPart(gtf, cds):    part["donorSites"]    += [int(gtf[4])]
	return part

def checkForStarts(res, dict_parts, cds, gtf, generegion, minintron, minexon, reject = []):
	# Get all the available initial parts
	parts = flattenLol([[copy.deepcopy(p) for p in q.values() if p["initial"] and p["generegion"] == generegion] for q in dict_parts.values()])
	labels_k = {}
	for opart in parts:
		part       = copy.deepcopy(opart)
		trialer    = copy.deepcopy(res[generegion])
		trialparts = [p for p in trialer if int(p["gtfline"][3]) > int(part["gtfline"][4])]
		if trialparts:
			trialparts[0]["leftFlat"]    = trialparts[0]["lbackupflat"]
			trialparts[0]["optionsLeft"] = trialparts[0]["lbackup"]
			part = initialiseEndOptions(part)
			part = createPart(part, "fresh", startframe = 0, pid = "-1")
			gtf  = part["gtfline"]
			# Add site type information
			part["initialSites"] = [int(gtf[3])]
			part = assignSiteTypes(part, gtf, cds)
			trialparts = [part] + trialparts
			labels_l = carefulPartCombinations(trialparts, [], cds, generegion, minintron, minexon, reject = reject)
			for l in labels_l: labels_k[l] = labels_l[l]
		# Also try merging the new part with any old parts that overlap it.
		part    = copy.deepcopy(opart)
		trialer = copy.deepcopy(res[generegion])
		# We want to merge the part with the last part that occupies the same frame
		trialparts = [p for p in trialer if int(p["gtfline"][3]) > int(part["gtfline"][4])]
		remaining =  [p for p in trialer if not (int(p["gtfline"][3]) > int(part["gtfline"][4]))]
		found = False
		for p in reversed(remaining):
			if gtfLinesOverlap(part["gtfline"], p["gtfline"]) and gtfCompatibleFrame(part["gtfline"], p["gtfline"]):
				pp = copy.deepcopy(p)
				pp["status"]         = "fresh"
				gtf                  = part["gtfline"]
				pp["gtfline"][3]     = gtf[3]
				pp["optionsLeft"][0] = [int(gtf[3])]
				pp["leftFlat"]       = [int(gtf[3])]
				pp["initialSites"]   = [int(gtf[3])]
				trialparts = [pp] + trialparts
				found = True
				break
			else:
				 trialparts.append(p)
		if found:
			labels_l = carefulPartCombinations(trialparts, [], cds, generegion, minintron, minexon, reject = reject)
			for l in labels_l: labels_k[l] = labels_l[l]
	return labels_k

def checkForStops(dict_cdsBases, dict_parts, labels_k, cds, generegion, minintron, minexon, reject = []):
	labels_r = {}
	parts = flattenLol([[copy.deepcopy(p) for p in q.values() if isTerminalPart(p["gtfline"], dict_cdsBases[generegion]) and p["generegion"] == generegion] for q in dict_parts.values()])
	for opart in parts:
		for label in labels_k:
			trialer = copy.deepcopy(labels_k[label])
			# Try removing any end parts and shoehorning in the probe part
			part = copy.deepcopy(opart)
			trialparts = [p for p in trialer["parts"] if int(p["gtfline"][4]) < int(part["gtfline"][3])]
			if trialparts:
				trialparts[-1]["rightFlat"]    = trialparts[-1]["rbackupflat"]
				trialparts[-1]["optionsRight"] = trialparts[-1]["rbackup"]
				trialparts[-1]["status"]       = "waitingH"
				part = initialiseEndOptions(part)
				gtf  = part["gtfline"]
				part["optionsLeft"][0]                 = [int(gtf[3])]
				part["optionsRight"][getEndFrame(gtf)] = [int(gtf[4])]
				part["leftFlat"]                       = [int(gtf[3])]
				part["rightFlat"]                      = [int(gtf[4])]
				part["status"]                         = "fresh"
                                part["lbackupflat"]                    = part["leftFlat"]
                                part["rbackupflat"]                    = part["rightFlat"]
				part["rbackup"]                        = copy.deepcopy(part["optionsRight"])
				part["lbackup"]                        = copy.deepcopy(part["optionsLeft"])
				# Add site type information
				if isInitialPart(gtf, cds): part["initialSites"]   += [int(gtf[3])]
				if isAcceptorPart(gtf, cds): part["acceptorSites"] += [int(gtf[3])]
				part["terminalSites"] = [int(gtf[4])]
				part["donorSites"]    = []
				part["id"] = max([t["id"] for t in trialparts]) + 1
				trialparts = trialparts + [part]
				labels_l   = carefulPartCombinations(trialparts, [], cds, generegion, minintron, minexon, reject = reject)
				for l in labels_l: labels_r[l] = labels_l[l]
			# Also try merging the new part with any old parts that overlap it.
			part    = copy.deepcopy(opart)
			trialer = copy.deepcopy(labels_k[label])
			# We want to merge the part with the last part that occupies the same frame
			trialparts = [p for p in trialer["parts"] if int(p["gtfline"][4]) < int(part["gtfline"][3])]
			remaining  = [p for p in trialer["parts"] if not (int(p["gtfline"][4]) < int(part["gtfline"][3]))]
			found = False
			for p in remaining:
				if gtfLinesOverlap(part["gtfline"], p["gtfline"]) and gtfCompatibleFrame(part["gtfline"], p["gtfline"]):
					pp = copy.deepcopy(p)
					pp["status"]        = "fresh"
					gtf                 = part["gtfline"]
					pp["gtfline"][4]    = gtf[4]
					pp["optionsRight"][getEndFrame(gtf)] = [int(gtf[4])]
					pp["rightFlat"]     = [int(gtf[4])]
					pp["terminalSites"] = [int(part["gtfline"][4])]
					trialparts.append(pp)
					found = True
					break
				else:
					 trialparts.append(p)
			if found:
				labels_l = carefulPartCombinations(trialparts, [], cds, generegion, minintron, minexon, reject = reject)
				for l in labels_l: labels_r[l] = labels_l[l]					
	return labels_r

def checkForStartsAndStops(res, dict_cdsBases, dict_parts, generegions, path_dblDir, minintron, minexon, d_reject = {}):
	# It is possible that the result of this checking process might lack a true initial or terminal exon.
	# We check at this stage whether this has taken place. We then shoehorn in all the options for initial
	# and terminal exons and pick the best ones.
	# In addition to the shoehorning, we accept extensions and contractions of existing parts.
	replace = False
	labels  = idict(generegions, {})
	for k in res:
		# Decide whether the gene model needs to be replaced.
		# If the result is blank (i.e. recommend remove gene), then skip this generegion.
		if not res[k]: continue
		# Check if any of the parts are initial
		partssorted = sorted(res[k], key = lambda x: int(x["gtfline"][3]))
		gtf         = [p["gtfline"] for p in partssorted]
		cds         = dict_cdsBases[k]
		gtfcds      = getCdsLive(gtf, cds)
		firstcodon  = gtfcds[:3]
		gtfcdsframe = int(gtf[0][7])
		endframe    = getEndFrame(gtf[-1])
		# If the first codon is not a start codon, try adding in some starts.
		labels_k = {}
		if not (isStartCodon(firstcodon) and gtfcdsframe == 0):
			replace = True
			labels_k = checkForStarts(res, dict_parts, cds, gtf, k, minintron, minexon, d_reject.get(k,[]))
		else:
			aaseq = seqFromPartsList(res[k], cds)
			labels_k = {"stillgood": {"parts": res[k], "aa": aaseq, "seq": makeSeq("stillgood", aaseq), "generegion": k}}
		# Check if the last part is terminal
		lastcodon = gtfcds[-3:]
		if not (isStopCodon(lastcodon) and endframe == 0):
			replace = True
			# Get all the available terminal parts
			labels_r = checkForStops(dict_cdsBases, dict_parts, labels_k, cds, k, minintron, minexon, d_reject.get(k,[]))
			for l in labels_r: labels_k[l] = labels_r[l]	
			if "stillgood" in labels_k: del labels_k["stillgood"]
		labels[k] = labels_k
	# If the gene model does not need to be replaced, simply return it.
	# If it does, then choose between the options for replacement.
	# It is possible that no viable gene model will be found.
	if replace:
		dict_opts   = dict((k, [l["seq"] for l in labels[k].values()])for k in generegions)
		path_chkDir = makeIfAbsent(path_dblDir + "/check_final")
		winners     = chooseWinnersSeeded(dict_opts, path_chkDir, seed = "auto")
		for generegion in labels:
			if not generegion in winners: continue
			w = labels[generegion][winners[generegion].id]
			for k in w: res[generegion] = w["parts"]
			writeCsv([p["gtfline"] for p in res[generegion]], path_chkDir + "/" + generegion + ".gtf")
		writeSeqs(winners.values(), path_chkDir +"/result.aln")
	return res

def carefulPartCombinations(plist, latestparts, cds, generegion, mi, mx, tryempty=False, reject = [], dr=""):
	# Try not adding in the part, if relevant.
	labels_e = getLabels(plist, cds, reject, generegion, mi, mx, True) if tryempty or not latestparts else {}
	# Then try adding in the part.
	labels_l = {}
	if latestparts:
		plist += latestparts
		labels_l = getLabels(plist, cds, reject, generegion, mi, mx, True)
		# Try to foresee problems with non-viable part sets
		if len(plist) > 1:
			last  = plist[-1]["gtfline"]
			penlt = plist[-2]["gtfline"]
			# If the last two lines overlap or the last is before the first,
			# try removing each of them. Note we don't need to remove the last
			# one if tryempty is on, that will already have been covered
			# Update: only try removing last part if the part in question is terminal.
			if gtfLinesOverlap(last, penlt) or int(last[4]) < int(penlt[3]):
				if (not tryempty) and plist[-1]["terminalSites"]:
					labels_c = getLabels(delIndices(plist,[-1]), cds, reject, generegion, mi, mx, True)
					for l in labels_c: labels_l[l] = labels_c[l]
				labels_d = getLabels(delIndices(plist,[-2]), cds, reject, generegion, mi, mx, True)
				for l in labels_d: labels_l[l] = labels_d[l]
	# Any time we have a start codon, be sure to squash all the previous stuff
	# and try out the initial part on its own. (The part on its own should indeed
	# suffice - we shouldn't need to retain anything after it.)
	for lp in [a for a in latestparts if a["initialSites"]]:
		labels_s = getLabels([lp], cds, reject, generegion, mi, mx, True)
		for l in labels_s: labels_l[l] = labels_s[l]
	return mergeDicts([labels_l, labels_e])

def getLabels(plistc, cds, reject, generegion, mi, mx, chst):
	opts, plist2 = getPartCombinationsSafe(plistc, cds, reject)
	return getCombinationsFasta(plist2, opts, cds, generegion, mi, mx, checkStart = chst)

def getPartCombinationsSafe(partslist, cds, reject = []):
	# The latest "waiting" part and any subsequent ones should be fully laid out.
	plist = copy.deepcopy(partslist)
	if [p["id"] for p in plist if p["status"].startswith("waiting")]:
		maxWaitingId = max([p["id"] for p in plist if p["status"].startswith("waiting")])
		for p in plist:
			if maxWaitingId == p["id"]:
				p["rightFlat"]    = copy.deepcopy(p["rbackupflat"])
				p["optionsRight"] = copy.deepcopy(p["rbackup"])
			if maxWaitingId < p["id"]:
				p["status"] = "doubleQ"
				p["rightFlat"]    = copy.deepcopy(p["rbackupflat"])
				p["optionsRight"] = copy.deepcopy(p["rbackup"])
				p["optionsLeft"] = copy.deepcopy(p["lbackup"])
				p["leftFlat"]    = copy.deepcopy(p["lbackupflat"])
	return getPartCombinations(plist, cds, checkInitTerm=True, reject = reject), plist


#####################################################
# Execution
#####################################################

def go(path_inf, path_ref, path_resultsDir, path_wDir, minintron, minexon, int_numcores, int_slopAmount):
	"""The main body of the program.
	"""
	#########################################################
	# Set up basic variables and data containers. Extract
	# information from the genome files.
	#########################################################
	dict_seqInfo    = readInputLocations(path_inf, path_ref)#qe
	dict_genomeInfo = getGenomeInfo(dict_seqInfo)#qe
	#########################################################
	# Get the base region for each gtf, slopping a certain
	# amount each side. Also get the fasta sequences.
	#########################################################
	sprint("Grabbing cds and aa sequences for inputted transcripts...")
	prepareSequences(path_wDir, dict_seqInfo, dict_genomeInfo, int_numCores, int_slopAmount)#qe
	#########################################################
	# For each genome, merge together overlapping tight base
	# regions. This will allow us to, bluntly, auto-detect
	# alternative transcripts in a region
	#########################################################
	dict_generegions = prepareGeneregions(dict_seqInfo, dict_genomeInfo, path_wDir, int_numCores,int_slopAmount)#qe
	################################################################
	# Grab the relative boundaries for the original sequences
	################################################################
	sprint("Relativising sequences...")
	relativiseSequences(dict_generegions, dict_seqInfo)#qe
	#############################################################
	# For each base, exonerate each protein sequence against it.
	# Then grab the gene parts.
	#############################################################
	goExonerate(path_wDir, dict_generegions, dict_seqInfo, int_numCores)#qe
	#############################################################
	# Align all options and pick the best set
	#############################################################
	winners = prepareOptions(dict_generegions)#qe
	######################################################
	# Align the original sequences to use as a reference
	######################################################
	path_orFa       = path_wDir + "/or.fasta"#qe
	path_orAln      = path_wDir + "/or.aln"#qe
	getOrAln(dict_generegions, dict_seqInfo, path_orFa, path_orAln)#qe
	####################################################
	# Align the winners and get the coordinates.
	####################################################
	path_winnersAln = prepareWinnerAlignment(path_wDir, winners, dict_generegions, path_orFa)#qe
	dict_partCoords = getWinnerAlnCoords(winners, path_winnersAln, dict_generegions, path_winnersAln+".flat")#qe
	###################################################
	# Get a connectivty graph for these bits.
	###################################################
	groups, adj = groupParts(dict_partCoords)#qe
	###############################################################
	# Fix each pair in turn, as directed by the connectivity graph
	###############################################################
	parts = dict((k, dict_generegions[k]["winner"]) for k in dict_generegions)#qe
	res   = fixIt(adj, parts, dict_generegions, path_wDir, minintron, minexon, path_winnersAln)#qe
	###############################################################
	# Production code to quickly jump to the doublecheck section.
	###############################################################
	path_pDir = makeIfAbsent(path_wDir +"/pickle")
	# Pickle.
	pickle.dump(res, open(path_pDir +"/res.p", "w"))#qe
	pickle.dump(dict_generegions, open(path_pDir +"/met.p", "w"))#qe
	pickle.dump(dict_seqInfo, open(path_pDir +"/seq.p", "w"))#qe
	# Unpickle
	res = pickle.load(open(path_pDir +"/res.p", "r"))#qg
	dict_generegions = pickle.load(open(path_pDir +"/met.p", "r"))#qg
	dict_seqInfo = pickle.load(open(path_pDir +"/seq.p", "r"))#qg
	##################################################
	# Compare with the original sequences
	##################################################
	resO = compareWithOriginal(res, dict_generegions, dict_seqInfo, path_wDir)#qg
	##################################################
	# If the process yields nothing, retain the
	# original gene model.
	##################################################
	resO = reviveCasualties(resO, dict_generegions, dict_seqInfo, path_wDir)
	resO = checkAqImprovement(resO, dict_generegions, dict_seqInfo, path_wDir, alignmentScore(readSeqs(path_orAln)))
	##################################################
	# Output the gtfs, fastas, etc. etc. etc.
	##################################################
	dict_generegions = writeOutput(dict_generegions, dict_seqInfo, resO, path_resultsDir)#qg
	getStats(dict_generegions, dict_seqInfo, res, path_resultsDir, path_wDir)

if __name__ == '__main__':
	"""OMGene - mutual improvement of gene models through gene orthology
	"""
	# Check the shell...
	checkShell()
	
	# Read in the command line arguments
	parser = argparse.ArgumentParser(description="Run OMGene")
	parser.add_argument("-i", "--infoFiles", metavar="infiles", dest="IN")
	parser.add_argument("-r", "--referenceFiles", metavar="reffiles", dest="RF", default="")
	parser.add_argument("-o", "--outDir", metavar="outdir", dest="OD")
	parser.add_argument("-g", "--gapPenalty", metavar="outdir", dest="GP")
	parser.add_argument("--minintron", metavar="minintron", dest="MNI", default=4)
	parser.add_argument("--minexon", metavar="minexon", dest="MNE", default=4)
	parser.add_argument("-c", "--cores", metavar="numcores", dest="CO", default=1)
	parser.add_argument("-b", "--buffer", metavar="bufferamount", dest="SL", default=600)
	parser.add_argument("--safealign", action="store_true", dest="SA", default=False)
	parser.add_argument("--donors", "-do", metavar="donors", dest="DO", default="gt,gc")
	parser.add_argument("--acceptors", "-ac", metavar="acceptors", dest="AC", default="ag")
	parser.add_argument("--startcodon", "-sc", metavar="startcodons", dest="SC", default="atg")
	args = parser.parse_args()
	#CUPCAKES: dependencies: pyfaidx, exonerate, blastp
	#Check existence and non-confliction of arguments
	path_inf       = args.IN
	path_ref       = args.RF
	path_outDir    = args.OD
	if not path_inf: sys.exit("Please specify an input file using the -i option")
	minexon	= int(args.MNE)
	minintron      = int(args.MNI)
	int_numCores   = int(args.CO)
	int_slopAmount = int(args.SL)
	static_safeAlign = args.SA
	
	static_do, static_ac, static_sc = processCodonOptions(args.DO, args.AC, args.SC)

	if args.GP:
		static_blosdict = goodBlosum(-1*int(args.GP), 0)
		static_blosmat  = blosMat(static_blosdict)
		static_blospos  = dict((e,i) for i,e in enumerate(static_aa))
	path_resultsDir, path_wDir = prepareOutputFolder(path_outDir)
	go(path_inf, path_ref, path_resultsDir, path_wDir, minintron, minexon, int_numCores, int_slopAmount)
