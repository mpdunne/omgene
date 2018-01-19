#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################
# Extracts genes from a gff file.
################################################

import csv
import sys
import re
import random
import string
import os

def readCsv(path_file, ignoreBlank = True, ignoreHash = True, grabHeader = False):
	with open(path_file, "r") as f:
		data = [line for line in csv.reader(f, delimiter="\t") if (''.join(line).strip() or not ignoreBlank) and (not line[0].startswith("#") or not ignoreHash)]
	if not grabHeader:
		return data
	with open(path_file, "r") as f:
		header = [line for line in csv.reader(f, delimiter="\t") if (''.join(line).strip() or not ignoreBlank) and line[0].startswith("#")]
	return data, header

def writeCsv(rows, path_file):
        with open(path_file, "w") as f:
                writer = csv.writer(f, delimiter="\t", quoting = csv.QUOTE_NONE, quotechar='')
                writer.writerows(rows)

args = sys.argv[1:]
path_gffin = args[0]
prefix = args[1]

if not os.path.exists(path_gffin):
	print "Path " + path_gffin + " does not exist!"

if not re.search("gff3?$", path_gffin):
	print "File does not appear to be in gff3 format"

gff = readCsv(path_gffin)

# We're only going to be interested in "gene" entries.
# Construct a family tree for the GFF entries.

# We're going to require that transcripts have the format gene > mrna > exon/cds

d_down = {}; d_up   = {}; d_ids  = {}
genes  = []; lmrnas = []; cdss   = []; exons  = []

for line in gff:
	# The line might not have an id tag. If it doesn't, make one up.
	eid = re.sub(r".*ID=([^=;]*).*", r"\1", line[8]) if "ID" in line[8] else "id_"+"".join([random.choice(string.letters) for i in range(0,10)])
	if "Parent" in line[8]:
		# Multiple parents can be listed.
		parents = re.sub(r".*Parent=([^=;]*).*", r"\1", line[8]).split(",")
		for p in parents:
			d_down[p] = list(set(d_down.get(p, []) + [eid]))
		d_up[eid] = list(set(d_up.get(eid, []) + parents))
	if line[2].lower() == "gene":  genes.append(eid)
	if line[2].lower() == "mrna":  lmrnas.append(eid)
	if line[2].lower() == "cds":   cdss.append(eid)
	if line[2].lower() == "exon": exons.append(eid)
	d_ids[eid] = d_ids.get(eid, []) + [line]

mrnas = {}
genenames = {}
mrnanames = {}
mrnagenes = {}
strands = {}

transcriptcount = {}
gffout = []

for g in genes:
	print "Processing GFF3 entry for " + g + "..."
	if g in d_down:
		children_mrna = list(set(d_down[g]))
		# Each of the genes direct children should be tagged mrna.
		# Expect one entry for each gene, one entry for each mrna, potential
		# multiple entries for each mnra child.
		g_lines = d_ids[g]
		if len(g_lines) !=1: continue
		# Sort out the name.
		if "Name=" in g_lines[0][8]:
			genename = re.sub(r".*Name=([^=;]*).*", r"\1", g_lines[0][8])
		else:
			genename = g
		# Ain't got no time for partial genes!!
		if "partial=true" in g_lines[0][8]: continue
		# These genes have got to be protein coding!!!
		if not "gene_biotype=protein_coding" in g_lines[0][8]: continue
		gffout += g_lines
		genenames[g] = genename
		g_strand = g_lines[0][6]
		for m in children_mrna:
			if not m in lmrnas: continue
			m_lines = d_ids[m]
			# Should only be one line per mrna, if not, abandon.
			if len(m_lines) != 1: continue
			m_strand = m_lines[0][6]
			if m_strand != g_strand: continue
			strands[m] = m_strand
			gffout += m_lines
			if not m in d_down: continue
			children_exonsandcds = d_down[m]
			# The transcript name is the gene name plus a number.	
			tcount = transcriptcount[genename] = transcriptcount.get(genename, 0) + 1
			mrnanames[m] = genename + ".t"+str(tcount)
			mrnagenes[m] = g
			for x in children_exonsandcds:
				if not ((x in cdss or x in exons) and m in lmrnas): continue
				x_lines = d_ids[x]
				gffout += x_lines
				mrnas[m] = mrnas.get(m, []) + x_lines

path_gffout = re.sub(r"gff3?", r"gene_exons.gff3", path_gffin)
writeCsv(sorted(gffout, key = lambda x: int(x[3])), path_gffout)

# Now need to add start and stop codons to the mrna annotations for the gtf.
gtf = {}
for m in mrnas:
	print "Processing GTF entry for " + m + "..."
	lines = mrnas[m]
	gid = genenames[mrnagenes[m]]
	tid = mrnanames[m]
	# Exon IDs are added per exon. CDS regions contained in each
	# exon inherit the relevant exon ID.
	exons = [l for l in lines if l[2].lower() == "exon"]
	cds = [l for l in lines if l[2].lower() == "cds"]
	# Starp codons are annotated according to beginnings and endings
	# of CDS regions. Do things differently depending on strand.
	strand = strands[m]
	# The starp codon annotations here don't account for spliced starts and stops.
	# ignore this for now.
	gtf[m] = []
	tagbase = "transcript_id \""+prefix+"."+tid+"\"; gene_id \""+prefix+"."+gid+"\"; gene_name \""+gid+"\""
	basis  = [cds[0][0], cds[0][1], "x", "y", "z", cds[0][5], strand, "w", tagbase]
	if strand == "+":
		cds = sorted(cds, key = lambda x: int(x[3]))
		exons = sorted(exons, key = lambda x: int(x[3]))
		startcoords = [cds[0][3], str(int(cds[0][3])+2)]
		stopcoords = [str(int(cds[-1][4])-2), str(int(cds[-1][4]))]
		exonids = dict((i, exons[i]) for i in range(len(exons)))
		# Construct lines.
		startline = basis[:]; stopline = basis[:]
		startline[3] = startcoords[0]; startline[4] = startcoords[1]; startline[2] = "start_codon"; startline[7] = 0
		stopline[3]  = stopcoords[0];  stopline[4]  = stopcoords[1]; stopline[2] = "stop_codon"; stopline[7] = 0
		gtf[m] = [startline, stopline]
		for c in cds:
			c[8] = tagbase
			gtf[m] += [c]
		for x in exons:
			x[8] = tagbase
			gtf[m] += [x]
	if strand == "-":
		cds = sorted(cds, key = lambda x: -int(x[3]))
		exons = sorted(exons, key = lambda x: int(x[3]))
		startcoords = [str(int(cds[-1][4])-2), str(int(cds[-1][4]))]
		stopcoords = [cds[0][3], str(int(cds[0][3])+2)]
		exonids = dict((i, exons[i]) for i in range(len(exons)))
		# Construct lines.
                startline = basis[:]; stopline = basis[:]
		startline[3] = startcoords[0]; startline[4] = startcoords[1]; startline[2] = "start_codon"; startline[7] = 0
		stopline[3]  = stopcoords[0];  stopline[4]  = stopcoords[1]; stopline[2] = "stop_codon"; stopline[7] = 0
		gtf[m] = [startline, stopline]
		for c in cds:
			c[8] = tagbase
			gtf[m] += [c]
		for x in exons:
			x[8] = tagbase
			gtf[m] += [x]

path_gtfout = re.sub(r"gff3?", r"gene_exons.gtf", path_gffin)

gtfout = []
for m in gtf: gtfout += sorted(gtf[m], key = lambda x: int(x[3]))

writeCsv(gtfout, path_gtfout)
