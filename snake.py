import os
import glob
import math 

localrules: final, all, split
configfile: "config.yaml"

prefix = os.path.abspath(config["output_dir"]) + "/"

lines = open(config["reads_fofn"]).readlines()
nlines = len(lines)
if("files_per_job" in config):
	files_per_job = int(config['files_per_job'])
else:
	files_per_job = 2000

maxidx = math.ceil( (nlines*1.0)/files_per_job ) 
IDS = list(range(0,maxidx))
print(IDS)

rule all:
	input:
		prefix + "final"

rule split:
	input:
		config["reads_fofn"]
	output:
		fofns = expand(prefix + "fofns/{ID}.fofn", ID=IDS),
	run:
		lines = open(config["reads_fofn"]).readlines()
		count = 0
		ID = 0
		out = ""
		for line in lines:
			count += 1
			out += line
			if(count == files_per_job):
				open(prefix + "fofns/{}.fofn".format(ID), "w+").write(out)
				count = 0
				ID += 1
				out = ""
		open(prefix + "fofns/{}.fofn".format(ID), "w+").write(out)


rule basecall:
	input:
		fofn = prefix + "fofns/{ID}.fofn"
	output:
		reads = prefix + "reads/{ID}.done"
	params:
		cluster=" -l mfree=2G -pe serial 4 "
	threads: 4
	shell:"""
module purge
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs/prod modules-eichler
module load python/3.6.4
module load albacore/2.3.3

outdir={prefix}reads/{wildcards.ID}/
mkdir -p $outdir

cat {input.fofn} | \
	read_fast5_basecaller.py \
	-s $outdir -t {threads} -k SQK-RAD004 -f FLO-MIN106 

touch {output}
"""




rule final:
	input:
		reads = expand(prefix + "reads/{ID}.done", ID=IDS),
	output:
		prefix + "final"
	shell:
		"touch {output}"


