#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 22 March 2019

# Description: This script performs a given blast process from a fasta query to a fasta-derived database.
# It splits the fasta file into several jobs and merges the outputs.
########################################################################
"""

################################################################################
"""---1.0 Import Modules---"""

################################################################################
"""---2.0 Define Functions---"""



import sys, argparse, os
import time, glob
import subprocess

"""--------------------1.0 Initialize Variables---------------------------"""


parser = argparse.ArgumentParser(description=
	'Global mandatory parameters: [query_file] [database_file] [program] [number_Jobs] [db_type]\n'
	'Optional Database Parameters: See '+sys.argv[0]+' -h')

parser.add_argument("-q", "--query", required=True, help="Query File")
parser.add_argument("-j", "--jobs", help="Number of Jobs to Run", type=int, default=None)
parser.add_argument("-p", "--program", help="Program to run [blastn, blastp, blastx, tblastx or blat]", default="blastn")
parser.add_argument("-e", "--output", required=True, help="Output Folder")
parser.add_argument("--step", help="Step to begin with: [1] FormatDB, [2] Split Query, [3] Blast Run,", type=int, default=1)
parser.add_argument("-f", "--force", help="Force overwriting database if already formatted", action='store_true')
parser.add_argument("-s", "--split", help="Force overwriting split query if already found", action='store_true')

"""--------------------1.1 Database Variables ----------------------------"""

parser.add_argument("-d", "--database", required=True, help="Database File to Format")
parser.add_argument("--dbtype", required=True, help="Molecule type of target db [nucl, prot]")
parser.add_argument("--db_opt", nargs='*', dest='db_opt',  help="""Enter the options as option_flag value with no dash "-" (e.g. [input_type fasta])
	[-input_type type] [-dbtype molecule_type]
	[-title database_title] [-parse_seqids] [-hash_index] [-mask_data mask_data_files]
	[-mask_id mask_algo_ids] [-mask_desc mask_algo_descriptions] [-gi_mask]
	[-gi_mask_name gi_based_mask_names] [-out database_name] [-max_file_sz number_of_bytes]
	[-taxid TaxID] [-taxid_map TaxIDMapFile] [-logfile File_Name]\n""")

"""--------------------1.2 Blast Variables -------------------------------"""

parser.add_argument("--blast_opt", dest='blast_opt', nargs='*', default="-outfmt 6", help="""Enter the options as option_flag value with no dash "-" (e.g. [outfmt 6])
	[-import_search_strategy filename] [-export_search_strategy filename] [-task task_name] [-db database_name]
	[-dbsize num_letters] [-gilist filename] [-seqidlist filename]
	[-negative_gilist filename] [-entrez_query entrez_query]
	[-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
	[-subject subject_input_file] [-subject_loc range] [-query input_file]
	[-out output_file] [-evalue evalue] [-word_size int_value]
	[-gapopen open_penalty] [-gapextend extend_penalty]
	[-perc_identity float_value] [-xdrop_ungap float_value]
	[-xdrop_gap float_value] [-xdrop_gap_final float_value]
	[-searchsp int_value] [-max_hsps int_value] [-sum_statistics]
	[-penalty penalty] [-reward reward] [-no_greedy]
	[-min_raw_gapped_score int_value] [-template_type type]
	[-template_length int_value] [-dust DUST_options]
	[-filtering_db filtering_database]
	[-window_masker_taxid window_masker_taxid]
	[-window_masker_db window_masker_db] [-soft_masking soft_masking]
	[-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
	[-best_hit_score_edge float_value] [-window_size int_value]
	[-off_diagonal_range int_value] [-use_index boolean] [-index_name string]
	[-lcase_masking] [-query_loc range] [-strand strand] [-parse_deflines]
	[-outfmt format] [-show_gis] [-num_descriptions int_value]
	[-num_alignments int_value] [-html] [-max_target_seqs num_sequences]
	[-num_threads int_value] [-remote]""")

args = parser.parse_args()

"""Parse all the inputs and extract the target and origin folders"""

# Run directory
rundir = os.getcwd()

# Query information
query = args.query
query_dir = os.path.abspath(os.path.dirname(query))
query_name = os.path.basename(query)
query_base = os.path.splitext(query_name)[0]
query_fullpath = os.path.join(query_dir, query_name)

# Number of jobs
jobs = args.jobs if args.jobs else input("Number of Jobs to Run: ")
jobs = str(jobs)

# Program to run
program = args.program

# Database information
database = args.database
database_dir = os.path.abspath(os.path.dirname(database))
database_name = os.path.basename(database)
database_base = os.path.splitext(database_name)[0]
database_fullpath = os.path.join(database_dir, database_name)

dbtype = args.dbtype
force = args.force

# Other options
db_opt = args.db_opt

blast_opt = args.blast_opt

# Step to run
step = args.step

# Output information
output = args.output
output_dir = os.path.abspath(os.path.dirname(output))
output_name = os.path.basename(output)
output_fullpath = os.path.join(output_dir, output_name)


"""--------------------2.0 Format Database-------------------------------"""
if step == 1:
	if not os.path.exists(output_fullpath):
		os.makedirs(output_fullpath)
	os.chdir(rundir)
	print("Formatting Database: " + database_name)

# Create the format DB script.
	def formatDatabase():
		text_file = open(output_fullpath + "/Format" + database_base + "_DB.pbs", "w")
		text_file.write("#!/bin/bash\n")
		text_file.write("#PBS -N Database_Format\n")
		text_file.write("#PBS -l nodes=1:ppn=1\n")
		text_file.write("#PBS -l mem=5gb\n")
		text_file.write("#PBS -l walltime=3:00:00\n")
		text_file.write("#PBS -q iw-shared-6\n")
		text_file.write("#PBS -j oe\n")
		text_file.write("#PBS -o " + database_base + "_Format.out\n")
		text_file.write("#PBS -e " + database_base + "_Format.err\n\n")
		text_file.write("\ncd " + rundir + "\n\n")
		text_file.write("\nmodule load boost/1.53.0\nmodule load python/2.7\nmodule load ncbi_blast/2.2.29\n")
		text_file.write("\nmakeblastdb -in " + database_fullpath + " -dbtype " + dbtype + " -title " + database_base + " -out " + database_dir + "/" + database_base + " -logfile " + database_dir + "/DB_Format_Log.txt\n")
		text_file.close()

	formatDatabase()

	# Check if database exists and if overwrite is active
	if os.path.exists(database_dir + "/" + database_base + ".phr"):
		if force == True:
			try:
				os.remove(rundir + "/" + database_base + "_Format.out")
			except OSError:
				pass
			os.system("qsub " + output_fullpath + "/Format" + database_base + "_DB.pbs")
			for i in range(1000):
				if os.path.exists(rundir + "/" + database_base +"_Format.out"):
					break
				else:
					print("Still Formatting Database")
					time.sleep(30)
					continue
			os.system("mv " + rundir + "/" + database_base + "_Format.out " + output_fullpath)
	else:
		os.system("qsub " + output_fullpath + "/Format" + database_base + "_DB.pbs")
		for i in range(1000):
			if os.path.exists(rundir + "/" + database_base +"_Format.out"):
				break
			else:
				print("Still Formatting Database")
				time.sleep(30)
				continue
		os.system("mv " + rundir + "/" + database_base + "_Format.out " + output_fullpath)
	step += 1

"""------------------3.0 Split Query ----------------------"""

if step == 2:
	if not os.path.exists(output_fullpath):
		os.makedirs(output_fullpath)
	print("Splitting Query: " + query_name)

	def SplitQuery():
		divide = open(output_fullpath + "/Split.py", "w")
		divide.write("#!/usr/bin/env python\nfrom Bio import SeqIO\nimport sys, math\n")
		divide.write("record_iter = SeqIO.parse(open(\"" + query_fullpath + "\", \"r\"),\"fasta\")\n")
		divide.write("number=0\n")
		divide.write("for record in record_iter:\n")
		divide.write("\tnumber += 1\n")
		divide.write("num_seqs = int(math.ceil(number/float(" + jobs + ")))\n")
		divide.write("def batch_iterator(iterator, batch_size) :\n")
		divide.write("\tentry = True\n")
		divide.write("\twhile entry :\n")
		divide.write("\t\tbatch = []\n")
		divide.write("\t\twhile len(batch) < batch_size :\n")
		divide.write("\t\t\ttry :\n")
		divide.write("\t\t\t\tentry = next(iterator)\n")
		divide.write("\t\t\texcept StopIteration :\n")
		divide.write("\t\t\t\tentry = None\n")
		divide.write("\t\t\tif entry is None :\n")
		divide.write("\t\t\t\tbreak\n")
		divide.write("\t\t\tbatch.append(entry)\n")
		divide.write("\t\tif batch :\n")
		divide.write("\t\t\tyield batch\n")
		divide.write("record_iter = SeqIO.parse(open(\"" + query_fullpath + "\", \"r\"),\"fasta\")\n")
		divide.write("for i, batch in enumerate(batch_iterator(record_iter, num_seqs)) :\n")
		divide.write("\tfilename = \"" + query_dir + "/" + query_base + "_Split/" + query_base + "_%i.fasta\" % (i+1)\n")
		divide.write("\thandle = open(filename, \"w\")\n")
		divide.write("\tcount = SeqIO.write(batch, handle, \"fasta\")\n")
		divide.write("\thandle.close()\n")
		divide.write("\tprint(\"Wrote %i records to %s\" % (count, filename)) \n")
		divide.close()


		split_pbs = open(output_fullpath + "/Split.pbs", "w")
		split_pbs.write("#!/bin/bash\n#PBS -N Split_Query\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=5gb\n#PBS -l walltime=12:00:00\n#PBS -q iw-shared-6\n#PBS -j oe\n#PBS -o Split_" + query_base + ".out\n#PBS -e Split_" + query_base + ".err\n\n")
		split_pbs.write("\nmodule load intel\ncd "+ output_fullpath + "\n")
		split_pbs.write("\npython Split.py\n")
		split_pbs.close()

	SplitQuery()


	if os.path.exists(query_dir + "/" + query_base + "_Split"):
		if force == True:
			try:
				os.removedirs(query_dir + "/" + query_base + "_Split")
				os.makedirs(query_dir + "/" + query_base + "_Split")
			except OSError:
				pass
			os.system("qsub " + output_fullpath + "/Split.pbs")
			for i in range(1000):
				if os.path.exists(rundir + "/Split_" + query_base + ".out"):
					break
				else:
					print("Spliting Query")
					time.sleep(30)
					continue
			os.system("mv " + rundir + "/Split_" + query_base + ".out " + output_dir + "/" + output_name + "/")
	else:
		os.makedirs(query_dir + "/" + query_base + "_Split")
		os.system("qsub " + output_fullpath + "/Split.pbs")
		for i in range(1000):
			if os.path.exists(rundir + "/Split_" + query_base + ".out"):
				break
			else:
				print("Spliting Query")
				time.sleep(30)
				continue
		os.system("mv " + rundir + "/Split_" + query_base + ".out " + output_dir + "/" + output_name + "/")
	step += 1


"""------------------4.0 Run Blast------------------------"""

if step == 3:
	if not os.path.exists(output_fullpath):
		os.makedirs(output_fullpath)
	print("Running Blast: " + query_name)
	queries = glob.glob(query_dir + "/" + query_base + "_Split/" + query_base + "_*.fasta")
	#queries = glob.glob(output_dir + "/" + output_name + "/query_*.fasta")
	lon = len(queries)
	i = 0
	while i < lon:
		query_i = queries[i]
		blast_pbs = open(output_fullpath + "/Blast-" + query_base + "_" + str(i+1) + ".pbs", "w")
		blast_pbs.write("#PBS -N Blast-" + query_base + "_" + str(i+1)+"\n#PBS -l nodes=1:ppn=2\n#PBS -l mem=10gb\n#PBS -l walltime=12:00:00\n#PBS -q iw-shared-6\n#PBS -j oe\n#PBS -o Blast-" + query_base + "_" + str(i+1)+".out\n#PBS -e Blast-" + query_base + "_" + str(i+1)+".err\n\n")
		blast_pbs.write("\nmodule load boost/1.53.0;\nmodule load python/2.7;\nmodule load ncbi_blast/2.2.29;\n")
		blast_pbs.write("\ncd "+ output_fullpath +"\n")
		blast_pbs.write("\npwd\n")
		blast_pbs.write("\n" + program + " -db " + database_dir + "/" + database_base + " -query " + query_i + " -out " + output_fullpath + "/" + query_base + "_" + str(i+1) + ".fasta.blast -num_threads 2 " + ''.join(blast_opt) + "\n")
		blast_pbs.close()
		os.chdir(output_fullpath)
		os.system("qsub "+ output_fullpath + "/Blast-" + query_base + "_" + str(i+1) +".pbs")
		time.sleep(1)
		i += 1

print("Jobs Submitted. Check outputs in a while!")
