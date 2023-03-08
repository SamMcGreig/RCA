import os
import sys
import subprocess
import argparse
import datetime
import statistics
import math

from pathlib import Path
from Bio import SeqIO

def main():

	# Script version
	rca_version = 1

	# Input arguments
	options = parse_arguments()

	### Nanopore Pipeline
	if(sys.argv[1] == "nanopore"):

		# Directory list to create
		dirs = ["Ashure",
				"isONclust",
				"RCA_fastq",
				"Final_consensus"
			]

		if options.create_dirs == "Y":
			create_dirs(options.output, dirs)

		if options.ashure == "Y":
			run_ashure(options.input, options.output, options.ashure_path, options.ictv_db)
			extract_fastq(f"{options.output}/Ashure/", f"{options.output}/RCA_fastq/")

		if options.isonclust == "Y":
			run_isonclust(f"{options.output}/RCA_fastq/", f"{options.output}/isONclust/", options.threads, options.quality_score)

		if options.medaka == "Y":
			get_rep_seq(f"{options.output}/isONclust/")
			run_medaka(f"{options.output}/isONclust/", f"{options.output}/Final_consensus/", options.threads)

		if options.doc_env == "Y":
			document_env("RCA_nanopore", rca_version, options.output, options)

	### Illumina Pipeline
	if(sys.argv[1] == "illumina"):

		# Directory list to create
		dirs = ["Bbduk",
				"Spades",
				"Virus_detection/1.Putative_contigs",
				"Virus_detection/2.Mummer",
				"Virus_detection/3.Final_sequences"
			]

		if options.create_dirs == "Y":
			create_dirs(options.output, dirs)

		if options.bbduk == "Y":
			run_bbduk(options.input, f"{options.output}/Bbduk/", options.min_len, options.adapters, options.min_quality)

		if options.spades == "Y":
			run_spades(f"{options.output}/Bbduk/", f"{options.output}/Spades/", options.threads, options.memory)

		if options.virus_detection == "Y":
			extract_viruses(f"{options.output}/Spades/", f"{options.output}/Virus_detection/1.Putative_contigs/", 700)
			run_mummer(f"{options.output}/Virus_detection/1.Putative_contigs/", f"{options.output}/Virus_detection/2.Mummer/", f"{options.output}/Virus_detection/3.Final_sequences/")

		if options.doc_env == "Y":
			document_env("RCA_illumina", rca_version, options.output, options)

		

################################################################################
def create_dirs(output_dir, dirs):

	# Create the directories if they do not exist
	for directory in dirs:
		Path(output_dir, directory).mkdir(exist_ok = True, parents = True)

	print("Directory creation completed.")

################################################################################
def run_ashure(input_dir, output_dir, ashure, database):

	# Run Ashure and move the output files as required
	for file in os.listdir(input_dir):
		sample_name = file.split(".")[0]
		subprocess.call(f"python {ashure} run -fq {Path(input_dir, file)} -db {database} -o1 {input_dir}/Ashure/{sample_name}_consensus.csv -r fgs,msa", shell = True)
		os.rename("frags", f"{output_dir}/Ashure/{sample_name}_frags/")
		os.rename("msa", f"{output_dir}/Ashure/{sample_name}_msa/")
		os.rename("ashure.log", f"{output_dir}/Ashure/{sample_name}.log")

################################################################################
def extract_fastq(input_dir, output_dir):

	# Extract the fastq files and zip them
	for file in os.listdir(input_dir):
		if file.endswith("_msa"):
			sample_name = file.split("_")[0]
			subprocess.call(f"cat {input_dir}/{file}/*.fq > {output_dir}/{sample_name}.fastq", shell = True)

################################################################################
def run_isonclust(input_dir, output_dir, threads, quality_score):

	for file in os.listdir(input_dir):
		sample_name = file.split(".")[0]
		subprocess.call(f"isONclust --fastq {input_dir}/{file} --t {threads} --q {quality_score} --outfolder {output_dir}/{sample_name}/", shell = True)

		# For each fasta file, work out what 10% of the reads is at set as a threshold
		min_cluster_reads = 0
		with open(f"{output_dir}/{sample_name}/final_clusters.tsv") as clusters_tsv:
			largest_cluster_count = 0
			for line in clusters_tsv:
				line = line.split("\t")[0]
				if line == "0":
					largest_cluster_count += 1
				else:
					min_cluster_reads = int(round(largest_cluster_count * 0.1))
					break

		subprocess.call(f"isONclust write_fastq --clusters {output_dir}/{sample_name}/final_clusters.tsv --N {min_cluster_reads} --fastq {input_dir}/{file} --outfolder {output_dir}/{sample_name}/Clusters_fastq", shell = True)

################################################################################
def get_rep_seq(input_dir):

	for file in os.listdir(input_dir):
		sample_name = file.split(".")[0]
		create_dirs(f"{input_dir}/{sample_name}/", ["Rep_seq"])

		# Read through each fasta file
		for cluster in os.listdir(f"{input_dir}/{sample_name}/Clusters_fastq/"):
			kmer_size = 4
			for record in SeqIO.parse(f"{input_dir}/{sample_name}/Clusters_fastq/{cluster}", "fastq"):
				kmer_dict = {}
				n_kmers = len(record.seq) - kmer_size + 1

				for i in range(n_kmers):
					kmer = record.seq[i:i + kmer_size]
					if str(kmer) in kmer_dict:
						kmer_dict[str(kmer)] += 1
					else:
						kmer_dict.setdefault(str(kmer), 1)

				# Calculate Kmers
				repseq_check = "Y"
				mean_kmer = round(statistics.mean(kmer_dict.values()))
				for k in kmer_dict:
					if kmer_dict[k] > (mean_kmer * 10):
						# print(f"FAIL:	{sample_name}	{cluster}")
						# print(kmer_dict.values())
						repseq_check = "N"
						break

				if repseq_check == "Y":
					# print(f"PASS:	{sample_name}	{cluster}")
					# print(record.description)
					# print(record.seq)

					with open(f"{input_dir}/{sample_name}/Rep_seq/{cluster.split('.')[0]}.fasta", "w") as rep_seq_out:
						SeqIO.write(record, rep_seq_out, "fasta")

					break

################################################################################
def run_medaka(input_dir, output_dir, threads):

	for file in os.listdir(input_dir):
		sample_name = file.split(".")[0]
		create_dirs(f"{input_dir}/{sample_name}/", ["Medaka"])

		# Generate consensus sequences
		for repseq in os.listdir(f"{input_dir}/{sample_name}/Rep_seq/"):
			draft_seq = f"{input_dir}/{sample_name}/Rep_seq/{repseq}"
			cluster_id = repseq.split(".")[0]
			input_cluster = f"{input_dir}/{sample_name}/Clusters_fastq/{cluster_id}.fastq"
			medaka_output = f"{input_dir}/{sample_name}/Medaka/{cluster_id}/"

			subprocess.call(f"medaka_consensus -i {input_cluster} -d {draft_seq} -t {threads} -o {medaka_output}", shell = True)

		# Rename and relocate consensus sequences
		with open(f"{output_dir}/{sample_name}.fasta", "w") as consensus_out:
			for consensus in sorted(os.listdir(f"{input_dir}/{sample_name}/Medaka/")):
				for record in SeqIO.parse(f"{input_dir}/{sample_name}/Medaka/{consensus}/consensus.fasta", "fasta"):
					consensus_out.write(f">cluster_{consensus}_consensus\n")
					consensus_out.write(f"{record.seq}\n")

################################################################################
def run_bbduk(input_dir, output_dir, min_len, adapters, min_quality):

	for file in os.listdir(input_dir):
		if("_R1" in file):
			file_R1 = file.split(".")[0]
			file_R2 = file.replace("R1", "R2").split(".")[0]
			infile_R1 = f"{input_dir}/{file_R1}.fastq.gz"
			infile_R2 = f"{input_dir}/{file_R2}.fastq.gz"
			outfile_R1 = f"{output_dir}/{file_R1.replace('_L001_R1_001', '_R1')}_trimmed.fastq.gz"
			outfile_R2 = f"{output_dir}/{file_R2.replace('_L001_R2_001', '_R2')}_trimmed.fastq.gz"
			
			subprocess.call(f"bbduk.sh in1={infile_R1} in2={infile_R2} out1={outfile_R1} out2={outfile_R2} minlen={min_len} ktrim=r k=23 mink=11 hdist=1 ref={adapters} qtrim=r trimq={min_quality} ziplevel=5", shell = True)

################################################################################
def run_spades(input_dir, output_dir, threads, memory):

	for file in os.listdir(input_dir):
		if("_R1" in file):
			sample_name = file.split("_")[0]
			file_R1 = file.split(".")[0]
			file_R2 = file.replace("R1", "R2").split(".")[0]
			infile_R1 = f"{input_dir}/{file_R1}.fastq.gz"
			infile_R2 = f"{input_dir}/{file_R2}.fastq.gz"
			outfile = f"{output_dir}/{sample_name}"

			subprocess.call(f"spades.py --meta -1 {infile_R1} -2 {infile_R2} -o {outfile} -t {threads} -m {memory}", shell = True)

################################################################################
def extract_viruses(input_dir, output_dir, min_virus_len):

	for file in os.listdir(input_dir):
		sample_name = file
		create_dirs(f"{output_dir}/", [sample_name])

		seq_dict = {}
		max_cov = 0
		for record in SeqIO.parse(f"{input_dir}/{sample_name}/contigs.fasta", "fasta"):
			seq_dict.setdefault(record.description, str(record.seq))
			if float(record.description.split("_")[5]) > max_cov:
				max_cov = float(record.description.split("_")[5])

		for seq in seq_dict:
			if float(seq.split("_")[5]) >= (max_cov * 0.1):
				if len(seq_dict[seq]) >= 700:
					with open(f"{output_dir}/{sample_name}/{seq}.fasta", "w") as out_fasta:
						out_fasta.write(f">{seq}\n")
						out_fasta.write(f"{seq_dict[seq]}\n")

################################################################################
def run_mummer(input_dir, mummer_output_dir, fasta_output_dir):

	for sample in os.listdir(input_dir):
		sample_name = sample.split(".")[0]
		create_dirs(f"{mummer_output_dir}/", [sample_name])
		create_dirs(f"{fasta_output_dir}/", [sample_name])

		for sequence in os.listdir(f"{input_dir}/{sample}"):
			sequence_name = sequence.split(".fasta")[0]
			input_file = f"{input_dir}/{sample}/{sequence}"
			mummer_output_file = f"{mummer_output_dir}/{sample}/{sequence_name}.txt"
			fasta_output_file = f"{fasta_output_dir}/{sample}/{sequence_name}.fasta"

			subprocess.call(f"mummer -maxmatch -F -L {input_file} {input_file} > {mummer_output_file} 2> log.txt", shell = True)

			final_mummer_dict = {}
			mummer_dict = {}
			sequential_ID = 1
			print(sample_name)
			with open(mummer_output_file) as mummer_out:
				for line in mummer_out:
					if line.startswith(">"):
						contig_id = line.split()[1]
						mummer_dict.setdefault(contig_id, {})
					else:
						contig_id = line.split()[0]
						ref_pos = int(line.split()[1])
						alt_pos = int(line.split()[2])

						if contig_id in mummer_dict:
							if ref_pos < alt_pos:
								mummer_dict[contig_id].setdefault(sequential_ID, [])
								mummer_dict[contig_id][sequential_ID].append(ref_pos)
								mummer_dict[contig_id][sequential_ID].append(alt_pos)
								sequential_ID += 1

				seq = [0, 0]
				for key in mummer_dict:
					for sub_key in mummer_dict[key]:
						if mummer_dict[key][sub_key][1] - mummer_dict[key][sub_key][0] > (seq[1] - seq[0]):
							seq = mummer_dict[key][sub_key]
					final_mummer_dict.setdefault(key, seq)

				# print(mummer_dict)
				# print(f"FINAL 	{seq}")
				print(final_mummer_dict)

				for record in SeqIO.parse(input_file, "fasta"):
					with open(fasta_output_file, "w") as fasta_out:
						if int(seq[0]) == 0:
							fasta_out.write(f">{record.description}\n")
							fasta_out.write(f"{str(record.seq)}\n")
						else:
							fasta_out.write(f">{record.description}\n")
							fasta_out.write(f"{str(record.seq)[seq[0]:seq[1]]}\n")





################################################################################
def document_env(script_name, script_version, output_dir, input_params):
	# Report the arguments used to run the program
	# Report the environemnt the program was run in

	print(f"Printing {script_name} pipeline version information")
	with open(f"{output_dir}/{script_name}_Pipeline_params.txt", "w") as log_output:
		log_output.write(f"{script_name} Pipeline Version: {script_version}\n")
		log_output.write(f"Datetime: {datetime.datetime.now()}\n")
		log_output.write(f"Parameters:\n")
		for arg in vars(input_params):
			log_output.write(f"{arg} {getattr(input_params, arg)}\n")
	subprocess.call(f"conda list > {output_dir}/{script_name}_env.txt", shell = True)

################################################################################
def parse_arguments():
	parser = argparse.ArgumentParser(prog = "RCA", description = "Runs the RCA pipeline.")
	subparsers = parser.add_subparsers(help = "sub-command help")
	nanopore = subparsers.add_parser("nanopore", help = "Runs the RCA nanopore pipeline.")
	illumina = subparsers.add_parser("illumina", help = "Runs the RCA illumina pipeline.")

	################################################################################
	### Nanopore Pipeline

	# Key arguments
	nanopore.add_argument("--input", help = "This is the location of the raw data directory.", required = True)
	nanopore.add_argument("--output", help = "This is where the output data will be generated.", required = True)
	nanopore.add_argument("--ictv_db", help = "This is the path to the ictv database fasta file.", default = "/data/bigbio_00/smcgreig/Blast_databases/ICTV/ICTV_viruses.fasta")
	nanopore.add_argument("--ashure_path", help = "This is the path to the ashure script.", default = "/home/smcgreig/Scripts/RCA/ashure/src/ashure.py")

	# Extra arguments, useful for if a specific job has failed and you don't want to start from scratch
	nanopore.add_argument("--create_dirs", help = "Creates the directory structure. Default Y.", choices = ["Y", "N"], default = "Y")
	nanopore.add_argument("--ashure", help = "Runs Ashure. Default Y.", choices = ["Y", "N"], default = "Y")
	nanopore.add_argument("--isonclust", help = "Runs isONclust. Default Y.", choices = ["Y", "N"], default = "Y")
	nanopore.add_argument("--medaka", help = "Runs medaka. Default Y.", choices = ["Y", "N"], default = "Y")
	nanopore.add_argument("--doc_env", help = "Creates environment and parameter details. Default Y.", choices = ["Y", "N"], default = "Y")

	# Tool specific parameters

	# isONclust
	nanopore.add_argument("--threads", help = "The number of threads to run isONclust and Medaka with. Default 120.", default = 120)
	nanopore.add_argument("--quality_score", help = "The quality score cutoff for isONclust. Default 7.", default = 7)

	################################################################################
	### Illumina Pipeline

	# Key arguments
	illumina.add_argument("--input", help = "This is the location of the raw data directory.", required = True)
	illumina.add_argument("--output", help = "This is where the output data will be generated.", required = True)

	illumina.add_argument("--create_dirs", help = "Creates the directory structure. Default Y.", choices = ["Y", "N"], default = "Y")
	illumina.add_argument("--bbduk", help = "Runs Bbduk. Default Y.", choices = ["Y", "N"], default = "Y")
	illumina.add_argument("--spades", help = "Runs isONclust. Default Y.", choices = ["Y", "N"], default = "Y")
	illumina.add_argument("--virus_detection", help = "Extract potential viruses. Default Y.", choices = ["Y", "N"], default = "Y")
	illumina.add_argument("--doc_env", help = "Creates environment and parameter details. Default Y.", choices = ["Y", "N"], default = "Y")

	# Tool specific parameters

	# Bbduk
	illumina.add_argument("--adapters", help = "bbduk adapter references.", default = "/home/smcgreig/Scripts/AnguaWGS/adapters.fa")
	illumina.add_argument("--min_len", help = "This is the mininum length of a read allowed by bbduk. Default 50.", default = 50)
	illumina.add_argument("--min_quality", help = "bbduk phred quality trim parameter. Default 10.", default = 10)

	# Spades
	illumina.add_argument("--threads", help = "The number of threads to run Spades with. Default 60.", default = 60)
	illumina.add_argument("--memory", help = "The amount of memory, in GB, to run Spades with. Default 500.", default = 500)

	# Virus detection
	illumina.add_argument("--hmm", help = "This is the path to the hmm file.", default = "/home/smcgreig/miniconda2/envs/RCA/Databases/nbc_hmms.hmm")


	return parser.parse_args()

################################################################################

if __name__ == '__main__':
	main()
