#!/usr/bin/env python
#
#	ACBC-Assemble
#
#	This script supervises assembly of paired-end reads.  It started as a fork of Matt Johnson's
#	HybPiper (https://github.com/mossmatters/HybPiper), which was licensed under the GNU General
#	Public License v3.0.  I have extensively modified it, using my own style, and adapted it to
#	facilitate analysis of the targeted-capture data generated in the Gavin Naylor lab (Florida
#	Shark Research Program).
#
#	Significant modifications to the original HybPiper code include:
#	  * Improved modularity/encapsulation; cleaned up code formatting and internal documentation.
#	  *	Code supporting the Velvet assembler has been removed.
#	  *	Support for unpaired reads has been eliminated.
#	  *	Major changes to the logic for driving the assemblies and iterating k-mer lists so as to
#		find contigs that overlap the exon reference sequences.
#
#	Command line options are passed to the other executables, which are run as separate commands
#	in order to take advantage of GNU Parallel.
#
#	David L. Swofford, 5 Feb 2019

from __future__ import print_function	# for python 2 vs 3 compatibility

import acbc
import argparse, os, sys
import importlib, shutil, glob
import subprocess
import logging

exonerate_genefilename = "exonerate_genelist.txt"
spades_genefilename = "spades_genelist.txt"

required_executables = [
	  "blastx"
	, "exonerate"
	, "parallel"
	, "makeblastdb"
	, "spades.py"
	, "bwa"
	, "samtools"
	]

required_python_packages = [
	  "Bio"
	]

required_other_scripts = [
	  "distribute_reads_to_targets.py"
	, "distribute_targets.py"
	, "exonerate_hits.py"
	]

def check_dependencies(run_dir):
	"""
	Checks for the presence of required executables and Python packages.
	"""
	ok = True
	print("\nChecking for executables...")
	for e in required_executables:
		e_loc = py_which(e)
		if e_loc:
			print("  {} found at {}".format(e, e_loc))
		else:
			print("  {} not found in your $PATH.".format(e))
			ok = False

	print("\nChecking for required Python packages...")
	for p in required_python_packages:
		try:
			i = importlib.import_module(p)
			print("  Package {} was successfully loaded.".format(p))
		except ImportError:
			print("  Package {} not found.".format(p))
			ok = False

	print("\nChecking for other pipeline scripts...")
	for s in required_other_scripts:
		print("  Script '{}' ".format(s), end="")
		if os.path.isfile(os.path.join(run_dir, s)):
			print("found.")
		else:
			print("not found.  Make sure it is in the same	directory as '{}'." \
			      .format(os.path.basename(__file__)))
			ok = False

	return ok

def blastx(readfiles, baitfile, evalue, basename, cpu=None, max_target_seqs=10):
	dna = set("ATCGN")
	if os.path.isfile(baitfile):
#		Quick detection of whether baitfile is DNA.
		with open(baitfile) as bf:
			header = bf.readline()
			seqline = bf.readline().rstrip().upper()
			if not set(seqline) - dna:
				print("\nError: Only ATCGN characters were found in first line of {}. You need " \
				      "a protein bait file for BLASTx!".format(baitfile))
				return None

		if os.path.isfile(os.path.split(baitfile)[0]+'.psq'):
			db_file = baitfile
		else:
			print("Making protein blastdb in current directory.")
			if os.path.split(baitfile)[0]:
				shutil.copy(baitfile, '.')
			db_file = os.path.split(baitfile)[1]
			makeblastdb_cmd = "makeblastdb -dbtype prot -in {}".format(db_file)
			print(makeblastdb_cmd)
			exitcode = subprocess.call(makeblastdb_cmd, shell=True)
			if exitcode:
				return None
	else:
		print(("Cannot find baitfile at: {}".format(baitfile)))
		return None

#	Remove previous blast results if they exist (because we will be appending)
	acbc.delete_file(basename+".blastx")

	for read_file in readfiles:

#		Piping commands for Fastq -> FASTA (curly braces must be doubled within a formatted string)
		pipe_cmd = "cat {} " \
		           "|  awk '{{if(NR % 4 == 1 || NR % 4 == 2) {{sub(/@/, \">\"); print; }} }}'" \
		           .format(read_file)
		blastx_command = "blastx -db {} -query - -evalue {} -outfmt 6 -max_target_seqs {}" \
		                 .format(db_file, evalue, max_target_seqs)
		if cpu:
			full_command = "time {} | parallel -j {} -k --block 200K --recstart '>' --pipe '{}' >> {}.blastx ".format(pipe_cmd, cpu, blastx_command, basename)
		else:
			full_command = "time {} | parallel -k --block 200K --recstart '>' --pipe '{}' >> {}.blastx ".format(pipe_cmd, blastx_command, basename)
		print(full_command)
		exitcode = subprocess.call(full_command, shell=True)
		if exitcode:
			#Concatenate the two blastfiles.

			return None

	return basename + '.blastx'

def distribute_blast(blastx_outputfile, readfiles, baitfile, run_dir, target=None, exclude=None):
	#NEED TO ADD SOMETHING ABOUT DIRECTORIES HERE.
	#print run_dir
	read_cmd = "time python {} {} {}".format(os.path.join(run_dir, "distribute_reads_to_targets.py"), blastx_outputfile, " ".join(readfiles))
	exitcode = subprocess.call(read_cmd, shell=True)
	if exitcode:
		print("ERROR: Something went wrong with distributing reads to gene directories.")
		return exitcode
	target_cmds = ["time python", os.path.join(run_dir, "distribute_targets.py"), baitfile, "--blastx", blastx_outputfile]
	if target:
		target_cmds.append("--target {}".format(target))
	if exclude:
		target_cmds.append("-- exclude {}".format(exclude))
	target_cmd = " ".join(target_cmds)
	exitcode = subprocess.call(target_cmd, shell=True)
	if exitcode:
		print("ERROR: Something went wrong distributing targets to gene directories.")
		return exitcode
	return None

def distribute_bwa(bamfile, readfiles, baitfile, run_dir, target=None, exclude=None):
	#NEED TO ADD SOMETHING ABOUT DIRECTORIES HERE.
	#print run_dir
	read_cmd = "time python {} {} {}".format(os.path.join(run_dir, "distribute_reads_to_targets_bwa.py"), bamfile, " ".join(readfiles))
	print(("[CMD] {}\n".format(read_cmd)))
	exitcode = subprocess.call(read_cmd, shell=True)

	if exitcode:
		print("ERROR: Something went wrong with distributing reads to gene directories.")
		return exitcode
	target_cmds = ["time python", os.path.join(run_dir, "distribute_targets.py"), baitfile, "--bam", bamfile]
	if target:
		target_cmds.append("--target {}".format(target))
	if exclude:
		target_cmds.append("--exclude {}".format(exclude))
	target_cmd = " ".join(target_cmds)
	print("[DISTRIBUTE]: {}".format(target_cmd))
	exitcode = subprocess.call(target_cmd, shell=True)
	if exitcode:
		print("ERROR: Something went wrong distributing targets to gene directories.")
		return exitcode
	return None

def make_basename(readfiles, prefix=None):
	"""Unless prefix is set, generate a directory based off the readfiles, using everything up to the first underscore.
		If prefix is set, generate the directory "prefix" and set basename to be the last component of the path.

		"""
	if prefix:
			if not os.path.exists(prefix):
				os.makedirs(prefix)

			prefixParentDir, prefix = os.path.split(prefix)
			if not prefix:
				# if prefix has a trailing /, prefixParentDir will have the / stripped and prefix will be empty.
				# so try again
				prefix = os.path.split(prefixParentDir)[1]


			return prefixParentDir, prefix

		## --prefix is not set on cmd line; write output to subdir in .
	basename = os.path.split(readfiles[0])[1].split('_')[0]

	if not os.path.exists(basename):
		os.makedirs(basename)
	return '.', basename

def spades(genes_to_assemble, run_dir, cov_cutoff=8, cpu=None, paired=True, kvals=None, timeout=None, suppress_rerun=False, rerun_only=False):
	"Run SPAdes on each gene separately using GNU parallel."""

	logger.debug("entering spades() with genes_to_assemble --> {}".format(genes_to_assemble))
	acbc.make_genelist_file_from_set(genes_to_assemble, spades_genefilename)

	acbc.delete_file("spades.log")
	acbc.delete_file("spades_redo.log")

	spades_runner_list = ["python", "{}/spades_runner.py".format(run_dir), spades_genefilename, "--cov_cutoff", str(cov_cutoff)]
	if cpu:
		spades_runner_list.append("--cpu")
		spades_runner_list.append(str(cpu))
	if not paired:
		spades_runner_list.append("--single")
	if timeout:
		spades_runner_list.append("--timeout")
		spades_runner_list.append("{}%".format(timeout))
	if kvals:
		spades_runner_list.append("--kvals")
		spades_runner_list.append("{}".format(", ".join(kvals)))
	if suppress_rerun:
		spades_runner_list.append("--suppress_rerun")
	if rerun_only:
		spades_runner_list.append("--redos_only")

	spades_runner_cmd = " ".join(spades_runner_list)
	print("\nCalling SPAdes with command: {}\n".format(spades_runner_cmd))
	exitcode = subprocess.call(spades_runner_cmd, shell=True)
	if exitcode:
		sys.stderr.write("WARNING: Something went wrong with the assemblies! Check for failed assemblies and re-run!\n")
		return None

	spades_failed = acbc.make_set_from_genelist_file("failed_spades.txt")
	spades_duds = acbc.make_set_from_genelist_file("spades_duds.txt")
	acbc.delete_file("failed_spades.txt")
#	acbc.delete_file("spades_duds.txt")
	logger.debug("spades: spades_failed = {}".format(spades_failed))
	logger.debug("spades: spades_duds = {}".format(spades_duds))

	return spades_failed, spades_duds

def exonerate(genes_to_exonerate, basename, run_dir, replace=True, cpu=None, thresh=55, depth_multiplier=0,
              length_pct=100, timeout=None):
	if replace:
		for g in genes_to_exonerate:
			if os.path.isdir(os.path.join(g, basename)):
				shutil.rmtree(os.path.join(g, basename))
	if len(genes_to_exonerate) == 0:  ### can still happen???
		print(("nError: No genes recovered for {}.".format(basename)))
		return 1

	print(("\nRunning Exonerate to generate sequences for {} genes".format(len(genes_to_exonerate))))

	parallel_cmd_list = ["time parallel", "--eta"]
	if cpu:
		parallel_cmd_list.append("-j {}".format(cpu))
	if timeout:
		parallel_cmd_list.append("--timeout {}%".format(timeout))

#	Build Exonerate command.
	acbc.make_genelist_file_from_set(genes_to_exonerate, exonerate_genefilename)
	exonerate_cmd_list = [
		  "python"
		, "{}/exonerate_hits.py".format(run_dir)
		, "{}/{}_baits.fasta"
		, "{{}}/{{}}_{}".format("contigs.fasta")
		, "--prefix {{}}/{}".format(basename)
		, "-t {}".format(thresh)
		, "--depth_multiplier {}".format(depth_multiplier)
		, "--length_pct {}".format(length_pct)
		, "::::"
		, exonerate_genefilename
		, "> exonerate_success.txt"
		]

	exonerate_cmd = " ".join(parallel_cmd_list) + " " + " ".join(exonerate_cmd_list)
	logger.info("calling '{}'".format(exonerate_cmd))
	exitcode = subprocess.call(exonerate_cmd, shell=True)
	if exitcode:
		print("\nError: Something went wrong with Exonerate.")
		return exitcode, None

	exonerate_success = acbc.make_set_from_genelist_file("exonerate_success.txt")
	acbc.delete_file("exonerate_success.txt")
	logger.debug("exonerate: exonerate_success --> {}".format(exonerate_success))
	return 0, exonerate_success

def bwa(readfiles, baitfile, basename, cpu):
	"""Conduct BWA search of reads against the baitfile.
	Returns an error if the second line of the baitfile contains characters other than ACTGN"""
	dna = set("ATCGN")
	if os.path.isfile(baitfile):
		#Quick detection of whether baitfile is DNA.
		with open(baitfile) as bf:
			header = bf.readline()
			seqline = bf.readline().rstrip().upper()
			if set(seqline) - dna:
				print("ERROR: characters other than ACTGN found in first line. You need a nucleotide bait file for BWA!")
				return None

		if os.path.isfile(os.path.split(baitfile)[0]+'.amb'):
			db_file = baitfile
		else:
			print("Making nucleotide bwa index in current directory.")
			baitfileDir = os.path.split(baitfile)[0]
			if baitfileDir:
				if os.path.realpath(baitfileDir) != os.path.realpath('.'):
					shutil.copy(baitfile, '.')
			db_file = os.path.split(baitfile)[1]
			make_bwa_index_cmd = "bwa index {}".format(db_file)
			print(("[CMD]: {}".format(make_bwa_index_cmd)))
			exitcode = subprocess.call(make_bwa_index_cmd, shell=True)
			if exitcode:
				return None
	else:
		print(("ERROR: Cannot find baitfile at: {}".format(baitfile)))
		return None

	if not cpu:
		import multiprocessing
		cpu = multiprocessing.cpu_count()

	if len(readfiles) < 3:
		bwa_fastq = " ".join(readfiles)
	else:
		bwa_fastq = readfiles

	bwa_commands = ["time bwa mem", "-t", str(cpu), db_file, bwa_fastq, " | samtools view -h -b -S - > "]
	bwa_commands.append(basename+".bam")
	full_command = " ".join(bwa_commands)
	print(("[CMD]: {}".format(full_command)))
	exitcode = subprocess.call(full_command, shell=True)
	if exitcode:
		return None

	return basename + '.bam'

def main():

	helptext = """
	ACBC-Assemble

	This script supervises assembly of paired-end reads.  Command line options are passed to the
	other executables, which are run as separate commands in order to take advantage of GNU
	Parallel.

	It can optionally check whether required dependencies are installed and that the other scripts
	called by this one are found in the same directory as this one (see --check-depend).

	Unless --prefix is set, output will be put within a directory named after your read files.
	"""
	parser = argparse.ArgumentParser(description=helptext, formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument("--check-depend",dest='check_depend',help="Check for dependencies (executables and Python packages) and exit. May not work at all on Windows.",action='store_true')
	parser.add_argument("--bwa",dest="bwa",action='store_true',help="Use BWA to search reads for hits to target. Requires BWA and a bait file that is nucleotides!",default=False)
	parser.add_argument("--no-blast",dest="blast",action="store_false",help="Do not run the blast step. Downstream steps will still depend on the *_all.blastx file. \nUseful for re-runnning assembly/exonerate steps with different options.")
	parser.add_argument("--no-distribute",dest="distribute",action="store_false",help="Do not distribute the reads and bait sequences to sub-directories.")
	parser.add_argument("--no-velvet",dest="velvet",action="store_false",help="Do not run the velvet stages (velveth and velvetg)")
	parser.add_argument("--no-exonerate",dest="exonerate",action="store_false",help="Do not run the Exonerate step, which assembles full length CDS regions and proteins from each gene")
	parser.add_argument("--force_exonerate_hit",dest="force_exonerate_hit",action="store_true",help="Iterate SPAdes assembly if necessary, reducing largest k-mer value until Exonerate finds a hit to the target sequence")
	parser.add_argument("--no-assemble",dest="assemble",action="store_false",help="Skip the SPAdes assembly stage.")
	parser.add_argument("--keep_failed_assemblies",action="store_true",help="keep directories for failed assemblies; default is to remove them", default = False)

	parser.add_argument('-r',"--readfiles",nargs='+',help="One or more read files to start the pipeline. If exactly two are specified, will assume it is paired Illumina reads.",default=[])
	parser.add_argument('-b','--baitfile',help="FASTA file containing bait sequences for each gene. If there are multiple baits for a gene, the id must be of the form: >Taxon-geneName",default=None)

	parser.add_argument('--cpu',type=int,default=0,help="Limit the number of CPUs. Default is to use all cores available.")
	parser.add_argument('--evalue',type=float,default=1e-10,help="e-value threshold for blastx hits, default: %(default)s")
	parser.add_argument('--max_target_seqs',type=int,default=10,help='Max target seqs to save in blast search, default: %(default)s')
	parser.add_argument('--cov_cutoff',type=int,default=8,help="Coverage cutoff for velvetg. default: %(default)s")
	parser.add_argument('--ins_length',type=int,default=200,help="Insert length for velvetg. default: %(default)s")
	parser.add_argument("--kvals",nargs='+',help="Values of k for velvet assemblies. Velvet needs to be compiled to handle larger k-values! Default auto-dectection by SPAdes.",default=None)
	parser.add_argument("--thresh",type=int,help="Percent Identity Threshold for stitching together exonerate results. Default is 55, but increase this if you are worried about contaminant sequences.",default=65)
	parser.add_argument("--length_pct",help="Include an exonerate hit if it is at least as long as X percentage of the reference protein length. Default = 90%%",default=90,type=int)
	parser.add_argument("--depth_multiplier",help="Accept any full-length exonerate hit if it has a coverage depth X times the next best hit. Set to zero to not use depth. Default = 10",default=10,type=int)

	parser.add_argument('--prefix',help="Directory name for pipeline output, default is to use the FASTQ file name.",default=None)
	parser.add_argument("--timeout",help="Use GNU Parallel to kill long-running processes if they take longer than X percent of average.",default=0)

	parser.add_argument("--target",help="Use this target to align sequences for each gene. Other targets for that gene will be used only for read sorting. Can be a tab-delimited file (one gene per line) or a single sequence name",default=None)
	parser.add_argument("--exclude",help="Do not use any sequence with the specified string as a target sequence for exonerate. The sequence will be used for read sorting.",default=None)
	parser.set_defaults(check_depend=False,blast=True,distribute=True,velvet=False,assemble=True,exonerate=True)
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)
	args = parser.parse_args()

	global logger
	logger = acbc.set_logger("acbc", debug=False)

	run_dir = os.path.realpath(os.path.split(sys.argv[0])[0])
	print(("HybPiper was called with these arguments:\n{}".format(" ".join(sys.argv))))

	if args.check_depend:
		if check_dependencies(run_dir):
			print("\nAll dependencies seem to be correctly installed.")
		else:
			print("\nError: One or more dependencies was not found.")
		return

	if args.baitfile:
		baitfile = os.path.abspath(args.baitfile)
	else:
		parser.print_help()
		return
	if len(args.readfiles) < 1:
		print("\nError: Please specify readfiles with -r")
		return
	elif len(args.readfiles) > 2:
		print("\nError: Please specify exactly two (paired) read files.")
		return
	if not args.baitfile:
		print("\nError: Please specify a FASTA file containing target sequences.")
		return

	readfiles = [os.path.abspath(x) for x in args.readfiles]
	print("\nBeginning assembly")
	print("  readfiles     = {}".format("\n                  ".join(readfiles)))

#	Generate directory
	basedir, basename = make_basename(args.readfiles, prefix=args.prefix)
	os.chdir(os.path.join(basedir, basename))


	if 1:	### need a start clean arg, and load from prev files if not starting clean
		gene_hits = set()
		spades_successes = set()
		exonerate_hits = set()

#	BWA...
	if args.bwa:
		if args.blast:
			args.blast=False
			bamfile = bwa(readfiles, baitfile, basename, cpu=args.cpu)
			if not bamfile:
				print("\nError: Something went wrong with the BWA step, exiting!")
				return
		else:
			bamfile = basename + ".bam"

#	BLAST...
	if args.blast:
		blastx_outputfile = blastx(readfiles, baitfile, args.evalue, basename, cpu=args.cpu, max_target_seqs=args.max_target_seqs)
		if not blastx_outputfile:
			print("\nError: Something is wrong with the Blastx step, exiting!")
			return
	else:
		blastx_outputfile = basename + ".blastx"

#	Distribute reads...
	if args.distribute:
		pre_existing_fastas = glob.glob("*/*_interleaved.fasta")
		### provide option to remove existing assembled files here (--always-assemble ?) in which case remove _contigs file and directories
		for fn in pre_existing_fastas:
			os.remove(fn)
		if args.bwa:
			exitcode = distribute_bwa(bamfile, readfiles, baitfile, run_dir, args.target, args.exclude)
		else:
			exitcode = distribute_blast(blastx_outputfile, readfiles, baitfile, run_dir, args.target, args.exclude)
		if exitcode:
			sys.exit(1)

	genes_with_reads = set(x for x in os.listdir(".") if os.path.isfile(os.path.join(x, x + "_interleaved.fasta")))
	logger.debug("genes_with_reads --> {}".format(genes_with_reads))
	if not genes_with_reads:
		print("\nError: No genes with BLAST hits were found.")
		return

	gene_hits |= genes_with_reads

#	Initialize 'assembled_genes' to genes that have a ""*_contigs.fasta" file.
	assembled_genes = set(x for x in os.listdir(".") if os.path.isfile(os.path.join(x, x + "_contigs.fasta")))
	exonerate_successes = set()		### get from prev file

	assembled_here = set()
	failures_here = set()
	genes_to_assemble = genes_with_reads - assembled_genes
	logger.debug("genes_to_assemble --> {}".format(genes_to_assemble))
	if not genes_to_assemble:
		print("\nAll genes for which reads are available have already been assembled. " \
		      "Remove their directories to re-assemble.")
		return

#	Start fresh, removing any files left over from earlier runs.
	### OR NOT???
	acbc.delete_file("spades_genelist.txt")
	acbc.delete_file("spades_duds.txt")

	print("\nAssembling {} genes".format(len(genes_with_reads)))

	iter = 0
	done = False
	while genes_to_assemble:
		doing_rerun = (iter != 0)
		if doing_rerun:
			acbc.make_genelist_file_from_set(spades_failed, "failed_spades.txt")
		spades_failed, spades_duds = spades(genes_to_assemble, run_dir, cov_cutoff=args.cov_cutoff, cpu=args.cpu,
		                                    kvals=args.kvals, timeout=args.timeout, suppress_rerun=True,
		                                    rerun_only=doing_rerun)
		spades_successes |= genes_to_assemble - (spades_failed | spades_duds)
		logger.debug("spades_failed from spades --> {}".format(spades_failed))
		logger.debug("spades_duds from spades   --> {}".format(spades_duds))
		logger.debug("spades_successes          --> {}".format(spades_successes))

		### use force_exonerate_hit

		genes_to_assemble -= spades_duds
		exonerate_failed = set()
		if args.exonerate:
			genes_to_exonerate = genes_to_assemble - spades_failed
			logger.debug("genes_to_exonerate --> {}".format(genes_to_exonerate))
			if genes_to_exonerate:
				exitcode, successes = exonerate(genes_to_exonerate, basename, run_dir, cpu=args.cpu,
				                                thresh=args.thresh, length_pct=args.length_pct,
				                                depth_multiplier=args.depth_multiplier, timeout=args.timeout)
				if exitcode:
					return

				exonerate_successes |= successes

				logger.debug("back from exonerate: successes   = {}".format(successes))
				exonerate_failed = genes_to_exonerate - successes
				assembled_here |= successes
				genes_to_assemble -= successes
				if args.force_exonerate_hit:
					genes_to_assemble |= exonerate_failed
				logger.debug("                     exonerate_failed    = {}".format(exonerate_failed))
				logger.debug("                     assembled_here    --> {}".format(assembled_here))
				logger.debug("                     assembled_genes   --> {}".format(assembled_genes))
				logger.debug("                     genes_to_assemble --> {}".format(genes_to_assemble))

		iter += 1
		assembled_genes |= successes
		logger.debug("assembled_here for next iteration:")
		logger.debug("    assembled_here    = {}".format(assembled_here))
		logger.debug("    assembled_genes   = {}".format(assembled_genes))
		logger.debug("    genes_to_assemble = {}\n".format(genes_to_assemble))

#	Remove directories for unsuccessful assemblies.
	failed_genes = genes_with_reads - assembled_genes
	if not args.keep_failed_assemblies:
		for g in failed_genes:
			shutil.rmtree(g)

	print("\nNew assemblies generated in this run for {} genes".format(len(assembled_here)))
	print("Total number of genes now successfully assembled = {}".format(len(assembled_genes)))
	print("Number of genes for which assembly was unsuccessful in this run = {}".format(len(failed_genes)))

	acbc.make_genelist_file_from_set(assembled_genes, "genes_with_seqs.txt")
	acbc.make_genelist_file_from_set(exonerate_successes, "exonerate_genelist.txt")
	acbc.make_genelist_file_from_set(failed_genes, "failed_genes.txt")

	paralog_warnings = [x for x in os.listdir(".") if os.path.isfile(os.path.join(x, basename, "paralog_warning.txt"))]
	acbc.make_genelist_file_from_set(paralog_warnings, "genes_with_paralog_warnings.txt")
	print("Number of genes for which potential paralogs were detected = {}\n".format(len(paralog_warnings)))

	acbc.delete_file("spades_genelist.txt")
#	acbc.delete_file("spades_duds.txt")

if __name__ == "__main__":
	main()
