#!/usr/bin/env python

from __future__ import print_function	# for python 2 vs 3 compatibility

import acbc
import argparse, os, sys, shutil, subprocess
import logging

def make_spades_cmd(genelist,cov_cutoff=8,cpu=None,paired=True,kvals=None,redo=False,timeout=None,unpaired=False):

    if kvals:
        kvals = ",".join(kvals)

    parallel_cmd_list = ["time","parallel","--eta"]
    if cpu:
        parallel_cmd_list.append("-j {}".format(cpu))
    if timeout:
        parallel_cmd_list.append("--timeout {}%".format(timeout))

    spades_cmd_list = ["spades.py --only-assembler --threads 1 --cov-cutoff",str(cov_cutoff)]
    if kvals:
        spades_cmd_list.append("-k {}".format(kvals))
    if unpaired:
        spades_cmd_list.append("-s {}/{}_unpaired.fasta")
    if paired:
        spades_cmd_list.append("--12 {}/{}_interleaved.fasta")
    else:
        spades_cmd_list.append("-s {}/{}_unpaired.fasta")

    spades_cmd_list.append("-o {{}}/{{}}_spades :::: {} > spades.log".format(genelist))

    spades_cmd = " ".join(parallel_cmd_list) + " " + " ".join(spades_cmd_list)
    return spades_cmd


def spades_initial(genelist_file,cov_cutoff=8,cpu=None,paired=True,kvals=None,timeout=None,unpaired=False):
    "Run SPAdes on each gene separately using GNU parallel."""
    if os.path.isfile("spades.log"):
        os.remove("spades.log")

    genes = [x.rstrip() for x in open(genelist_file)]
    #print paired
    spades_cmd = make_spades_cmd(genelist_file, cov_cutoff, cpu, paired=paired, kvals=kvals, unpaired=unpaired)

    sys.stderr.write("Running SPAdes on {} genes\n".format(len(genes)))
    sys.stderr.write(spades_cmd + "\n")
    exitcode = subprocess.call(spades_cmd,shell=True)

    if exitcode:
        sys.stderr.write("\nError: One or more genes had an error with SPAdes assembly. This may be due to low coverage. No contigs found for the following genes:\n")

    spades_failed = []

    for gene in genes:
        gene_failed = False
        if os.path.isfile("{}/{}_spades/contigs.fasta".format(gene,gene)):
            contig_file_size = os.stat("{}/{}_spades/contigs.fasta".format(gene,gene)).st_size
            if  contig_file_size> 0:
                shutil.copy("{}/{}_spades/contigs.fasta".format(gene,gene),"{}/{}_contigs.fasta".format(gene,gene))
            else:
                gene_failed = True
        else:
            gene_failed = True

        if gene_failed:
            sys.stderr.write("{}\n".format(gene))
            spades_failed.append(gene)
    return spades_failed

def rerun_spades(genelist_file,cov_cutoff=8,cpu=None, paired = True):
    genes = [x.rstrip() for x in open(genelist_file)]

    redo_cmds_file = open("redo_spades_commands.txt",'w')

#   The 'spades_failed' set contains genes for which a good assembly has (still) not been obtained for the current
#   k-mer set.  The 'spades_duds' set contains genes for which the assembly failed for even the shortest k-mer.
    spades_failed = []
    spades_duds = []

    genes_redos = []
    restart_ks = []
    for gene in genes:
        all_kmers = [int(x[1:]) for x in os.listdir(os.path.join(gene,"{}_spades".format(gene))) if x.startswith("K")]
        all_kmers.sort()

        if len(all_kmers) < 2:
            sys.stderr.write("WARNING: All Kmers failed for {}!\n".format(gene))
            spades_duds.append(gene)
            continue
        else:
            genes_redos.append(gene)
        redo_kmers = [str(x) for x in all_kmers[:-1]]
        logger.debug("redo_kmers --> {}\n".format(redo_kmers))
        restart_k = "k{}".format(redo_kmers[-1])

#       Remove data for longest k-mer so that restart will start from shorter maximum k-mer length:
        last = max(all_kmers)
        dir_to_remove = os.path.join(gene,"{}_spades".format(gene), "K"+str(last))
        logger.debug("rerun_spades deleting dir {}\n".format(dir_to_remove))
        shutil.rmtree(dir_to_remove)

        kvals = ",".join(redo_kmers)
        spades_cmd = "spades.py --restart-from {} -k {} --cov-cutoff {} -o {}/{}_spades".format(restart_k, kvals, cov_cutoff, gene, gene)
        redo_cmds_file.write(spades_cmd + "\n")

    redo_cmds_file.close()
    if cpu:
        redo_spades_cmd = "parallel -j {} --eta --timeout 400% :::: redo_spades_commands.txt > spades_redo.log".format(cpu)
    else:
        redo_spades_cmd = "parallel --eta --timeout 400% :::: redo_spades_commands.txt > spades_redo.log"

    sys.stderr.write("Re-running SPAdes for {} genes\n".format(len(genes_redos)))
    sys.stderr.write(redo_spades_cmd+"\n")
    exitcode = subprocess.call(redo_spades_cmd,shell=True)
    logger.debug("redo_spades_cmd returned exitcode = {}".format(exitcode))
    if exitcode:
        sys.stderr.write("ERROR: One or more genes had an error with SPAdes assembly. This may be due to low coverage. No contigs found for the following genes:\n")

    for gene in genes_redos:
        logger.debug("checking gene {} for success:\n".format(gene))
        gene_failed = False
        if os.path.isfile("{}/{}_spades/contigs.fasta".format(gene, gene)):
            gene_file = "{}/{}_spades/contigs.fasta".format(gene, gene)
            logger.debug("    file {} is present, size = {}:".format(gene_file, os.stat(gene_file).st_size))
            if os.stat("{}/{}_spades/contigs.fasta".format(gene, gene)).st_size > 0:
                shutil.copy("{}/{}_spades/contigs.fasta".format(gene, gene),"{}/{}_contigs.fasta".format(gene, gene))
            else:
                gene_failed = True
        else:
            gene_failed = True

        if gene_failed:
            spades_failed.append(gene)

    acbc.delete_file("redo_spades_commands.txt")
    return spades_failed, spades_duds

def main():

    helptext = """Runs the SPAdes assembler, optionally with re-dos if any of the k-mers are
                  unsuccessful (these re-runs are attempted by removing the largest k-mer and
                  re-running spades)."""
    parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('genelist', help="Text file containing the name of each gene to conduct SPAdes assembly. One gene per line, should correspond to directories within the current directory.")### rename genefile?
    parser.add_argument('--cpu', type=int, default=0, help="Limit the number of CPUs. Default is to use all cores available.")
    parser.add_argument('--cov_cutoff', type=int, default=8, help="Coverage cutoff for SPAdes. default: %(default)s")
    parser.add_argument("--kvals", nargs='+', help="Values of k for SPAdes assemblies. Default is to use SPAdes auto detection based on read lengths (recommended).", default=None)
    parser.add_argument("--redos_only", action="store_true", default=False, help="Continue from previously assembled SPAdes assemblies and only conduct redos from failed_spades.txt")
    parser.add_argument("--single", help="Reads are single end. Default is paired end.", action='store_true', default=False)
    parser.add_argument("--timeout", help="Use GNU Parallel to kill processes that take longer than X times the average.", default=0)
    parser.add_argument("--unpaired", help="For assembly with both paired (interleaved) and unpaired reads", action="store_true", default=False) ####
    parser.add_argument("--suppress_rerun", help="Don't rerun spades with reduced k-mer set after failure", action="store_true", default=False)
    args = parser.parse_args()

    global logger
    logger = acbc.set_logger("spades_runner", debug=False)

    is_paired = not args.single

    failed_spadesfile_name = "failed_spades.txt"

    if os.path.isfile(failed_spadesfile_name) and args.redos_only:
        logger.debug("calling rerun_spades")###
        spades_failed, spades_duds = rerun_spades(failed_spadesfile_name, cpu=args.cpu, paired=is_paired)
        logger.debug("rerun_spades returned spades_failed = {}:\n".format(spades_failed))
    else:
        logger.debug("calling spades_initial with args.genelist = {}:\n".format(args.genelist))
        spades_failed = spades_initial(args.genelist, cov_cutoff=args.cov_cutoff, cpu=args.cpu, kvals=args.kvals,
                                       paired=is_paired, timeout=args.timeout, unpaired=args.unpaired)
        logger.debug("spades_initial returned spades_failed = {}:\n".format(spades_failed))

    if len(spades_failed) == 0:
        acbc.delete_file(failed_spadesfile_name)
    else:
        acbc.make_genelist_file_from_set(spades_failed, failed_spadesfile_name)

if __name__ == "__main__":
    main()
