import os
import sys
import re
import subprocess
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

def create_kmer_database(fastafile:str, outputdir:str, prefix:str, kmersize=40):
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = os.getcwd()
    kmerdbroot = outputdir + "/" + prefix + ".kmers.k" + str(kmersize)
    if not os.path.exists(kmerdbroot + ".ktab"):
       command = "FastK -k" + str(kmersize) + " -T2 -N" + kmerdbroot + " -p -t1 -v " + fastafile
       print("Running: " + command)
       logger.info("Running: " + command)
       proc = subprocess.Popen(command, shell=True, env=env)
       proc.wait()

def remove_kmer_database(outputdir:str, prefix:str, kmersize=40):
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = os.getcwd()
    kmerdbroot = outputdir + "/" + prefix + ".kmers.k" + str(kmersize)
    if os.path.exists(kmerdbroot + ".ktab"):
       command = "Fastrm " + kmerdbroot
       print("Running: " + command)
       logger.info("Running: " + command)
       proc = subprocess.Popen(command, shell=True, env=env)
       proc.wait()

def find_anotb_kmers(db1prefix:str, db2prefix:str, outputdir:str, prefix:str, kmersize=40):
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = os.getcwd()
    kmerdbroot = outputdir + "/" + prefix + ".kmers.k" + str(kmersize)
    if not os.path.exists(kmerdbroot + ".ktab"):
       command = "Logex -T2 -h \"" + kmerdbroot + " = A-B\" " + outputdir + "/" + db1prefix + ".kmers.k" + str(kmersize) + " " + outputdir + "/" + db2prefix + ".kmers.k" + str(kmersize)
       print("Running: " + command)
       logger.info("Running: " + command)
       proc = subprocess.Popen(command, shell=True, env=env)
       proc.wait()

    return kmerdbroot

def map_kmer_markers_onto_fasta(fastafile:str, markerfilelist:list, outputdir:str):
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = os.getcwd()

    outputbedlist = []
    tmpdir = outputdir + "/tmp"
    path = Path(tmpdir)
    logger.info("Creating temporary directory " + tmpdir + " for output")
    path.mkdir(exist_ok=True)

    mergefilestring = ""
    for markerfile in markerfilelist:
        pattern = r"\.fastq$|\.fasta$|\.fq$|\.fa$|\.fastq.gz$|\.fasta.gz$|\.fq.gz$|\.fa.gz$"
        assemblyroot = re.sub(pattern, "", fastafile)

        # KmerMap constructs the output file name from the third argument with ".assemblyroot" + ".kmers.merge.bed"
        kmermapinputprefix = markerfile
        kmermapoutputprefix = markerfile + "." + assemblyroot
        mergefilestring = mergefilestring + " " + kmermapoutputprefix + ".kmers.merge.bed"
        markercommand = "KmerMap -v -m -T2 -P" + tmpdir + " " + markerfile + " " + fastafile + " " + kmermapinputprefix
        
        if not os.path.exists(kmermapoutputprefix + ".kmers.merge.bed"):
            logger.info("Running: " + markercommand)
            proc = subprocess.Popen(markercommand, shell=True, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            proc.wait()

            if proc.returncode == 0:
                logger.info(markercommand + " completed successfully")
            else:
                logger.info(stderr.decode())
        else:
            logger.info("Skipping run of " + markercommand + ": output already exists")
    
    outputfile = outputdir + "/" + assemblyroot + ".kmers.merge.bed"

    if not os.path.exists(outputfile):
        logger.debug("Merging output files into " + outputfile)
        mergesortcommand = "cat" + mergefilestring + " | sort -k1,1n -k2,2n -k3,3n " + " > " + outputfile
        logger.debug(mergesortcommand)
        os.system(mergesortcommand)
    else:
        logger.debug("Skipping merging of kmer files into " + outputfile + ": output already exists")
   
    return outputfile

def find_extreme_kmers(fastkdbroot:str):
    if os.path.exists(fastkdbroot + ".hist"):
        env = os.environ.copy()
        env['LD_LIBRARY_PATH'] = os.getcwd()
        command = "Tabex -t2 " + fastkdbroot + " LIST"
        print("Running: " + command)
        logger.info("Running: " + command)
        proc = subprocess.Popen(command, shell=True, env=env)
        while True:
            line = proc.stdout.readline()
            if not line:
                break
            print(line.strip())
        proc.wait()
    else:
        logger.critical("No such fastk database: " + fastkdbroot)

# returns the total count of base1 and base2 on the two strands of a kmer string
def twobase_kmer_comp(base1:str, base2:str, kmerseq:str):
    kmerseqlc = kmerseq.tolower()
    base1 = base1.tolower()
    base2 = base2.tolower()
    compbases = {'a':'t', 't':'a', 'g':'c', 'c':'g'}

    fwdcount = kmerseq.count(base1) + kmerseq.count(base2)
    revcount = kmerseq.count(compbase[base1]) + kmerseq.count(compbase[base2])

    return [fwdcount, revcount]

