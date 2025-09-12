import sys
import argparse
import pysam
import re
import random
from collections import namedtuple

# create namedtuple for bed intervals:
bedinterval = namedtuple('bedinterval', ['chrom', 'start', 'end', 'name', 'rest'])

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Search BAM files for sequence spanning specified reference regions and report it, and optionally report variants in the alignment"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 1.0.0"
    )
    parser.add_argument('-b', '--bams', required=True, help='comma-delimited list of assembly bam files to search for spanning sequences (required)')
    parser.add_argument('-f', '--fastas', type=str, required=True, help='comma-delimited list of fasta files for the assemblies whose bam files are given by --bams (required)')
    parser.add_argument('-r', '--reffasta', type=str, required=True, help='fasta file for the common reference for the bam files (optional)')
    parser.add_argument('--region', type=str, required=False, help='single region of reference to search (optional)')
    parser.add_argument('--bed', type=str, required=False, help='bed file of reference regions to search (optional)')
    parser.add_argument('--prefix', type=str, required=False, default="assemblypatch", help='prefix to use in names of output files (optional, default assemblypatch)')
    return parser

def create_interval_from_string(interval:str):
    regmat = re.compile(r'(\S+):(\d+)\-(\d+)')
    regionmatch = regmat.match(interval)
    if regionmatch:
        bedint = bedinterval(chrom=regionmatch.group(1), start=int(regionmatch.group(2)), end=int(regionmatch.group(3)), name = "", rest="rest")
    else:
        bedint = None
    return bedint

def create_intervallist_from_bedfile(bedfile:str):
    intlist = []
    intnumber = 1
    with open(bedfile, "r") as bfh:
        intervalline = bfh.readline()
        while intervalline:
            intervalline = intervalline.rstrip()
            fields = intervalline.split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if len(fields) >= 4:
                name = fields[3]
            else:
                name = "Region" + str(intnumber)
                intnumber = intnumber + 1
            intlist.append(bedinterval(chrom=chrom, start=start, end=end, name=name, rest=""))
            intervalline = bfh.readline()

    return intlist

def compare_sequence_from_bams(bamobjects, fastaobjects, reffastaobj, regions, args):

    patchfilename = args.prefix + ".patchcoords.txt"
    with open(patchfilename, "w") as pfh:
        for region in regions:
            chrom = region.chrom
            start = region.start
            # because we'll be looking at flanks, we need to inch in from ends of chrom--should do this on the right as well
            if start == 0:
                start = start + 1
            end = region.end
            refseq = reffastaobj.fetch(reference=chrom, start=start, end=end)
            print("REGION " + chrom + ":" + str(start) + "-" + str(end))
            correctionseqs = []
            # correction contig and positions will contain just the first bamfile's coordinates/strand
            correctioncontig = None
            correctionstrand = None
            correctionstart = None
            correctionend = None
            for bamindex in range(len(bamobjects)):
                bamobj = bamobjects[bamindex]
                fastaobj = fastaobjects[bamindex]
                # find_flank_positions routine takes one-based ref coordinates as arguments:
                [contig, leftcoord, rightcoord, strand, spanningseq] = find_flank_positions(bamobj, fastaobj, chrom, start, end+1)
                if leftcoord is not None and rightcoord is not None:
                    # trim extra base on each side:
                    spanningseq = spanningseq[1:-1]
                    correctionseqs.append(spanningseq)
                    if bamindex == 0:
                        # still 1-based, just switching from left/right to actual query coords for reverse strand
                        if strand == 'F':
                            correctionstart = leftcoord + 1
                            correctionend = rightcoord - 1
                        else:
                            correctionstart = rightcoord + 1
                            correctionend = leftcoord - 1
                        correctioncontig = contig
                        correctionstrand = strand

            if len(correctionseqs) == len(bamobjects):
                match = True
                potentialcorrection = correctionseqs[0]
                #print("0:" + potentialcorrection)
                for index in range(len(correctionseqs)):
                    if index == 0:
                        continue
                    #print(str(index) + ":" + correctionseqs[index])
                    if correctionseqs[index] != potentialcorrection:
                        print(potentialcorrection + "\n" + correctionseqs[index])
                        print("ASSEMBLYMISMATCH")
                        match = False
                        break
                if match and refseq != correctionseqs[0]:
                    print("ALLMATCH " + correctioncontig + ":" + str(correctionstart) + "-" + str(correctionend) + "/" + correctionstrand)
                    print("REFSEQ " + refseq)
                    print("CORREC " + correctionseqs[0])
                    # check for zero length alleles:
                    if correctionstart > correctionend or start > end - 1:
                        start = start - 1
                        if correctionstrand == 'F':
                            correctionstart = correctionstart - 1
                        else:
                            correctionend = correctionend + 1
                    if correctionstrand == 'F':
                        plusminusstrand = '+'
                    else:
                        plusminusstrand = '-'
                    pfh.write(chrom + "\t" + str(start + 1) + "\t" + str(end) + "\t" + correctioncontig + "\t" + str(correctionstart) + "\t" + str(correctionend) + "\t" + correctionstrand + "\n")
                elif match:
                    print("AGREEWITHBENCH")
            else:
                print("UNSPANNED")

    return 0

# goal of this routine is to return 1-based coordinates along the query corresponding to the 1-based ref coords passed to the routine
def find_flank_positions(bamobj, fastaobj, chrom, start, end):

    leftcoord = None
    rightcoord = None
    for assemblyalign in bamobj.fetch(contig=chrom, start=start-1, stop=end):
        if assemblyalign.is_secondary:
            continue
        # retrieve one-based coordinates of the alignment (querystart < queryend regardless of strand)
        [contigname, querystart, queryend, ref, refstart, refend, strand] = retrieve_align_data(assemblyalign)
        print(contigname + ":" + str(querystart) + "-" + str(queryend) + "/" + strand)
        print(ref + ":" + str(refstart) + "-" + str(refend))

        cigartuples = assemblyalign.cigartuples
        if cigartuples[-1][0] == 5:
            righthardclip = cigartuples[-1][1]
        else:
            righthardclip = 0
        if cigartuples[0][0] == 5:
            lefthardclip = cigartuples[0][1]
        else:
            lefthardclip = 0
        if cigartuples[-1][0] == 4:
            rightsoftclip = cigartuples[-1][1]
        else:
            rightsoftclip = 0
        if cigartuples[0][0] == 4:
            leftsoftclip = cigartuples[0][1]
        else:
            leftsoftclip = 0

        #print("Hard clipping: " + str(lefthardclip) + "/" + str(righthardclip) + " Soft clipping: " + str(leftsoftclip) + "/" + str(rightsoftclip))
        pairs = assemblyalign.get_aligned_pairs()
        if len(pairs) < 1000000:
            continue

        # for F strand, left coord returned by find_querypos_in_pairs is 0-based from left end of read, not counting hard clipping (at least in this version of pysam!)
        if refstart <= start and refend >= start:
            leftcontig = contigname
            #leftcoordzb = find_querypos_in_pairs(pairs, start-1)
            leftcoordzb = find_readpos_in_pairs(pairs, start-1)
            if leftcoordzb is not None:
                if strand == 'F':
                    leftcoord = leftcoordzb + lefthardclip + 1
                if strand == 'R':
                    leftcoord = queryend - leftcoordzb + leftsoftclip
        if refstart <= end and refend >= end:
            rightcontig = contigname
            #rightcoordzb = find_querypos_in_pairs(pairs, end-1)
            rightcoordzb = find_readpos_in_pairs(pairs, end-1)
            if rightcoordzb is not None:
                if strand == 'F':
                    rightcoord = rightcoordzb + lefthardclip + 1
                if strand == 'R':
                    rightcoord = queryend - rightcoordzb + leftsoftclip

    if leftcoord is not None and rightcoord is not None and leftcontig==rightcontig:
        print(chrom + ":" + str(start) + "-" + str(end) + " aligns to " + contigname + ":" + str(leftcoord) + "-" + str(rightcoord))
        if strand == 'F' and leftcoord < rightcoord+1:
            spanningseq = fastaobj.fetch(reference=leftcontig, start=leftcoord-1, end=rightcoord)
        elif rightcoord < leftcoord+1:
            spanningseq = fastaobj.fetch(reference=leftcontig, start=rightcoord-1, end=leftcoord)
            spanningseq = revcomp(spanningseq)
        else:
            flankpositions = [None, None, None, None, None]
            return flankpositions

        flankpositions = [leftcontig, leftcoord, rightcoord, strand, spanningseq]
        return flankpositions
    else:
        flankpositions = [None, None, None, None, None]
        return flankpositions


def retrieve_align_data(align)->list:
    if align.is_reverse:
        strand = 'R'
    else:
        strand = 'F'
    query = align.query_name

    # number of hard clipped bases:
    cigartuples = align.cigartuples
    # will use actual one-based positions, so I don't go crazy, then report BED format zero-based half open later
    # pysam's align.query_alignment_start is the 0-based coordinate within the uncomplemented hard clipped query sequence, so here we add hardclipping from the
    # appropriate end to get coordinates within the entire sequence (so low coordinates are towards the start of the original query sequence, not the left end
    # of the alignment)
    if strand == 'F':
        if cigartuples[0][0] == 5:
            hardclip = cigartuples[0][1]
        else:
            hardclip = 0
        querystart = hardclip + align.query_alignment_start + 1 # first unclipped base (was 0-based, now 1-based)
        queryend = hardclip + align.query_alignment_end # was 0-based index just past the last unclipped base now 1-based last unclipped base
    elif strand == 'R':
        if cigartuples[-1][0] == 5:
            hardclip = cigartuples[-1][1]
        else:
            hardclip = 0
        querystart = hardclip + align.query_length - align.query_alignment_end + 1 # was 0-based index just past end of reversed query, now 1-based 5' end of query
        queryend = hardclip + align.query_length - align.query_alignment_start # was 0-based index of first unclipped base, now 1-based 3' unclipped end of query

    ref = align.reference_name

    refstart = align.reference_start + 1 # was zero-based, now one-based
    refend = align.reference_end #  was zero-based, adding one for bed format "half open", now one-based

    aligndata = [query, querystart, queryend, ref, refstart, refend, strand]

    return aligndata

# "pairs" is a structure created by pysam for an alignment, and contains a list
#   of tuples "consisting of the 0-based offset from the start of the read 
#   sequence followed by the 0-based reference position
#
# This routine searches the list of tuples for the one corresponding to the 
#   desired reference position ("pos"), and returns the read position
#   for that tuple, assuming it exists
#
# Why would a tuple have a ref position of "None"? These tuples are insertions
#   in the read sequence (and deletions from the reference have read pos of 
#   None). In the case of deletions from the reference, the "ifnone" argument
#   can tell the routine whether to report the read position that is one to 
#   the left of the deleted ref location ("lower") or one to the right ("higher")
#
# Note: this routine takes a *zero-based* position as its "pos" argument!
#
def find_readpos_in_pairs(pairs, pos, ifnone="lower"):
    # number of tuples:
    alignlength = len(pairs)
    ilow = 0
    ihigh = alignlength - 1

    hiref = pairs[ihigh][1]
    while hiref is None and ihigh > 0:
        ihigh = ihigh - 1
        hiref = pairs[ihigh][1]
    lowref = pairs[ilow][1]
    while lowref is None and ilow < len(pairs):
        ilow = ilow + 1
        lowref = pairs[ilow][1]

    # these shouldn't really happen
    if lowref is None or hiref is None or lowref > pos or hiref < pos or lowref > hiref:
        return None

    # need to be careful here to avoid an infinite loop:
    lastimid = None
    while ihigh - 1 > ilow and hiref >= pos and lowref <= pos:
        imid = int((ilow + ihigh)/2)
        midref = pairs[imid][1]
        #print("Pos: " + str(pos) + " (" + str(lowref) + "/" + str(hiref) + ") " + str(imid) + "/" + str(ilow) + "/" + str(ihigh))
        # if imid tuple has an inserted sequence
        while midref is None:
            if ifnone == "lower":
                imid = imid - 1
            else:
                imid = imid + 1
            midref = pairs[imid][1]
            if (ifnone == "lower" and imid <= ilow) or (ifnone=="higher" and imid >= ihigh):
                break
        if midref is None:
            return None
        if lastimid is not None and lastimid == imid:
            for i in range(ilow, ihigh):
                if pairs[i][1] is not None and pairs[i][1] == pos:
                    return pairs[i][0]
            return None
        
        if midref == pos:
            return pairs[imid][0] # might be None if nothing is aligned here
        elif midref > pos:
            ihigh = imid
            hiref = pairs[ihigh][1]
        elif midref < pos:
            ilow = imid
            lowref = pairs[ilow][1]
        if (pairs[ilow][1] is not None and pairs[ilow][1] > pos) or (pairs[ihigh][1] is not None and pairs[ihigh][1] < pos):
            return None
        lastimid = imid

    return None

def find_querypos_in_pairs(pairs, pos, ifnone="lower"):
    alignlength = len(pairs)
    low = 0
    high = alignlength - 1

    hiref = pairs[high][1]
    while hiref is None and high > 0:
        high = high - 1
        hiref = pairs[high][1]
    lowref = pairs[low][1]
    while lowref is None and low < len(pairs):
        low = low + 1
        lowref = pairs[low][1]

    if lowref is None or hiref is None or lowref > pos or hiref < pos:
        return None

    while high > low and hiref >= pos and lowref <= pos:
        mid = int((low + high)/2)
        if random.random() >= 0.5:
            while pairs[mid][1] is None and mid <= high:
                mid = mid + 1
        else: # take first non-None refpos to the left
            while pairs[mid][1] is None and mid >= low:
                mid = mid - 1
        midval = pairs[mid][1]
        if midval is None:
            return None
        #print(str(low) + " " + str(high) + " " + str(midval) + " " + str(pos) + " " + str(pairs[low][1]) + " " + str(pairs[high][1]) + " " + str(len(pairs)))
        if midval == pos:
            return pairs[mid][0] # might be None if nothing is aligned here
        elif midval > pos:
            high = mid
        elif midval < pos:
            low = mid
        if (pairs[low][1] is not None and pairs[low][1] > pos) or (pairs[high][1] is not None and pairs[high][1] < pos):
            return None

    return None

def revcomp(seq:str) -> str:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'M': 'N', 'K': 'N', 'R': 'N', 'W': 'N', 'Y': 'N'}
    bases = list(seq)
    bases = bases[::-1]
    bases = [complement[base] for base in bases]

    return ''.join(bases)

def main() -> None:
    parser = init_argparse()
    args = parser.parse_args()

    regions = []
    if args.region:
        regions.append(create_interval_from_string(args.region))
    elif args.bed:
        regions = create_intervallist_from_bedfile(args.bed)

    bamfiles = args.bams.split(",")
    bamobjects = []
    for bam in bamfiles:
        samobject = pysam.AlignmentFile(bam, "rb")
        bamobjects.append(samobject)

    fastafiles = args.fastas.split(",")
    fastaobjects = []
    for fasta in fastafiles:
        fastaobj = pysam.FastaFile(fasta)
        fastaobjects.append(fastaobj)
    reffastaobj = pysam.FastaFile(args.reffasta)

    if len(bamobjects) != len(fastaobjects):
        print("The number of bamfiles given by --bams must be the same as the number of fasta files given by --fastas")
        sys.exit()
    else:
        print(str(len(bamobjects)) + " bam files passed as --bams")
        print(str(len(regions)) + " regions")

    compare_sequence_from_bams(bamobjects, fastaobjects, reffastaobj, regions, args)

if __name__ == "__main__":
    main()
