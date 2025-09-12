import sys
import argparse
import pysam

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="From a tab-delimited file of coordinates, create a VCF file to patch a reference fasta"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 1.0.0"
    )
    parser.add_argument('-c', '--coords', required=True, help='bed-like file of 1-based interval locations (ref, then patch) to include in vcf')
    parser.add_argument('-r', '--reffasta', type=str, required=True, help='fasta file for assembly reference')
    parser.add_argument('-p', '--patchfasta', type=str, required=True, help='fasta file for patches')
    parser.add_argument('--vcf', type=str, required=True, help='name of vcf file to output')
    parser.add_argument('--header', type=str, required=False, help='name of vcf file header to include (optional)')
    return parser

def write_vcf_file(args) -> None:
    reffastafile = args.reffasta
    refobj = pysam.FastaFile(reffastafile)
    patchfastafile = args.patchfasta
    patchobj = pysam.FastaFile(patchfastafile)
    vcffile = args.vcf
    headerfile = args.header

    coordfile = args.coords
    varcount = 1
    ofh=open(vcffile, 'w')
    if (headerfile is not None):
        with open(headerfile) as hfh:
            for line in hfh.readlines():
                ofh.write(line)

    with open(coordfile) as fh:
        for line in fh.readlines():
            fields = line.rstrip().split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            patchchrom = fields[3]
            patchstart = int(fields[4])
            patchend = int(fields[5])
            strand = fields[6]
            # assume passed coordinate file is one-based (pysam considers start 0-based)

            if start > end or patchstart > patchend:
                print("Ignoring patch line with starting point larger than ending point")
                continue

            refseq = refobj.fetch(chrom, start - 1, end).upper()
            
            patchseq = patchobj.fetch(patchchrom, patchstart - 1, patchend).upper()

            if strand == "-" or strand == "R":
                complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
                patchseq = patchseq[::-1]
                patchseq_rc = "".join(complement.get(base, base) for base in patchseq)
                patchseq = patchseq_rc
                [refseq, patchseq] = trimend(refseq, patchseq) 
                [refseq, patchseq, start] = trimstart(refseq, patchseq, start)
            else:
                [refseq, patchseq] = trimend(refseq, patchseq) 
                [refseq, patchseq, start] = trimstart(refseq, patchseq, start)

            ofh.write(chrom + "\t" + str(start) + "\t" + str(varcount) + "\t" + refseq + "\t" + patchseq + "\t0\tPASS\t.\tGT\t1/1\n")
            varcount = varcount+1

    return

def trimstart(ref, patch, pos):
    # trim starts of alleles where they are the same
    charsremoved = 0
    for r, p in zip(ref, patch):
        if r == p and len(ref)>1 and len(patch)>1:
            ref = ref[1:]
            patch = patch[1:]
            pos = pos + 1
            charsremoved = charsremoved + 1
        else:
            break
    print("Removed " + str(charsremoved) + " bases from start of ref and patch")

    changedvals = [ref, patch, pos]
    return changedvals

def trimend(ref, patch):
    # trim identical ends of alleles as well
    charsremoved = 0
    revref = ref[::-1]
    revpatch = patch[::-1]
    for r, p in zip(revref, revpatch):
        if r == p and len(ref)>1 and len(patch)>1:
            ref = ref[:-1]
            patch = patch[:-1]
            charsremoved = charsremoved + 1
        else:
            break
    print("Removed " + str(charsremoved) + " bases from end of ref and patch")

    changedvals = [ref, patch]
    return changedvals

def main() -> None:
    parser = init_argparse()
    args = parser.parse_args()

    write_vcf_file(args)

if __name__ == "__main__":
    main()
