#!/usr/bin/env python
# TODO: add an errorcheck make sure the barcodes are unique
# TODO: add verbosity option
# TODO: catch error if samtools quickcheck fails
# Usage: Merge BAM files from mulitple 10x runs, making unique CB tags across runs

import pysam
import gzip
import argparse
import multiprocessing 
import os
from posixpath import dirname
from pathlib import Path

def reader(f):
    """ read lines from a file """
    for line in f:
        yield line.rstrip()

def pool_barcodes(bc_files):
    """
    Update cell barcodes, making them unique across 10x runs
    """

    # initialize list of dictionaries for storing old and new barcodes 
    barcodes = [{} for _ in range(n)]

    # write all new barcodes to file
    pool_file = open(out_dir + "/" + base_name + "_barcodes.tsv", "w")
    
    # write map of new and old barcodes to file
    info_file = open(out_dir + "/" + base_name + "_info.tsv", "w")
    info_file.writelines("CB_pool\tCB_origin\tRun_id\n")

    # iterate through each sample
    for s in range(n):
        file = gzip.open(bc_files[s]) if Path(bc_files[s]).suffix == '.gz' else open(bc_files[s])
        with file as infile:
            for old_bc in reader(infile):
                new_bc = old_bc[:-1]+str(s+1) # make new barcode
                barcodes[s][old_bc] = new_bc # save to dictionary
                
                pool_file.writelines(new_bc + "\n") # save to file

                info_out = [new_bc, old_bc, str(s + 1)]
                info_file.writelines("\t".join(info_out) + "\n") # save new/old map to file
                
    pool_file.close()
    info_file.close()

    return barcodes

def retag_bam(barcodes, inBAM, outBAM):
    """
        read the reads from inBAM 
        update each read's cell barcode in the CB tag
        barcodes is a dictionary of old (keys) and new (values) barcodes
        write updated reads to outBAM 
    """

    reads = pysam.AlignmentFile(inBAM, "rb")
    head = reads.header.to_dict()
    out = pysam.AlignmentFile(outBAM, "wb", header=head)

    # iterate through each read
    for read in reads:
        # check to see whether the CB tag needs to be changed
        if read.has_tag('CB') and read.get_tag('CB') in barcodes:
            # set the new CB tag
            read.set_tag('CB', barcodes[read.get_tag('CB')])
        out.write(read)

def main(bc_files, bam_files, outBAM):
    """
        run the script
    """
    # 1. pool barcodes together
    print("Creating new cell barcodes...")
    barcodes = pool_barcodes(bc_files)

    # 2. retag each bam file
    # setup pool of workers
    p = multiprocessing.Pool(processes=n)
    temp_name = "retagged_temp"
    for s in range(n):
        print("Retagging {}...".format(bam_files[s]))
        tempBAM = out_dir + "/" + temp_name + "%d.bam" %(s)
        p.apply_async(retag_bam, (barcodes[s], bam_files[s], tempBAM))
    p.close()
    p.join()

    # 3. merge retagged bam files
    print("Merging retagged bam files...")
    file_list = [out_dir + "/" + temp_name + "%d.bam" %(s) for s in range(len(bam_files))]
    merge_params = [outBAM] + ["-f"] + file_list
    pysam.merge(*merge_params)
    [os.remove(file) for file in file_list]

    # 4. check the output bam file for integrity
    print("Running samtools quickcheck {} ...".format(outBAM))
    pysam.quickcheck(outBAM)

    # 5. sort and index the output bam file
    print("Running samtools sort {}...".format(outBAM))
    sortedBAM = out_dir + "/" + base_name + "_sorted.bam"
    sort_params = ["-o", sortedBAM, outBAM]
    pysam.sort(*sort_params)
    os.remove(outBAM)
    
    print("Running samtools index {}...".format(outBAM))
    pysam.index(sortedBAM)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge BAM files from mulitple 10x runs, making unique CB tags across runs")
    parser.add_argument(
        "-o", "--out", type=Path, help="the filename to save the merged BAM file"
    )
    parser.add_argument(
        "--barcodes-files", type=str, nargs="+", help="input barcodes files"
    )
    parser.add_argument(
        "--bam-files", type=str, nargs="+", help="input BAM files"
    )
    args = parser.parse_args()

    # error checks
    if len(args.barcodes_files) != len(args.bam_files):
        raise ValueError("Number of barcodes files and BAM files must be the same.")

    n = len(args.bam_files)
    
    for s in range(n):
        if not os.path.exists(args.barcodes_files[s]):
            raise FileNotFoundError("Barcodes file {} does not exist".format(args.barcodes_files[s]))
        if not os.path.exists(args.bam_files[s]):
            raise FileNotFoundError("BAM file {} does not exist".format(args.bam_files[s]))

    if os.path.exists(args.out):
        raise ValueError("Output file {} already exists.".format(args.out))

    # make global variables
    out_dir = str(args.out.parent)
    base_name = args.out.stem

    # read barcodes
    print("Merging {} runs: ".format(n) + " ".join(args.bam_files))
    main(args.barcodes_files, args.bam_files, str(args.out))

