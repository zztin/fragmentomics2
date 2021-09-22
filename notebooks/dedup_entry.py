import pysam
import singlecellmultiomics
from singlecellmultiomics.molecule import MoleculeIterator
import os
from singlecellmultiomics.utils.sequtils import hamming_distance
from singlecellmultiomics.bamProcessing import sorted_bam_file
from singlecellmultiomics.fragment import FragmentWithoutUMI
from singlecellmultiomics.fragment import Fragment

import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Deduplicate base on molecular location.')
    parser.add_argument('--read_bam', type=str,
                    help='bam file path')
    parser.add_argument('--out_path', type=str,
                    help='out path')
    parser.add_argument('--prefix', type=str, help='prefix of the bam file name')
    parser.add_argument('--SM', type=str, help="value to put in SM tag")
    args = parser.parse_args()
    SM_bam = f"{args.out_path}/{args.prefix}_test.SMtagged.sorted.bam"
    t_bam = f"{args.out_path}/{args.prefix}_deduplicated_test.tagged.bam"
    with pysam.AlignmentFile(args.read_bam, "r", ignore_truncation=True) as g:
        with sorted_bam_file(SM_bam, origin_bam=g) as f:
            for i, read in enumerate(g):
                read.set_tag("SM", args.SM)
                f.write(read)

    with pysam.AlignmentFile(SM_bam) as f:
        with sorted_bam_file(t_bam, origin_bam=f, ) as target_bam:
            for i,m in enumerate(MoleculeIterator(
                            alignments=f,
                            moleculeClass=singlecellmultiomics.molecule.Molecule,
                            fragmentClass=Fragment,
                            every_fragment_as_molecule=False,
                            perform_qflag=False
                            )):
                read_name = f'consensus_{m.get_a_reference_id()}_{i}',
                m.write_tags()
                m.write_pysam(target_bam)
    pysam.index(t_bam, f'{t_bam}.bai')
    #os.remove(SM_bam)
