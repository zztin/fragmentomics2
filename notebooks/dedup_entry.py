import pysam
import singlecellmultiomics
from singlecellmultiomics.molecule import MoleculeIterator
import os
from singlecellmultiomics.utils.sequtils import hamming_distance
from singlecellmultiomics.bamProcessing import sorted_bam_file
from singlecellmultiomics.fragment import FragmentWithoutUMI
from singlecellmultiomics.fragment import CHICFragment
from singlecellmultiomics.utils.prefetch import UnitialisedClass
from singlecellmultiomics.fastaProcessing import CachedFastaNoHandle

import argparse




if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Deduplicate base on molecular location.')
    parser.add_argument('--read_bam', type=str,
                    help='bam file path')
    parser.add_argument('--out_path', type=str,
                    help='out path')
    parser.add_argument('--prefix', type=str, help='prefix of the bam file name')
    parser.add_argument('--SM', type=str, help="value to put in SM tag")
    parser.add_argument('--ref', type=str, help='reference which bam file is mapped to. (str)')
    args = parser.parse_args()
    SM_bam = f"{args.out_path}/{args.prefix}_test.SMtagged.sorted.bam"
    t_bam = f"{args.out_path}/{args.prefix}_deduplicated_test.tagged.bam"

    # autodetect reference:
    reference = None
    if args.ref is None:
        args.ref = get_reference_from_pysam_alignmentFile(input_bam)

    if args.ref is not None:
        try:
            reference = UnitialisedClass(CachedFastaNoHandle, args.ref)
            print(f'Loaded reference from {args.ref}')
        except Exception as e:
            print("Error when loading the reference file, continuing without a reference")
            reference = None


    with pysam.AlignmentFile(args.read_bam, "r", ignore_truncation=True) as g:
        with sorted_bam_file(SM_bam, origin_bam=g) as f:
            for i, read in enumerate(g):
                read.set_tag("SM", args.SM)
                f.write(read)

    with pysam.AlignmentFile(args.read_bam) as f:
        with sorted_bam_file(t_bam, origin_bam=f, ) as target_bam:
            for i, m in enumerate(MoleculeIterator(
                            alignments=f,
                            moleculeClass=singlecellmultiomics.molecule.chic.CHICMolecule,
                            fragmentClass=CHICFragment,
                            every_fragment_as_molecule=False,
                            perform_qflag=False,
                            molecule_class_args={"reference": reference},
#                            fragment_class_args=,

            )):
                read_name = f'consensus_{m.get_a_reference_id()}_{i}'
                m.write_tags()
                m.write_pysam(target_bam,
                              consensus=True,
                              consensus_name=read_name)
    pysam.index(t_bam, f'{t_bam}.bai')
    os.remove(SM_bam)
    os.remove(f'{SM_bam}.bai')
