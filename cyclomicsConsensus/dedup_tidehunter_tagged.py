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
from singlecellmultiomics.bamProcessing.bamFunctions import get_reference_from_pysam_alignmentFile
import argparse
import time
import linecache
import os
import tracemalloc

def display_top(snapshot, filename, i, key_type='lineno', limit=10):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)
    records = []
    records.append("Top %s lines by molecule no.%s" % (limit, i))
    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        print("#%s: %s:%s: %.1f KiB" % (index, frame.filename, frame.lineno, stat.size / 1024))
        records.append("#%s: %s:%s: %.1f KiB"
              % (index, frame.filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            records.append('    %s' % line)
            print('   %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        records.append("%s other: %.1f KiB" % (len(other), size / 1024))
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    records.append("Total allocated size: %.1f KiB" % (total / 1024))
    print("Total allocated size: %.1f KiB" % (total / 1024))
    with open(f"{args.out_path}/../{filename}.tracemalloc.log", "w+") as f:
        f.write("\n".join(records))



if __name__=="__main__":
    pid = os.getpid()
    parser = argparse.ArgumentParser(description='Deduplicate base on molecular location.')
    parser.add_argument('--read_bam', type=str,
                    help='bam file path')
    parser.add_argument('--out_path', type=str,
                    help='out path')
    parser.add_argument('--prefix', type=str, help='prefix of the bam file name')
    parser.add_argument('--SM', type=str, help="value to put in SM tag")
    parser.add_argument('--ref', type=str, help='reference which bam file is mapped to. (str)')
    parser.add_argument('--merge', type=str, help='merge reads from different nanopore reads but covering the same start,'
                                                  'end sites within max 2bp range. Direction of the read is set by the first added read.')

    args = parser.parse_args()
    SM_bam = f"{args.out_path}/{args.prefix}_{pid}_test.SMtagged.sorted.bam"
    t_bam = f"{args.out_path}/{args.prefix}_{pid}_deduplicated_samecoor.tagged.bam"

    # autodetect reference:
    reference = None
    if args.ref is None:
        args.ref = get_reference_from_pysam_alignmentFile(args.read_bam)

    if args.ref is not None:
        try:
            reference = UnitialisedClass(CachedFastaNoHandle, args.ref)
            print(f'Loaded reference from {args.ref}')
        except Exception as e:
            print("Error when loading the reference file, continuing without a reference")
            reference = None

    if args.SM is not None:
        with pysam.AlignmentFile(args.read_bam, "r", ignore_truncation=True) as g:
            with sorted_bam_file(SM_bam, origin_bam=g) as f:
                for i, read in enumerate(g):
                    read.set_tag("SM", args.SM)
                    f.write(read)
    else:
        SM_bam = args.read_bam
    timeA = time.time()
    print('SM tag written.', timeA)



    tracemalloc.start()
    with pysam.AlignmentFile(SM_bam) as f:
        
       
        with sorted_bam_file(t_bam, origin_bam=f, ) as target_bam:
            
#            snapshot2 = tracemalloc.take_snapshot()
            for i, m in enumerate(MoleculeIterator(
                            alignments=f,
                            moleculeClass=singlecellmultiomics.molecule.chic.CHICMolecule,
                            fragmentClass=CHICFragment,
                            every_fragment_as_molecule=False,
                            perform_qflag=False,
                            molecule_class_args={"reference": reference, "max_associated_fragments": 100},
                            fragment_class_args={"assignment_radius": 4},
                            max_buffer_size=1000000,
                            yield_overflow=False,
    
            )):
                
#                snapshot3 = tracemalloc.take_snapshot()
                read_name = f'consensus_{m.get_a_reference_id()}_{i}'
                # write tags to all fragments associated with the molecule

                snapshot1 = tracemalloc.take_snapshot()

                m.write_tags()

                m.write_pysam(target_bam,
                              consensus=True,
                              consensus_name=read_name,
                              no_source_reads=True,
                              )

#                top_stats = snapshot3.compare_to(snapshot2, 'lineno')
#                print("[ Top 10 differences ]")
#                for stat in top_stats[:10]:
#                    print(stat)

                if i > 3855:
                        print('snapshot',i, )
                        print(m.span)
                        print(len(m.read_names))
                        print(m.overflow_fragments)
                        print(len(m.fragments))
                        snapshot2 = tracemalloc.take_snapshot()
                        display_top(snapshot2, filename=f'{args.prefix}_{pid}', i = i)
                        top_stats = snapshot2.compare_to(snapshot1, 'lineno')
                        print("[ Top 10 differences ]")
                        for stat in top_stats[:10]:
                            print(stat)


    pysam.index(t_bam, f'{t_bam}.bai')
    if args.SM is not None:
        os.remove(SM_bam)
        os.remove(f'{SM_bam}.bai')

print((time.time() - timeA)/60, 'min')
