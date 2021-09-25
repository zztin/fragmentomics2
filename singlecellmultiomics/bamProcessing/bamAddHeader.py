from singlecellmultiomics.bamProcessing import sorted_bam_file
import pysam


def bam_add_header(bam_file, SM_tag, out_path):
    with pysam.AlignmentFile(bam_file, "r", ignore_truncation=True) as g:
        with sorted_bam_file(out_path, origin_bam=g) as f:
            for i, read in enumerate(g):
                read.set_tag("SM", SM_tag)
                f.write(read)

