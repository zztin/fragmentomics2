import os
import pysam
import time
import contextlib
from shutil import which

def get_reference_from_pysam_alignmentFile(pysam_AlignmentFile, ignore_missing=False):
    """Extract path to reference from pysam handle

    Args:
        pysam_AlignmentFile (pysam.AlignmentFile)
        ignore_missing(bool) : Check if the file exists, if not return None
    Returns:
        path : path to bam file (if exists or ignore_missing is supplied) or None
    """
    try:
        for x in pysam_AlignmentFile.header.as_dict()['PG']:
            if  x.get('ID')!='bwa':
                continue
            for argument in x.get('CL').split():
                if (argument.endswith('.fa') or argument.endswith('.fasta') or argument.endswith('.fasta.gz') or argument.endswith('.fa.gz')) and (ignore_missing or os.path.exists(argument)):
                    return argument
    except Exception as e:
        pass


def get_reference_path_from_bam(bam, ignore_missing=False ):
        """Extract path to reference from bam file

        Args:
            bam(str) or pysam.AlignmentFile : path to bam file
            ignore_missing(bool) : Check if the file exists, if not return None
        Returns:
            path : path to bam file (if exists or ignore_missing is supplied) or None
        """
        if type(bam) is str:
            with pysam.AlignmentFile(bam) as pysam_AlignmentFile_handle:
                return  get_reference_from_pysam_alignmentFile( pysam_AlignmentFile_handle, ignore_missing=ignore_missing)
        elif type(bam) is pysam.AlignmentFile:
            return  get_reference_from_pysam_alignmentFile(bam, ignore_missing=ignore_missing)

        else:
            raise ValueError('Supply either a path to a bam file or pysam.AlignmentFile object')


@contextlib.contextmanager
def sorted_bam_file( write_path,origin_bam=None, header=None,read_groups=None):
    """ Get writing handle of a sorted bam file
    Args:
        write_path (str) : write to a  bam file at this path
        origin_bam (pysam.AlignmentFile ) : bam file to copy header for or
        header (dict) : header for the bam file to write
        read_groups(set/dict) : set or dictionary which contains read groups. The dictionary should have the format { read_group_id (str)
                { 'ID': ID, 'LB':library,
                'PL':platform,
                'SM':sampleLib,
                'PU':readGroup }

    """
    try:
        unsorted_path = f'{write_path}.unsorted'
        if header is not None:
            pass
        elif origin_bam is not None:
            header = origin_bam.header.copy()
        else:
            raise ValueError("Supply a header or origin_bam object")
        unsorted_alignments = pysam.AlignmentFile(unsorted_path, "wb", header=header)
        yield unsorted_alignments
    finally:
        unsorted_alignments.close()
        if read_groups is not None:
            add_readgroups_to_header( unsorted_path, read_groups )

        # Write, sort and index
        sort_and_index(unsorted_path,  write_path, remove_unsorted=True)

def write_program_tag(input_header,
    program_name,
    command_line,
    description,
    version
    ):
    """Write Program Tag to bam file header
    Args:
        input_header  (dict): header to write PG tag to
        program_name (str) : value to write to PN tag
        command_line (str) : value to write to CL tag
        version (str) : value to write to VN tag
        description (str) : value to write to DS tag
    """
    if not 'PG' in input_header:
        input_header['PG'] = []

    input_header['PG'].append({
        'PN' : program_name,
        'CL' : command_line,
        'VN' : version,
        'DS' : description
    })


def sort_and_index(unsorted_path, sorted_path, remove_unsorted=False):
    """ Sort and index a bam file
    Args:
        unsorted_path (str) : path to unsorted bam file
        sorted_path (str) : write sorted file here
        remove_unsorted (bool) : remove the unsorted file

    Raises:
        SamtoolsError when sorting or indexing fails
    """
    pysam.sort("-o", sorted_path, unsorted_path)
    pysam.index(sorted_path)
    if remove_unsorted:
        os.remove(unsorted_path)

def GATK_indel_realign( origin_bam, target_bam,
    contig, region_start, region_end,
    known_variants_vcf_path,
    #realignerTargetCreatorArgs=None,
    #indelRealignerArgs=None,
    gatk_path = 'GenomeAnalysisTK.jar',
    interval_path=None,
    java_cmd = 'java -jar -Xmx8G -Djava.io.tmpdir=./gatk_tmp',
    reference = None,
    interval_write_path=None
     ):
    """
    Re-align a specified region in a bam file using GenomeAnalysisTK

    origin_bam (str) :  path to extract region from to re-align
    target_bam(str) : path to write re-aligned reads to
    contig (str) : contig of selected region
    region_start (int) : start coordinate of region to re align (1 based)
    region_end (int) :end coordiante of  selected region (1 based)
    known_variants_vcf_path (str) : path to vcf containing reference variants
    interval_path (str) : Use this intervals to perform realignment, when not specified intervals are generated using RealignerTargetCreator
    interval_write_path (str) : when interval_path is not supplied, write the interval file here
    java_cmd (str) : Command to open java
    gatk_path (str) : path to GenomeAnalysisTK.jar
    reference (str) : path to reference Fasta file
    """
    if not os.path.exists('./gatk_tmp'):
        try:
            os.path.makedirs('./gatk_tmp')
        except Exception as e:
            pass


    # Detect reference from bam file
    if reference  is None:
        reference = get_reference_path_from_bam(origin_bam)
    if reference is None:
        raise ValueError('Supply a path to a reference Fasta file (reference)')

    if interval_path is None:
        if interval_write_path is None:
            interval_write_path = f"{target_bam.replace('.bam','')}.intervals"

        target_creator_cmd = f'{java_cmd} {gatk_path} \
        -T RealignerTargetCreator \
        -R {reference} \
        -L {contig}:{region_start}-{region_end} \
        -known {known_variants_vcf_path} \
        -I {origin_bam} \
        -o {interval_path}'

        # Create the intervals file
        os.system(target_creator_cmd)
        interval_path = interval_write_path

    # Perform realignment
    realign_cmd = f'{java_cmd} {gatk_path} \
    -T IndelRealigner \
    -R {reference} \
    -targetIntervals {interval_path} \
    -known {known_variants_vcf_path} \
    -L {contig}:{region_start}-{region_end} \
    -I {origin_bam} \
    -o {target_bam}'
    os.system(realign_cmd)
    return target_bam



def add_readgroups_to_header( origin_bam_path, readgroups_in, target_bam_path=None, header_write_mode='auto' ):
    """ Add the readgroups in the set readgroups to the header of origin_bam_path.

    This function first loads the header of the origin to memory.
    The supplied readgroups are added to this header.
    The new header is then exported to a SAM file. The SAM file is then
    concatenated to the original bam file.

    Args:

        origin_bam_path(str) : path to bam file to which to add readgroups to header

        readgroups_in(set/dict) : set or dictionary which contains read groups. The dictionary should have the format { read_group_id (str)
                { 'ID': ID, 'LB':library,
                'PL':platform,
                'SM':sampleLib,
                'PU':readGroup }

        target_bam_path(str) : path to write bam file including the readgrouped header to. When not supplied the output is written to the input bam file

    """

    if target_bam_path is None:
        target_bam_path = origin_bam_path

    # Create a read group dictionary
    if type(readgroups_in) is set:
        readGroupsDict={}
        for readGroup in readgroups_in:
            flowCell,lane,sampleLib = readGroup.split('.')
            try:
                library,_ = sampleLib.rsplit('_',1)
            except Exception as e:
                # the library is not part of the sample name:
                library = 'undefinedLibrary'
            readGroupsDict[readGroup] = {
                    'ID':readGroup,
                    'LB':library,
                    'PL':'ILLUMINA',
                    'SM':sampleLib,
                    'PU':readGroup}
    elif type(readgroups_in) is dict:
        readGroupsDict = readgroups_in
    else:
        raise ValueError("supply a set or dict for readgroups_in")

    # Write the re-headered bam file to this path
    complete_temp_path = origin_bam_path.replace('.bam','')+'.rehead.bam'


    # When header_write_mode is auto, when samtools is available, samtools will be used, otherwise pysam
    if header_write_mode=='auto':
        if which('samtools') is None:
            header_write_mode='pysam'
        else:
            header_write_mode='samtools'

    if header_write_mode=='pysam':
        with pysam.AlignmentFile(origin_bam_path, "rb") as origin:
            header = origin.header.copy()
            hCopy = header.to_dict()
            hCopy['RG'] = list(readGroupsDict.values())
            with pysam.AlignmentFile(complete_temp_path, "wb", header=hCopy) as out:
                for read in origin:
                    out.write(read)

        os.rename(complete_temp_path,target_bam_path)

    elif header_write_mode=='samtools':
        with pysam.AlignmentFile(origin_bam_path, "rb") as origin:
            header = origin.header.copy()

            # Write the new header to this sam file:
            headerSamFilePath = origin_bam_path.replace('.bam','')+'.header.sam'

            hCopy = header.to_dict()
            hCopy['RG'] = list(readGroupsDict.values())

            # Write the sam file with the complete header:
            with pysam.AlignmentFile(headerSamFilePath,'w',header=hCopy):
                pass

        # Concatenate and remove origin
        rehead_cmd = f"""{{ cat {headerSamFilePath}; samtools view {origin_bam_path}; }} | samtools view -b > {complete_temp_path} ;
                mv {complete_temp_path } {target_bam_path};rm {headerSamFilePath};
                """
        os.system(rehead_cmd)
    else:
        raise ValueError('header_write_mode should be either, auto, pysam or samtools')
