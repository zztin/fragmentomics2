#!/usr/bin/env python3
import pysam
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import singlecellmultiomics.pyutils
import pysamiterators
import collections
import glob
import pickle
import pandas as pd
from colorama import Fore, Style
from singlecellmultiomics.bamProcessing.bamFunctions import sorted_bam_file, sort_and_index, get_reference_from_pysam_alignmentFile, add_readgroups_to_header, write_program_tag
from singlecellmultiomics.pyutils import meanOfCounter
import argparse

f'Please use an environment with python 3.6 or higher!'

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="""
Known variant locations extraction tool
""")
argparser.add_argument('bamfiles', nargs='+')
argparser.add_argument('-ssnv', help="sSNV bed file", type=str, required=True)
argparser.add_argument('-gsnv', help="gSNV bed file", type=str, required=True)
argparser.add_argument('-reference', help="reference fasta file", type=str, required=True)
argparser.add_argument('-head', type=int)
argparser.add_argument('-min_read_obs', type=int, default=2)


args=  argparser.parse_args()
WINDOW_RADIUS = 250

paths = args.bamfiles


def obtain_variant_statistics(
    alignment_file_paths,
    cell_obs, statistics,cell_call_data,
    reference,
    chromosome,
    ssnv_position,gsnv_position,haplotype_scores,
    WINDOW_RADIUS, out , min_read_obs,read_groups

    ):


    sSNV_ref_base = reference.fetch(chromosome,ssnv_position,ssnv_position+1)
    gSNV_ref_base = reference.fetch(chromosome,gsnv_position,gsnv_position+1)

    window_molecules = []

    cell_read_obs = collections.defaultdict( collections.Counter ) # sample -> tuple -> amount of reads
    for pathi,path in enumerate(alignment_file_paths):
        with pysam.AlignmentFile(path) as alignments:

            for molecule_id,molecule in enumerate(
                singlecellmultiomics.molecule.MoleculeIterator(alignments,
                    fragment_class_args={
                        'umi_hamming_distance':1,

                    },
                    molecule_class_args={
                        'reference':reference
                    },
                   moleculeClass=singlecellmultiomics.molecule.NlaIIIMolecule,
                   fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
                    start = ssnv_position-WINDOW_RADIUS,
                    end = ssnv_position+ WINDOW_RADIUS,
                   contig=chromosome
                )):

                # For every molecule obtain the consensus from which to extract the gSNV and sSNV:
                try:
                    consensus = molecule.get_consensus()
                except Exception as e:
                    if str(e) == 'Could not extract any safe data from molecule':
                        statistics[(chromosome, ssnv_position)]['R2_unmapped'][True]+=1
                    continue

                # Extract the gSNV and sSNV:
                ssnv_state = consensus.get((chromosome, ssnv_position))
                gsnv_state = consensus.get((chromosome, gsnv_position) )

                # Store all used molecules in the window for inspection:
                window_molecules.append((molecule,ssnv_state,gsnv_state))

                # If both the ssnv and gsnv are none there is no information we can use.
                if ssnv_state is None and gsnv_state is None:
                    continue

                # Store the observation
                # the amount of reads of evidence is len(molecule)
                cell_obs[(chromosome, ssnv_position)][molecule.get_sample()][(ssnv_state,gsnv_state)] += 1
                cell_read_obs[molecule.get_sample()][(ssnv_state,gsnv_state)] += len(molecule)

                # Store statistics
                statistics[(chromosome, ssnv_position)]['max_mapping_quality'][ molecule.get_max_mapping_qual() ] +=1
                statistics[(chromosome, ssnv_position)]['fragment_size'][molecule.get_safely_aligned_length()]+=1
                statistics[(chromosome, ssnv_position)]['ivt_dups'][ len(  molecule.get_rt_reactions() )] +=1
                statistics[(chromosome, ssnv_position)]['undigested'][ molecule.get_undigested_site_count()]+=1
                statistics[(chromosome, ssnv_position)]['reads'][ len(molecule) ] +=1
                statistics[(chromosome, ssnv_position)]['molecules'][1] +=1

                # Store alignment statistics:
                for operation, per_bp in molecule.get_alignment_stats().items():
                    statistics[(chromosome, ssnv_position)][operation][per_bp] +=1

                try:
                    statistics[(chromosome, ssnv_position)]['ssnv_ref_phred'][ molecule.get_mean_base_quality(chromosome, ssnv_position, sSNV_ref_base ) ] +=1
                except:
                    pass

                try:
                    statistics[(chromosome, ssnv_position)]['ssnv_alt_phred'][ molecule.get_mean_base_quality(chromosome, ssnv_position, not_base = sSNV_ref_base )]+=1
                except:
                    pass

                try:
                    statistics[(chromosome, ssnv_position)]['gsnv_ref_phred'][ molecule.get_mean_base_quality(chromosome, gsnv_position, gSNV_ref_base)]+=1
                except:
                    pass

                try:
                    statistics[(chromosome, ssnv_position)]['gsnv_any_alt_phred'][ molecule.get_mean_base_quality(chromosome, gsnv_position, not_base= gSNV_ref_base)]+=1
                except:
                    pass

    # After finishing iteration over all molecules assign genotypes
    chrom, pos = chromosome, ssnv_position #
    obs_for_cells = cell_obs[(chrom, pos)]

    sSNV_alt_base = None
    gSNV_alt_base = None

    genotype_obs = collections.Counter()
    complete_genotype_obs = collections.Counter()
    sSNV_obs_phased = collections.Counter()
    gSNV_obs_phased = collections.Counter()

    sSNV_obs = collections.Counter()
    gSNV_obs = collections.Counter()

    for cell, cell_data in obs_for_cells.items():
        for ssnv,gsnv in cell_data:
            genotype_obs[(ssnv,gsnv)]+=1

            gSNV_obs[gsnv] +=1
            sSNV_obs[ssnv] += 1

            if ssnv is not None and gsnv is not None:
                complete_genotype_obs[(ssnv,gsnv)]+=1
                gSNV_obs_phased[gsnv] +=1 # Only count these when the germline variant is detected
                sSNV_obs_phased[ssnv] += 1

    print(Style.BRIGHT + f'Genotype observations for variant {chrom}:{pos}' + Style.RESET_ALL)
    print('som\tgerm\tobs')
    for (ssnv,gsnv), obs  in complete_genotype_obs.most_common():
        print(f' {ssnv}\t{gsnv}\t{obs}')

    if len(complete_genotype_obs)<=2:
        print(f'not enough genotype observations for a variant call (<=2)')

    ### Conbase algorithm : ###
    #
    # determine if there is an alternative base in the first place
    # a fraction of the reads in a cell need to vote for a tuple,
    # this fraction is stored in the alpha parameter , or a minimum amount of reads, stored in the beta parameter
    # determine tp*, the alleles we expect observe
    # ϴ τ α γ κ λ ν ξ ρ ϕ
    α = 0.2 # minimum relative abundance of sSNV voting reads in single sample
    β = 10 # minimum amount of sSSNV reads in cell, or in total if α is exceeded
    γ = 0.9 # minimum amount of votes for sSNV
    ε = 2 #minimum amount of cells voting for sSNV
    ω = 0.9 # gsnv majority

    sSNV_votes = collections.Counter() # { sSNV_alt_base : votes }
    total_samples_which_voted = 0
    for sample, observed_tuples in cell_read_obs.items():
        # First we create a Counter just counting the amount of evidence per base for this sample :
        evidence_total_reads = collections.Counter()
        total_reads = 0
        for (sSNV_state,gSNV_state),reads in observed_tuples.most_common():
            if sSNV_state is None:
                continue
            evidence_total_reads[sSNV_state] += reads
            total_reads += reads

        ref_sSNV_reads = evidence_total_reads[sSNV_ref_base] # this is the amount of reads which contain evidence for the reference base
        votes_for_this_sample = set() # the alternative bases this sample votes for
        for sSNV_state, sSNV_supporting_reads in evidence_total_reads.most_common():
            # The reference base does not vote.
            if sSNV_state == sSNV_ref_base or sSNV_state is None:
                continue

            # check if at least alpha reads vote for the sSNV
            alpha_value = sSNV_supporting_reads/ref_sSNV_reads

            vote = (1 if
                    (alpha_value >= α and (sSNV_supporting_reads+ref_sSNV_reads) >= β ) or \
                    (alpha_value < α and ref_sSNV_reads >= β )
                    else 0
                )

            if vote:
                votes_for_this_sample.add(sSNV_state)
                sSNV_votes[sSNV_state] += 1
                total_samples_which_voted += 1

    # done voting.
    # the most probable variant base is at least 90% voted for (lambda parameter)
    # and at least ε cells need to vote for it
    statistics[(chromosome, ssnv_position)]['total_samples_voted'] =  total_samples_which_voted
    if total_samples_which_voted < ε:
        # We don't have enough votes
        return

    sSNV_alt_base, sSNV_alt_obs = sSNV_votes.most_common()[0]
    statistics[(chromosome, ssnv_position)]['sSNV_alt_vote_ratio'] =  (sSNV_alt_obs / total_samples_which_voted)
    if (sSNV_alt_obs / total_samples_which_voted) < γ:
        # The ratio of votes is below threshold
        return

    print(f'sSNV alt is {sSNV_alt_base}')

    ### Here the "Stats" part of Conbase ends ###
    #############################################

    # Now we determined the sSNV alt base,
    # now determine the linked gSNV
    gSNV_alt_base=None # Lazy not defined before
    for basecall, obs in gSNV_obs_phased.most_common():
        if basecall!=gSNV_ref_base:
            gSNV_alt_base = basecall
            break

    if sSNV_alt_base is None or gSNV_alt_base is None:
        # No phased alt base found ...
        return

    # Determine the phase (most common genotypes)
    sSNV_phase = None
    wt_allele_gSNV= None
    sSNV_phased_votes = sum( (obs
        for (sSNV_state, gSNV_state), obs
        in complete_genotype_obs.most_common()
        if sSNV_state == sSNV_alt_base
        ))

    for (sSNV_state, gSNV_state), obs in complete_genotype_obs.most_common():
        if sSNV_state == sSNV_alt_base:
            sSNV_phase = (sSNV_state, gSNV_state)
            phased_gSNV = gSNV_state
            if gSNV_state==gSNV_ref_base:
                # the reference allele is alt
                wt_allele_gSNV = gSNV_alt_base
                snv_allele_gSNV = gSNV_ref_base
            else:
                wt_allele_gSNV = gSNV_ref_base
                snv_allele_gSNV = gSNV_alt_base
                # the reference allele is ref
            break


        # Verify that at least ω votes are cast for the currently selected allele :
        if (sSNV_phased_votes / obs) < ω:
            return

        # The valid tuples are thus:
        uninformative_allele = (sSNV_ref_base, wt_allele_gSNV)
        informative_allele_wt = (sSNV_ref_base, snv_allele_gSNV)

        valid_tuples = [ sSNV_phase,  # mutated
                        informative_allele_wt, #wt
                        uninformative_allele


        # As we have umi's we just have a threshold for the least amount of reads
        # we want to observe for a molecule to be taken into account
        # Count how often we found valid and invalid genotypes
        valid = 0
        invalid = 0
        valid_var = 0
        invalid_var = 0
        for (ssnv,gsnv), tuple_obs in complete_genotype_obs.most_common():
            if ssnv == sSNV_alt_base: # variant:
                if (ssnv,gsnv) in valid_tuples:
                    valid_var+=tuple_obs
                else:
                    invalid_var+=tuple_obs

            if (ssnv,gsnv) in valid_tuples:
                valid += tuple_obs
            else:
                invalid +=tuple_obs

        phase_ratio = 0
        if valid_var+invalid_var > 0  :
            phase_ratio = valid_var/(valid_var+invalid_var)

        # Score Tuples with evidence for variant
        haplotype_scores[(chrom,pos)] = {
            'valid_tuples':valid,
            'invalid_tuples':invalid,
            'valid_var_tuples':valid_var,
            'invalid_var_tuples':invalid_var,
            'phasing_ratio' : phase_ratio,
            'gSNV_allelic_bias':gSNV_obs[gSNV_ref_base]/(gSNV_obs[gSNV_ref_base]+gSNV_obs[gSNV_alt_base])
        }

        print(f'Germline variant obs: {gSNV_ref_base} {gSNV_alt_base}')
        print(f'sSNV obs: {sSNV_ref_base} {sSNV_alt_base}')
        if sSNV_phase is not None:
            print(f'sSNV variant is phased with {phased_gSNV}')
        print(Style.BRIGHT +'Valid tuples:'+Style.RESET_ALL)
        for g,s in valid_tuples:
            print(f' {g}\t{s}')

        print(Style.BRIGHT +'Scores:'+Style.RESET_ALL)
        for name, obs in haplotype_scores[(chrom,pos)].items():
            print(f' {name}\t{obs}')

        # Create the cell call dictionary

        for cell, observed_tuples in cell_read_obs.items():
            total_reads = 0
            phased_variant_support_reads = 0
            unphased_variant_support_reads = 0
            variant_neg_support_reads = 0
            uninformative_reads = 0
            conflict_reads = 0
            for  (sSNV_state,gSNV_state),reads in observed_tuples:
                if sSNV_state is None:
                    continue
                total_reads+=reads
                if sSNV_state == sSNV_alt_base:
                    if gSNV_state == wt_allele_gSNV:
                        conflict_reads += reads
                    elif gSNV_state == snv_allele_gSNV:
                        phased_variant_support_reads += reads # reads containing the sSNV and gSNV as expected
                    elif gSNV_state is None:
                        unphased_variant_support_reads += reads  # reads containing sSNV but not overlapping with gSNV

                elif sSNV_state == sSNV_ref_base:
                    if gSNV_state == snv_allele_gSNV:
                        variant_neg_support_reads += reads # reads on informative allele where we found evidence of the sSNV not being present


            if conflict_reads/(total_reads) > 0.2:
                cell_call_data[(chrom, pos)][cell] =  -1 # invalid

            if (unphased_variant_support_reads+phased_variant_support_reads)/total_reads > 0.1:
                cell_call_data[(chrom, pos)][cell] = 1
                if unphased_variant_support_reads+phased_variant_support_reads>=3:
                    cell_call_data[(chrom, pos)][cell] = 10

            if (phased_variant_support_reads)/total_reads > 0.1:
                cell_call_data[(chrom, pos)][cell] = 2
                if phased_variant_support_reads>=3:
                    cell_call_data[(chrom, pos)][cell] = 20

            if variant_neg_support_reads/total_reads > 0.1:
                cell_call_data[(chrom, pos)][cell] += 0.1 # 0.1 for ref allele obs


        # Annotate every molecule...

        for molecule_id,(m,ssnv_state, gsnv_state) in enumerate(window_molecules):
            m.set_meta('mi',molecule_id)
            if gsnv_state is None:
                m.set_meta('gv','?')
            else:
                m.set_meta('gv',gsnv_state)

            if ssnv_state is None:
                m.set_meta('sv','?')
            else:
                m.set_meta('sv',ssnv_state)


            if ssnv_state is None:
                m.set_meta('VD','NO_SNV_OVERLAP')
                continue

            if gsnv_state is not None and not (ssnv_state,gsnv_state) in valid_tuples:
                m.set_meta('VD','INVALID_PHASE')
                continue
            if ssnv_state == sSNV_alt_base:
                m.set_meta('VD','SNV_ALT')
                continue

            if ssnv_state == sSNV_ref_base and gsnv_state==phased_gSNV:
                m.set_meta('VD','SNV_REF')
                continue
            if gsnv_state!=phased_gSNV:
                m.set_meta('VD','UNINFORMATIVE_ALLELE')
                continue

            m.set_meta('VD','REJECTED')
        # write

        for m,ssnv_state, gsnv_state in window_molecules:
            m.write_tags()
            m.write_pysam(out)

            # Update read groups
            for fragment in m:
                read_groups.add(fragment.get_read_group())


# Load probed variants
probed_variants = {}
with open(args.ssnv) as s, \
     open(args.gsnv) as g:
    for i,(ssnv_line,gsnv_line) in enumerate(zip(s, g)):
        if ssnv_line.startswith('track name'):
            continue
        chrom, snv_pos, _ = ssnv_line.strip().split()
        _, gsnv_pos,__ = gsnv_line.strip().split()
        snv_pos, gsnv_pos = int(snv_pos), int(gsnv_pos)
        probed_variants[ (chrom,snv_pos) ]  = gsnv_pos

reference = pysamiterators.CachedFasta(pysam.FastaFile(args.reference))


cell_obs = collections.defaultdict(lambda: collections.defaultdict( collections.Counter) )
statistics = collections.defaultdict(lambda: collections.defaultdict( collections.Counter) )
cell_call_data =collections.defaultdict(dict) #location->cell->haplotype
haplotype_scores = {}


read_groups = set() # Store unique read groups in this set
with sorted_bam_file('evidence.bam', origin_bam=pysam.AlignmentFile(paths[0]), read_groups=read_groups ) as out:

    for variant_index, ((chromosome, ssnv_position),potential_gsnv_position) in enumerate(probed_variants.items()):

        obtain_variant_statistics(
            alignment_file_paths = paths,
            cell_obs = cell_obs,
            cell_call_data=cell_call_data,
            statistics = statistics,
            reference=reference,
            chromosome = chromosome,
            ssnv_position = ssnv_position,
            gsnv_position = potential_gsnv_position,
            WINDOW_RADIUS=WINDOW_RADIUS,
            haplotype_scores=haplotype_scores,
            out=out,min_read_obs=args.min_read_obs,
            read_groups=read_groups
        )

        if args.head and (variant_index>args.head-1):
            print(f'Stopping at variant {variant_index+1} because head was supplied ')
            break


lambda_free_dict = {}
for key, stats in statistics.items():
    lambda_free_dict[key] = {

        'mean_clip_pbp' : meanOfCounter(stats['clips_per_bp']),
        'mean_ins_pbp' : meanOfCounter(stats['inserts_per_bp']),
        'mean_del_pbp' : meanOfCounter(stats['deletions_per_bp']),
        'mean_matches_pbp' : meanOfCounter(stats['matches_per_bp']),
        'mean_alt_mapping_per_read' : meanOfCounter(stats['alt_per_read']),

        'ssnv_ref_phred' : meanOfCounter(stats['ssnv_ref_phred']),
        'ssnv_alt_phred' : meanOfCounter(stats['ssnv_alt_phred']),
        'gsnv_ref_phred' : meanOfCounter(stats['gsnv_ref_phred']),
        'gsnv_alt_phred' : meanOfCounter(stats['gsnv_alt_phred']),

        'mean_max_mapping_quality' : meanOfCounter(stats['max_mapping_quality']),
        'mean_ivt_dups' : meanOfCounter(stats['ivt_dups']),
        'mean_undigested' : meanOfCounter(stats['undigested']),
        'R2_unmapped' : stats['R2_unmapped'][True],


        'mean_fragment_size' : meanOfCounter(stats['fragment_size']),
        'mean_reads' : meanOfCounter(stats['reads']),
        'total_reads' : sum( (amount*frequency for amount, frequency in stats['reads'].most_common() )),
        'total_molecules' : sum( (amount*frequency for amount, frequency in stats['molecules'].most_common() ))

    }



print('Writing final site table')
site_stats = pd.DataFrame( lambda_free_dict ).T.join(pd.DataFrame(haplotype_scores).T)
site_stats.to_pickle('site_stats.pickle.gz')
site_stats.to_csv('site_stats.csv')

print('Writing final cell table')
cell_call_df = pd.DataFrame(cell_call_data)
cell_call_df.to_pickle('cell_calls.pickle.gz')
cell_call_df.to_csv('cell_calls.csv')