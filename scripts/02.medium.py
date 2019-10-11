#!/usr/bin/env python3
"""
Filter all funfams by:
    - high dops (>70)
    - at least 5 sequences

Then include funfams with the following characteristics:
    - lots of sequences
    - not many sequences
    - with structural reps
    - without structural reps

Within all those groups, select:
    - variety of unique folds

"""

# core
import json
import logging
import os
import os.path as path

# pip
import pandas as pd

# local
from cathbaker import loaders, models 

# params
MIN_DOPS_SCORE = 70
MIN_SEQUENCES = 5
LOW_SEQUENCE_COUNT = 50
HIGH_SEQUENCE_COUNT = 5000

# filenames
DATASET_DIR = path.relpath(path.join(path.dirname(__file__), '..', 'dataset', '02.medium'))
FF_ALN_DIR = path.join(DATASET_DIR, 'alignments')
FF_ALN_DATAFILE = path.join(DATASET_DIR, '02.medium.{category}.tsv')

# logging
logging.basicConfig(level='INFO', format='%(message)s')
LOG = logging.getLogger(__name__)

def log_title(msg, *args):
    """Formats 'title' log entries"""
    LOG.info('')
    LOG.info(msg, *args)
    LOG.info('-' * len(msg))

def log_kv(key, value):
    """Formats key/value log entries"""
    LOG.info('{:<60s} {}'.format(key, value))

def main():
    """Generates the benchmark datasets"""

    if not path.exists(FF_ALN_DIR):
        LOG.info(f"Creating alignment directory '{FF_ALN_DIR}'")
        os.makedirs(FF_ALN_DIR)

    log_title('LOADING FUNFAM DATA')
    all_funfams = loaders.load_all_funfam_info(cath_version='v4_2_0')

    # for ff in all_funfams[:3]:
    #     LOG.info("ff: %s", ff)

    # unique superfamily ids
    all_sfam_ids = {ff.superfamily_id for ff in all_funfams}
    all_top_ids = {ff.topology_id for ff in all_funfams}
    LOG.info('Found %s funfams from %s unique superfamilies and %s unique folds', 
        len(all_funfams), len(all_sfam_ids), len(all_top_ids))

    # apply global filters
    log_title('APPLYING GLOBAL FILTERS')
    filtered_funfams = [ff for ff in all_funfams if ff.seed_dops_score >= MIN_DOPS_SCORE and ff.num_members_in_seed_aln > MIN_SEQUENCES]
    log_kv(f'High sequence diversity (seed_dops_score >= {MIN_DOPS_SCORE}):', 
        f'{len(filtered_funfams)} (from {len(all_funfams)})')

    # create categories (high seqs, low seqs, with structural reps, without structural reps)
    log_title('CREATING DATASET CATEGORIES')
    low_seq_ffs = [ff for ff in filtered_funfams if ff.num_members_in_seed_aln < LOW_SEQUENCE_COUNT]
    log_kv(f'Small Funfams (num_members < {LOW_SEQUENCE_COUNT}):',
        f'{len(low_seq_ffs)} (from {len(filtered_funfams)})')

    high_seq_ffs = [ff for ff in filtered_funfams if ff.num_members_in_seed_aln >= HIGH_SEQUENCE_COUNT]
    log_kv(f'Large Funfams (num_members >= {HIGH_SEQUENCE_COUNT}):', 
        f'{len(high_seq_ffs)} (from {len(filtered_funfams)})')

    with_str_ffs = [ff for ff in filtered_funfams if ff.rep_source_id == 'cath']
    log_kv(f'Funfams with structure (rep_source_id="cath"):',
        f'{len(with_str_ffs)} (from {len(filtered_funfams)})')

    without_str_ffs = [ff for ff in filtered_funfams if ff.rep_source_id == 'uniprot']
    log_kv(f'Funfams without structure (rep_source_id="uniprot"):',
        f'{len(without_str_ffs)} (from {len(filtered_funfams)})')

    log_title('WRITING FUNFAM ALIGNMENTS')

    found_topologies = set()

    df_all = None
    for dataset_title_value in [
            ['low sequences', 'low_seq', low_seq_ffs],
            ['high sequences', 'high_seq', high_seq_ffs],
            ['structure', 'with_str', with_str_ffs],
            ['no structure', 'without_str', without_str_ffs],
        ]:
        title, category, ffs = dataset_title_value

        log_title(f'Working on {len(ffs)} funfams with {title} ...')
        ff_by_topology = {}
        for ff in ffs:
            if ff.topology_id not in found_topologies:
                found_topologies.add(ff.topology_id)
                ff_by_topology[ff.topology_id] = ff
        topology_ffs = list(ff_by_topology.values())
        log_kv('Funfams with unique topologies:', 
            f'{len(topology_ffs)} (from {len(all_top_ids)})')

        df = retrieve_funfam_alignments(topology_ffs)

        # add 'category' column
        df['category'] = category

        # datafile = FF_ALN_DATAFILE.format(category=category)
        # LOG.info(f"Writing TSV datafile '{datafile}'")
        # df.to_csv(datafile, sep='\t', index=False)

        if df_all is None:
            df_all = df
        else:
            df_all = pd.concat([df_all, df], ignore_index=True)

    log_title('WRITING DATASETS')
    log_kv(f'Total number of funfams', f'{len(df_all.index)} (from {len(all_funfams)})' )
    datafile = FF_ALN_DATAFILE.format(category='all')
    LOG.info(f"Writing TSV datafile '{datafile}'")
    df_all.to_csv(datafile, sep='\t', index=False)

    LOG.info('DONE')


def retrieve_funfam_alignments(funfams, dirpath=FF_ALN_DIR):
    """
    Retrieves funfam alignments from API
    
    Returns:
        df (pandas.DataFrame): dataframe with a summary of aln data
    """

    colnames = ['ff_id', 'sfam_id', 'seq_count', 'dops_score', 'gap_per']
    rows = []
    for ff in funfams:
        aln = loaders.get_funfam_alignment(**(ff.__dict__))
        aln_meta = aln.get_meta_summary()

        ff_id = f'{ff.superfamily_id}-ff-{ff.funfam_number}'
        ff_fname = f'cath.{ff.cath_version}.{ff.superfamily_id}-ff-{ff.funfam_number}.seed.sto'

        if aln_meta.dops_score < MIN_DOPS_SCORE:
            LOG.warning(f'WARNING: {ff_id}: alignment has dops_score={aln_meta.dops_score} (API reported: {ff.seed_dops_score}) (SKIPPING)')
            continue

        if aln_meta.seq_count is not ff.num_members_in_seed_aln:
            LOG.warning(f'WARNING: {ff_id}: alignment has seq_count={aln_meta.seq_count} (API reported: {ff.num_members_in_seed_aln})')

        local_file = path.join(dirpath, ff_fname)
        LOG.info(f"Saving Funfam alignment {ff_id}")
        aln.write_sto(local_file)
        
        gap_per = 100 * aln.total_gap_positions / aln.total_positions

        row = [ff_id, ff.superfamily_id, aln_meta.seq_count, aln_meta.dops_score, gap_per]
        rows.extend([row])

    df = pd.DataFrame(data=rows, columns=colnames)

    return df


if __name__ == '__main__':
    main()