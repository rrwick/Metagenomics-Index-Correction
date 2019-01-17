#!/usr/bin/env python3
"""
This script will dereplicate GTDB assemblies to a user-specified threshold,
copying the dereplicated assemblies to a new directory.

Usage:
    dereplicate_assemblies.py --threshold 0.005 assemblies derep_assemblies bac_and_arc_taxonomy_r86.tsv
"""

import argparse
import collections
import gzip
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import textwrap


def get_arguments():
    parser = argparse.ArgumentParser(description='Cluster assemblies in each taxon')
    parser.add_argument('in_dir', type=str,
                        help='Directory containing all GTDB assemblies')
    parser.add_argument('out_dir', type=str,
                        help='Directory where dereplicated assemblies will be copied')
    parser.add_argument('tax_file', type=str,
                        help='GTDB taxonomy file')
    parser.add_argument('--threshold', type=float, default=0.002,
                        help='Mash distance clustering threshold')
    parser.add_argument('--threads', type=int, default=16,
                        help='Number of threads (for Mash)')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    all_assemblies = find_all_assemblies(args.in_dir)
    os.makedirs(args.out_dir, exist_ok=True)
    classifications = load_classifications(args.tax_file)
    for taxon in sorted(classifications.keys()):
        process_one_taxon(taxon, classifications[taxon], all_assemblies, args.out_dir, args.threads, args.threshold)
    print()


def load_classifications(tax_file):
    classifications = collections.defaultdict(list)
    with open(tax_file, 'rt') as tax:
        for line in tax:
            parts = line.strip().split('\t')
            accession = parts[0]
            if accession.startswith('RS_') or accession.startswith('GB_'):
                accession = accession[3:]
            taxon = parts[1]
            classifications[taxon].append(accession)
    return classifications


def process_one_taxon(classification, accessions, all_assemblies, out_dir, threads, threshold):
    accessions = sorted(accessions)
    print()
    print(classification)
    acc_to_assemblies = find_assemblies_for_accessions(accessions, all_assemblies)
    if len(acc_to_assemblies) == 0:
        return
    if len(acc_to_assemblies) == 1:
        only_assembly = list(acc_to_assemblies.values())[0]
        print('Only one assembly for this species, copying to output directory:')
        print('    {} -> {}'.format(only_assembly, out_dir))
        shutil.copy(only_assembly, out_dir)
    else:
        print('{:,} assemblies for this species, clustering to dereplicate.'.format(len(acc_to_assemblies)))
        derep_assemblies = dereplicate(acc_to_assemblies, threads, threshold)
        print('Copying dereplicated assemblies to output directory:')
        for assembly in derep_assemblies:
            print('    {} -> {}'.format(assembly, out_dir))
            shutil.copy(assembly, out_dir)


def dereplicate(acc_to_assemblies, threads, threshold):
    all_assemblies = sorted(acc_to_assemblies.values())
    derep_assemblies = []
    with tempfile.TemporaryDirectory() as temp_dir:
        mash_sketch = build_mash_sketch(all_assemblies, threads, temp_dir)
        pairwise_distances = pairwise_mash_distances(mash_sketch, threads)
        assemblies, graph = create_graph_from_distances(pairwise_distances, threshold)
        clusters = cluster_assemblies(assemblies, graph)
        for num, assemblies in clusters.items():
            noun = 'assembly' if len(assemblies) == 1 else 'assemblies'
            print('cluster {:,} ({:,} {}): '.format(num, len(assemblies), noun), end='')
            if len(assemblies) == 1:
                representative = assemblies[0]
                print(representative)
            else:
                n50, representative = sorted([(get_assembly_n50(a), a) for a in assemblies])[-1]
                print('{} (N50 = {:,})'.format(representative, n50))
            derep_assemblies.append(representative)
    return derep_assemblies


def build_mash_sketch(assemblies, threads, temp_dir):
    mash_command = ['mash', 'sketch', '-p', str(threads), '-o', temp_dir + '/mash',
                    '-s', '10000'] + assemblies
    subprocess.run(mash_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return temp_dir + '/mash.msh'


def pairwise_mash_distances(mash_sketch, threads):
    mash_command = ['mash', 'dist', '-p', str(threads), mash_sketch, mash_sketch]
    mash_out = subprocess.run(mash_command, stdout=subprocess.PIPE).stdout.decode()
    return mash_out.splitlines()


def find_all_assemblies(in_dir):
    print('\nLooking for files in {}:'.format(in_dir))
    all_assemblies = [str(x) for x in sorted(pathlib.Path(in_dir).glob('**/*'))
                      if x.is_file()]
    print('found {:,} files'.format(len(all_assemblies)))
    return all_assemblies


def find_assemblies_for_accessions(accessions, all_assemblies):
    acc_to_assemblies = {}
    found_count, total_count = 0, 0
    not_found = []
    for accession in accessions:
        total_count += 1
        assembly_filename = get_assembly_filename(accession, all_assemblies)
        if assembly_filename is not None:
            found_count += 1
            acc_to_assemblies[accession] = assembly_filename
        else:
            not_found.append(accession)
        print('\r{:,} / {:,} assemblies found'.format(found_count, total_count),
              end='', flush=True)
    print()
    if not_found:
        print('    failed to find assemblies for the following accessions:')
        wrapper = textwrap.TextWrapper(initial_indent='    ', subsequent_indent='    ', width=100)
        print(wrapper.fill(', '.join(not_found)))
    return acc_to_assemblies


def get_assembly_filename(accession, all_assemblies):
    if accession.startswith('GCF_') or accession.startswith('GCA_'):
        accession = accession.split('.')[0]
        assert len(accession) == 13
    accession_dot = accession + '.'
    accession_under = accession + '_'
    assembly_filenames = [x for x in all_assemblies
                          if x.rpartition('/')[2].startswith(accession_dot)
                          or x.rpartition('/')[2].startswith(accession_under)]
    if len(assembly_filenames) == 0:
        return None
    elif len(assembly_filenames) > 1:
        sys.exit('\nError: ambiguous assembly filenames, accession={}, '
                 'filenames={}'.format(accession, assembly_filenames))
    return assembly_filenames[0]


def create_graph_from_distances(pairwise_distances, threshold):
    """
    Builds an undirected graph where nodes are assemblies and edges connect assemblies which have
    a pairwise Mash distance below the threshold.
    """
    print('Loading distances...', end='', flush=True)
    assemblies = set()
    graph = collections.defaultdict(set)
    all_connections = collections.defaultdict(set)
    count = 0
    for line in pairwise_distances:
        parts = line.split('\t')
        assembly_1 = parts[0]
        assembly_2 = parts[1]
        distance = float(parts[2])
        assemblies.add(assembly_1)
        assemblies.add(assembly_2)
        if assembly_1 == assembly_2:
            continue
        all_connections[assembly_1].add(assembly_2)
        all_connections[assembly_2].add(assembly_1)
        if distance < threshold:
            graph[assembly_1].add(assembly_2)
            graph[assembly_2].add(assembly_1)
        count += 1
    print(' done ({})'.format(count))
    assemblies = sorted(assemblies)
    assembly_count = len(assemblies)    
    for assembly in assemblies: # sanity check: make sure we have all the connections
        assert len(all_connections[assembly]) == assembly_count - 1
    return assemblies, graph


def cluster_assemblies(assemblies, graph):
    visited = set()
    i = 0
    clusters = {}
    for assembly in assemblies:
        if assembly in visited:
            continue
        i += 1
        connected = dfs(graph, assembly)
        clusters[i] = sorted(connected)
        visited |= connected
    return clusters


def dfs(graph, start):
    visited, stack = set(), [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(graph[vertex] - visited)
    return visited


def get_assembly_n50(filename):
    contig_lengths = sorted(get_contig_lengths(filename), reverse=True)
    total_length = sum(contig_lengths)
    target_length = total_length * 0.5
    length_so_far = 0
    for contig_length in contig_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            return contig_length
    return 0


def get_contig_lengths(filename):
    lengths = []
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    lengths.append(len(sequence))
                    sequence = ''
                name = line[1:].split()[0]
            else:
                sequence += line
        if name:
            lengths.append(len(sequence))
    return lengths


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


if __name__ == '__main__':
    main()
