#!/usr/bin/env python3
"""
This script will convert a GTDB taxonomy into the style of an NCBI taxonomy, appropriate for usage
when building a Centrifuge or Kraken database.
"""

import argparse
import gzip
import itertools
import os
import pathlib
import textwrap
import sys
from Bio import Phylo


def get_arguments():
    parser = argparse.ArgumentParser(description='Build an NCBI-style taxonomy from GTDB')

    input_args = parser.add_argument_group('Input')
    input_args.add_argument('--gtdb', type=str, required=True,
                            help='GTDB bac_taxonomy_r*.tsv file (required)')
    input_args.add_argument('--assemblies', type=str,
                            help='directory containing GTDB assembly files (required to make '
                                 'contig-to-taxid conversion file, concatenated FASTA or Kraken '
                                 'directory)')
    input_args.add_argument('--tree', type=str,
                            help='If provided, this script will only use assemblies that are in '
                                 'the tree (for building with the dereplicated genomes)')

    output_args = parser.add_argument_group('Output')
    output_args.add_argument('--nodes', type=str, default='ex.tree',
                             help='Filename of nodes file (default: ex.tree)')
    output_args.add_argument('--names', type=str, default='ex.name',
                             help='Filename of names file (default: ex.name)')
    output_args.add_argument('--conversion', type=str,
                             help='Filename of contig name to tax id conversion table')
    output_args.add_argument('--cat_fasta', type=str,
                             help='Filename of concatenated FASTA')
    output_args.add_argument('--kraken_dir', type=str,
                             help='Directory to put Kraken-formatted assemblies')
    return parser.parse_args()


def main():
    args = get_arguments()
    check_args(args)

    all_taxa, parents, accession_to_species = load_taxa(args.gtdb)
    id_to_taxon, taxon_to_id, max_id = set_tax_ids(all_taxa)
    not_unique_names = get_not_unique_names(id_to_taxon, max_id)
    write_nodes_file(args.nodes, id_to_taxon, taxon_to_id, parents, max_id)
    write_names_file(args.names, id_to_taxon, not_unique_names, max_id)

    if args.assemblies is not None:
        accessions = sorted(accession_to_species.keys())
        all_assemblies = [str(x) for x in sorted(pathlib.Path(args.assemblies).glob('**/*'))
                          if x.is_file()]
        if args.tree is not None:
            whitelist = load_tree(args.tree)
            accessions = [x for x in accessions if x in whitelist]
        acc_to_assemblies = find_assemblies_for_accessions(accessions, all_assemblies)

        if args.conversion is not None:
            write_conversion_file(args.conversion, accession_to_species, taxon_to_id,
                                  acc_to_assemblies, accessions)
        if args.cat_fasta is not None:
            make_cat_fasta(args.cat_fasta, acc_to_assemblies, accessions)
        if args.kraken_dir is not None:
            make_kraken_dir(args.kraken_dir, accession_to_species, taxon_to_id, acc_to_assemblies,
                            accessions)
    print()


def check_args(args):
    if args.assemblies is not None and not pathlib.Path(args.assemblies).is_dir():
        sys.exit('Error: {} is not a directory'.format(args.assemblies))
    if args.kraken_dir is not None:
        kraken_dir = pathlib.Path(args.kraken_dir)
        if kraken_dir.is_file():
            sys.exit('Error: {} is a file (must be a directory)'.format(kraken_dir))
        if kraken_dir.is_dir():
            if len(list(kraken_dir.iterdir())) > 0:
                sys.exit('Error: {} is not empty'.format(kraken_dir))
        else:
            os.makedirs(str(kraken_dir))


def load_taxa(gtdb_taxonomy_filename):
    print('\nLoading GTDB taxonomy:')
    all_taxa = set()
    parents = {'': ''}
    accession_to_species = {}
    unknown_taxa_counter = itertools.count(start=0)
    with open(gtdb_taxonomy_filename, 'rt') as tax_file:
        for line in tax_file:
            parts = line.strip().split('\t')
            assembly_accession = parts[0]
            if assembly_accession.startswith('RS_') or assembly_accession.startswith('GB_'):
                assembly_accession = assembly_accession[3:]
            tax_levels = parts[1].split(';')
            for i, taxon in enumerate(tax_levels):
                if len(taxon) == 3:
                    assert taxon == 's__'
                    genus = tax_levels[i-1]
                    assert genus.startswith('g__')
                    new_species_name = 's__{} unknown_{}'.format(genus[3:],
                                                                 next(unknown_taxa_counter))
                    taxon = new_species_name
                if taxon.startswith('s__'):
                    accession_to_species[assembly_accession] = taxon
                parent = '' if i == 0 else tax_levels[i-1]
                if taxon in parents:
                    if parents[taxon] != parent:
                        sys.exit('Error: duplicate taxon name found ({})'.format(taxon))
                else:
                    parents[taxon] = parent
                all_taxa.add(taxon)
    print('    {:,} taxa found'.format(len(all_taxa)))
    return all_taxa, parents, accession_to_species


def set_tax_ids(all_taxa):
    print('\nAssigning tax IDs:')
    id_to_taxon = {1: ''}
    taxon_to_id = {'': 1}
    tax_id_counter = itertools.count(start=2)
    for level in ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']:
        level_taxa = sorted(x for x in all_taxa if x.startswith(level))
        level_name = {'d__': 'domains: ', 'p__': 'phyla:   ', 'c__': 'classes: ',
                      'o__': 'orders:  ', 'f__': 'families:', 'g__': 'genera:  ',
                      's__': 'species: '}[level]
        for taxon in level_taxa:
            tax_id = next(tax_id_counter)
            id_to_taxon[tax_id] = taxon
            if taxon in taxon_to_id:
                assert taxon_to_id[taxon] == tax_id
            else:
                taxon_to_id[taxon] = tax_id
        print('    {}{:>7,}'.format(level_name, len(level_taxa)))
    max_id = next(tax_id_counter)
    return id_to_taxon, taxon_to_id, max_id


def get_rank_from_taxon_name(taxon):
    if taxon == '':
        return 'no rank'
    elif taxon.startswith('d__'):
        return 'domain'
    elif taxon.startswith('p__'):
        return 'phylum'
    elif taxon.startswith('c__'):
        return 'class'
    elif taxon.startswith('o__'):
        return 'order'
    elif taxon.startswith('f__'):
        return 'family'
    elif taxon.startswith('g__'):
        return 'genus'
    elif taxon.startswith('s__'):
        return 'species'


def get_embl_code(taxon):
    if taxon.startswith('s__'):
        g_letter = taxon[3].upper()
        s_letter = taxon.split(' ')[1][0].upper()
        return g_letter + s_letter
    else:
        return ''


def get_not_unique_names(id_to_taxon, max_id):
    used_names = set()
    not_unique_names = set()
    for tax_id in range(2, max_id):
        taxon = id_to_taxon[tax_id]
        taxon_name = taxon[3:]
        if taxon_name in used_names:
            not_unique_names.add(taxon_name)
        else:
            used_names.add(taxon_name)
    return not_unique_names


def write_nodes_file(nodes_out_file, id_to_taxon, taxon_to_id, parents, max_id):
    print('\nWriting taxonomy nodes to {}'.format(nodes_out_file))
    with open(nodes_out_file, 'wt') as nodes:
        for tax_id in range(1, max_id):
            taxon = id_to_taxon[tax_id]
            parent_tax_id = taxon_to_id[parents[taxon]]
            rank = get_rank_from_taxon_name(taxon)
            embl_code = get_embl_code(taxon)
            division_id = 8 if taxon == '' else 0
            inherited_div_flag = 0 if taxon == '' else 1
            genetic_code_id = 1 if taxon == '' else 11
            inherited_gc_flag = 0 if taxon == '' else 1
            mitochondrial_genetic_code_id = 0
            inherited_mgc_flag = 0 if taxon == '' else 1
            genbank_hidden_flag = 1
            hidden_subtree_root_flag = 0
            comments = taxon
            nodes.write('{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t'
                        '{}\t|\t{}\t|\t{}'
                        '\t|\n'.format(tax_id, parent_tax_id, rank, embl_code, division_id,
                                       inherited_div_flag, genetic_code_id, inherited_gc_flag,
                                       mitochondrial_genetic_code_id, inherited_mgc_flag,
                                       genbank_hidden_flag, hidden_subtree_root_flag, comments))


def write_names_file(names_out_file, id_to_taxon, not_unique_names, max_id):
    print('Writing taxonomy names to {}'.format(names_out_file))
    with open(names_out_file, 'wt') as names:
        names.write('1\t|\tall\t|\t\t|\tsynonym\t|\n')
        for tax_id in range(1, max_id):
            taxon = id_to_taxon[tax_id]
            if taxon == '':
                taxon_name = 'root'
            else:
                taxon_name = taxon[3:]
            if taxon_name in not_unique_names:
                unique_name = taxon
            else:
                unique_name = ''
            names.write('{}\t|\t{}\t|\t{}\t|\tscientific name\t|\n'.format(tax_id, taxon_name,
                                                                           unique_name))


def load_tree(tree_filename):
    if tree_filename is not None:
        print('\nLoading tree:')
        trees = Phylo.parse(tree_filename, 'newick')
        whitelist_assemblies = set()
        for tree in trees:
            whitelist_assemblies |= set(get_tip_names(tree.root))
        print('    found {:,} assemblies in {}'.format(len(whitelist_assemblies), tree_filename))
    else:
        whitelist_assemblies = None
    return whitelist_assemblies


def find_assemblies_for_accessions(accessions, all_assemblies):
    print('\nSearching for an assembly for each accession:')
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
        print('\r    {:,} / {:,} assemblies '
              'found'.format(found_count, total_count), end='', flush=True)
    print()
    if not_found:
        print('    failed to find assemblies for the following accessions:')
        wrapper = textwrap.TextWrapper(initial_indent='    ', subsequent_indent='    ', width=100)
        print(wrapper.fill(', '.join(not_found)))
    return acc_to_assemblies


def write_conversion_file(conversion_filename, accession_to_species, taxon_to_id, acc_to_assemblies,
                          accessions):
    print('\nWriting conversion table:')
    found_count, total_count = 0, len(acc_to_assemblies)
    with open(conversion_filename, 'wt') as conversion_file:
        for i, accession in enumerate(accessions):
            if accession in acc_to_assemblies:
                assembly_filename = acc_to_assemblies[accession]
                contig_names = load_contig_names(assembly_filename)
                tax_id = taxon_to_id[accession_to_species[accession]]
                for contig_name in contig_names:
                    conversion_file.write('{}_{}\t{}\n'.format(accession, contig_name, tax_id))
                found_count += 1
                print('\r    {:,} / {:,} assemblies'.format(found_count,
                                                            total_count), end='', flush=True)
    print()


def make_cat_fasta(cat_fasta_filename, acc_to_assemblies, accessions):
    print('\nWriting concatenated FASTA:')
    found_count, total_count = 0, len(acc_to_assemblies)
    with open(cat_fasta_filename, 'wt') as cat_fasta:
        for i, accession in enumerate(accessions):
            if accession in acc_to_assemblies:
                assembly_filename = acc_to_assemblies[accession]
                found_count += 1
                contigs = load_fasta(assembly_filename)
                for contig_name, contig_seq in contigs:
                    cat_fasta.write('>{}_{}\n'.format(accession, contig_name))
                    cat_fasta.write('{}\n'.format(contig_seq))
            print('\r    {:,} / {:,} assemblies'.format(found_count,
                                                        total_count), end='', flush=True)
    print()


def make_kraken_dir(kraken_dir, accession_to_species, taxon_to_id, acc_to_assemblies, accessions):
    print('\nMaking Kraken assembly directory:')
    kraken_dir = pathlib.Path(kraken_dir)
    found_count, total_count = 0, len(acc_to_assemblies)
    for accession in accessions:
        if accession in acc_to_assemblies:
            assembly_filename = acc_to_assemblies[accession]
            found_count += 1
            contigs = load_fasta(assembly_filename)
            tax_id = taxon_to_id[accession_to_species[accession]]
            new_filename = kraken_dir / (accession + '.fa')
            with open(str(new_filename), 'wt') as kraken_fasta:
                for contig_name, contig_seq in contigs:
                    kraken_fasta.write('>{}_{}|kraken:taxid|{}\n'.format(accession, contig_name,
                                                                         tax_id))
                    kraken_fasta.write('{}\n'.format(contig_seq))
        print('\r    {:,} / {:,} assemblies'.format(found_count,
                                                    total_count), end='', flush=True)
    print()


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


def get_compression_type(filename):
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    with open(filename, 'rb') as unknown_file:
        file_start = unknown_file.read(max_len)
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_function(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def load_contig_names(filename):
    try:
        contig_names = set()
        open_func = get_open_function(filename)
        with open_func(filename, 'rt') as fasta_file:
            for line in fasta_file:
                line = line.strip()
                if not line:
                    continue
                if line[0] == '>':  # Header line = start of new contig
                    contig_name = line[1:].split()[0]
                    if contig_name in contig_names:
                        sys.exit('Error: duplicate contig names in {}'.format(filename))
                    contig_names.add(contig_name)
        return sorted(contig_names)
    except EOFError:
        print('\n    Warning: {} seems to be corrupt'.format(filename))
        return []


def load_fasta(filename):
    try:
        fasta_seqs = []
        open_func = get_open_function(filename)
        with open_func(filename, 'rt') as fasta_file:
            name = ''
            sequence = []
            for line in fasta_file:
                line = line.strip()
                if not line:
                    continue
                if line[0] == '>':  # Header line = start of new contig
                    if name:
                        contig_name = name.split()[0]
                        fasta_seqs.append((contig_name, ''.join(sequence)))
                        sequence = []
                    name = line[1:]
                else:
                    sequence.append(line)
            if name:
                contig_name = name.split()[0]
                fasta_seqs.append((contig_name, ''.join(sequence)))
        return fasta_seqs
    except EOFError:
        print('\n    Warning: {} seems to be corrupt'.format(filename))
        return []


def get_tip_names(clade):
    tip_names = []
    if clade.name is not None and len(clade.clades) == 0:
        tip_names.append(clade.name)
    for child in clade:
        tip_names += get_tip_names(child)
    tip_names = [x[3:] if x.startswith('RS_') or x.startswith('GB_') else x for x in tip_names]
    return sorted(tip_names)


if __name__ == '__main__':
    main()
