#!/usr/bin/env python3

import collections
import sys


def main():
    try:
        sam_filename = sys.argv[1]
        tree_filename = sys.argv[2]
    except IndexError:
        print_header()
        quit()

    tax_id_to_parent = load_tax_id_to_parent(tree_filename)
    tax_id_to_rank = load_tax_id_to_rank(tree_filename, tax_id_to_parent)
    count_per_rank = collections.defaultdict(int)

    current_read_name, current_tax_ids = '', set()
    with open(sam_filename, 'rt') as sam_file:
        for line in sam_file:
            parts = line.strip().split('\t')
            read_name = parts[0]
            if read_name == 'readID':  # header
                continue
            tax_id = int(parts[2])
            if read_name != current_read_name:
                add_rank_count(current_read_name, count_per_rank, current_tax_ids, tax_id_to_rank, tax_id_to_parent)
                current_read_name, current_tax_ids = read_name, set()
            current_tax_ids.add(tax_id)
    add_rank_count(read_name, count_per_rank, current_tax_ids, tax_id_to_rank, tax_id_to_parent)

    print_result(sam_filename, count_per_rank)



def load_tax_id_to_parent(tree_filename):
    tax_id_to_parent = {}
    with open(tree_filename, 'rt') as tree_file:
        for line in tree_file:
            parts = line.strip().split('\t')
            tax_id = int(parts[0])
            parent = int(parts[2])
            tax_id_to_parent[tax_id] = parent
    return tax_id_to_parent


def load_tax_id_to_rank(tree_filename, tax_id_to_parent):
    tax_id_to_rank = {0: 'unclassified'}
    acceptable_ranks = {'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'}

    # The first time we go through the tax IDs, we save any which has an acceptable rank.
    with open(tree_filename, 'rt') as tree_file:
        for line in tree_file:
            parts = line.strip().split('\t')
            tax_id, rank = int(parts[0]), parts[4]
            if rank in acceptable_ranks:
                tax_id_to_rank[tax_id] = rank
            elif tax_id == 1:  # special case for the root node
                assert tax_id == int(parts[2])
                tax_id_to_rank[tax_id] = 'root'

    # Now we go through a second time to deal with tax IDs that didn't get an acceptable rank the
    # first time.
    with open(tree_filename, 'rt') as tree_file:
        for line in tree_file:
            parts = line.strip().split('\t')
            tax_id, rank = int(parts[0]), parts[4]
            if tax_id in tax_id_to_rank:
                continue
            assert rank not in acceptable_ranks
            ancestors = get_all_ancestors(tax_id, tax_id_to_parent)
            for ancestor in ancestors:
                if ancestor in tax_id_to_rank:
                    rank = tax_id_to_rank[ancestor]
                    tax_id_to_rank[tax_id] = rank
                    break
            assert tax_id in tax_id_to_rank

    return tax_id_to_rank


def add_rank_count(read_name, count_per_rank, tax_ids, tax_id_to_rank, tax_id_to_parent):
    if len(tax_ids) == 0:
        return
    print('{}\t{}'.format(read_name, tax_ids), file=sys.stderr, end='\t')
    if len(tax_ids) == 1:
        (tax_id,) = tax_ids
    else:
        tax_ids.discard(0)
        tax_id = find_lca(tax_ids, tax_id_to_parent)
    rank = tax_id_to_rank[tax_id]
    print('{}\t{}'.format(tax_id, rank), file=sys.stderr)
    # print('{}\t{}\t{}'.format(read_name, tax_id, rank))
    count_per_rank[rank] += 1


def find_lca(tax_ids, tax_id_to_parent):
    common_taxa = set()
    first_ancestor_list = None
    for tax_id in tax_ids:
        ancestors = get_all_ancestors(tax_id, tax_id_to_parent)
        if first_ancestor_list is None:
            first_ancestor_list = ancestors
        if not common_taxa:
            common_taxa = set(ancestors)
        else:
            common_taxa &= set(ancestors)
    best_choice = None
    for t in common_taxa:
        i = first_ancestor_list.index(t)
        if best_choice is None:
            best_choice = t, i
        else:
            if i < best_choice[1]:
                best_choice = t, i
    return best_choice[0]


def get_all_ancestors(tax_id, tax_id_to_parent):
    ancestors = [tax_id]
    while tax_id != 1:
        tax_id = tax_id_to_parent[tax_id]
        ancestors.append(tax_id)
    return ancestors


def print_header():
    header = 'file'
    for rank in ['unclassified', 'root', 'domain', 'phylum', 'class', 'order', 'family', 'genus',
                 'species']:
        header += '\t{}_count\t{}_percent'.format(rank, rank)
    print(header)


def print_result(sam_filename, count_per_rank):
    total = sum(count_per_rank.values())
    result = sam_filename
    for rank in ['unclassified', 'root', 'domain', 'phylum', 'class', 'order', 'family', 'genus',
                 'species']:
        count = count_per_rank[rank]
        percent = 100 * count / total
        result += '\t{}\t{:.4f}'.format(count, percent)
    print(result)


if __name__ == '__main__':
    main()
