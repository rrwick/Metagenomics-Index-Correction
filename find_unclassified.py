#!/usr/bin/env python3
"""
This script takes two Centrifuge SAM files as input. It finds the reads in the first SAM file
which are unclassified and then outputs the lines from the second SAM file which correspond to the
same reads.
"""

import argparse
import collections
import gzip
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Find unclassified reads in another SAM')
    parser.add_argument('sam_1', type=str,
                        help='the SAM file for which unclassified reads will be found')
    parser.add_argument('sam_2', type=str,
                        help='the SAM file for which those reads will be outputted')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()

    unclassified_reads = set()
    with get_open_func(args.sam_1)(args.sam_1, 'rt') as sam_1_file:
        for line in sam_1_file:
            parts = line.strip().split('\t')
            read_name, classification = parts[0], parts[1]
            if read_name == 'readID' or classification == 'unclassified':  # include the header as well as unclassified reads
                unclassified_reads.add(read_name)

    with get_open_func(args.sam_2)(args.sam_2, 'rt') as sam_2_file:
        for line in sam_2_file:
            line = line.strip()
            parts = line.split('\t')
            read_name = parts[0]
            if read_name in unclassified_reads:
                print(line)


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
