#!/usr/bin/env python3
"""
Counts the number of reads in a read set that are mostly (more than half) N.

Input:  one or more fastq files (can be gzipped)
Output: four tab-delimited columns:
          1) filename
          2) number of reads
          3) number of N reads
          4) percentage of N reads

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

import gzip
import sys


def main():
    n_count, all_count = 0, 0
    for filename in sys.argv[1:]:
        n, a = count_fastq(filename)
        n_count += n
        all_count += a
    n_frac = 100.0 * n_count / all_count
    print('{}\t{}\t{}\t{:.3f}%'.format(filename, all_count, n_count, n_frac))


def get_compression_type(filename):
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for filetype, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = filetype
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def count_fastq(fastq_filename):
    if get_compression_type(fastq_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    reads = []
    n_count, all_count = 0, 0
    with open_func(fastq_filename, 'rt') as fastq:
        for line in fastq:
            sequence = next(fastq).strip().lower()
            next(fastq)
            next(fastq)
            all_count += 1
            if sequence.count('n') / len(sequence) > 0.5:
                n_count += 1
    return n_count, all_count


if __name__ == '__main__':
    main()
