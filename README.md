# Scythe - A very simple adapter trimmer (version 0.93 BETA)

Contact: Vince Buffalo <vsbuffaloAAAAA@gmail.com> (with the poly-A tail removed)

Copyright (c) 2011 The Regents of University of California, Davis Campus.

## License: The MIT License

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

## About

Scythe uses a Naive Bayesian approach to classify contaminant
substrings in sequence reads. It considers quality information, which
can make it robust in picking out 3'-end adapters, which often include
poor quality bases.

The Bayesian approach compares to likelihood models, or the
probability of contaminant (C) given a sequence (S) and quality
(Q). Given that the read is contaminated, the probability of seeing a
certain number of matches and mistmatches is a function of the quality
of the sequence. Given the read is not contaminated (and is thus
random sequence), the probability of seeing a certain number of
matches and mismatches is a function of chance. The posterior is
calculated across both these likelihood models, and the class
(contaminated or not contaminted) with the maximum posterior
probability is the class selected.

## Requirements

Scythe can be compiled using GCC or Clang; compilation during
development used the latter. Scythe relies on Heng Li's kseq.h, which
is bundled with the source.

Scythe requires Zlib, which can be obtained at <http://www.zlib.net/>.

## Building and Installing Scythe

To build Scythe, cd into the code/ directory and enter:

    make

Then, copy or move "scythe" to a directory in your $PATH.

## Usage

*Note: this is outdated (and will be updated soon). Please refer to
 scythe --help for the newest version.*

Scythe can be run minimally with:

    scythe -a adapter_file.fasta sequences.fastq > trimmed_sequences.fasta

By default, the prior contamination rate is 0.05. This can be changed
(and one is encouraged to do so!) with:

    scythe -a adapter_file.fasta -p 0.1 -o trimmed_sequences.fasta sequences.fastq

If you'd like to not use standard out, an output file can be
explicitly specified with:

    scythe -a adapter_file.fasta -o trimmed_sequences.fasta sequences.fastq 

Also, more detailed output about matches can be obtained with:

    scythe -a adapter_file.fasta -o trimmed_sequences.fasta -m matches.txt sequences.fastq

Lastly, for debugging, everything can be printed to standard out with:

    scythe -a adapter_file.fasta -d sequences.fastq
