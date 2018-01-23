# Scythe - A Bayesian adapter trimmer (version 0.994 BETA)

Scythe and all supporting documentation
Copyright (c) Vince Buffalo, 2011-2014

Contact: Vince Buffalo <vsbuffaloAAAAA@gmail.com> (with the poly-A tail removed)

If you wish to report a bug, please open an issue on Github
(http://github.com/vsbuffalo/scythe/issues) so that it can be
tracked. You can contact me as well, but please open an issue first.

## About

Scythe uses a Naive Bayesian approach to classify contaminant
substrings in sequence reads. It considers quality information, which
can make it robust in picking out 3'-end adapters, which often include
poor quality bases.

Most next generation sequencing reads have deteriorating quality
towards the 3'-end. It's common for a quality-based trimmer to be
employed before mapping, assemblies, and analysis to remove these poor
quality bases. However, quality-based trimming could remove bases that
are helpful in identifying (and removing) 3'-end adapter
contaminants. Thus, it is recommended you run Scythe *before*
quality-based trimming, as part of a read quality control pipeline.

The Bayesian approach Scythe uses compares two likelihood models: the
probability of seeing the matches in a sequence given contamination,
and not given contamination. Given that the read is contaminated, the
probability of seeing a certain number of matches and mismatches is a
function of the quality of the sequence. Given the read is not
contaminated (and is thus assumed to be random sequence), the
probability of seeing a certain number of matches and mismatches is
chance. The posterior is calculated across both these likelihood
models, and the class (contaminated or not contaminated) with the
maximum posterior probability is the class selected.

## Requirements

Scythe can be compiled using GCC or Clang; compilation during
development used the latter. Scythe relies on Heng Li's kseq.h, which
is bundled with the source.

Scythe requires Zlib, which can be obtained at <http://www.zlib.net/>.

## Building and Installing Scythe

To build Scythe, enter:

    make all

Then, copy or move "scythe" to a directory in your $PATH.

## Usage

Scythe can be run minimally with:

    scythe -a adap.fa -o trimmed_sequences.fasta sequences.fastq

By default, the prior contamination rate is 0.3. This can be changed
(and one is encouraged to do so!) with:

    scythe -a adap.fa -p 0.1 -o trimmed_sequences.fastq sequences.fastq

If you'd like to use standard out, it is recommended you use the
`--quiet` option:

    scythe -a adap.fa --quiet sequences.fastq > trimmed_sequences.fastq

Also, more detailed output about matches can be obtained with:

<<<<<<< HEAD
    scythe -a adap.fa -o trimmed_sequences.fasta -m matches.txt sequences.fastq

By default, the Sanger fastq quality encoding (phred+33; pipeline >= 1.8) is used. 
Illumina (phred+64; pipelines 1.3 - 1.7) or Solexa ("Solexa"+64; pipelines < 1.3) 
qualities can be specified with -q:

    scythe -a adap.fa -q solexa -o trimmed_sequences.fasta sequences.fastq

Lastly, one can specify the minimum match length argument with `-n
<integer>` and the minimum length of sequence (discarded less than or
equal to this parameter) to keep after trimming with `-M <integer>`:

    scythe -a adapter_file.fasta -n 0 -M 10 -o trimmed_sequences.fastq sequences.fastq

The default is 5. If this pre-processing is upstream of assembly on a
very contaminated lane, decreasing this parameter could lead to *very*
liberal trimming, i.e. of only a few bases.

## Notes

<<<<<<< HEAD
Note that the provided adapter sequence files (*_adapters.fa) contain
non-FASTA characters to denote the locations of barcode sequences,
which always appear in TruSeq adapters, and may or may not appear in
forward and/or reverse reads using the original ("Solexa") Illumina
adapter sequences, depending on library preparation. You'll need to
modify the adapter sequence files in order to use them.

An example adapters file (adap.fa) is included for ease of use. It
omits barcodes, so it can be used on all samples of an indexed pool or
set of files (the first ~30 bp are sufficient to identify adapter
contamination). However, Scythe will check reads against all adapter
sequences in the file, so a file with 6 adapter sequences will cause
roughly 6x runtime. Since your samples will never include adapters
from several types of kits, you're encouraged to omit everything but
the adapter(s) that will be found in your sequences.

Scythe only checks for 3'-end contaminants. As of commit `7f49366`,
the algorithm has changd, and Scythe now matches 3' contaminants up to
the 5'-end. Still, do not use Scythe for 5'-trimming as (1) this is a
trivial problem for most Illumina sequences, as the quality is high in
this region, and (2) Scythe does not allow for adapters overlap the
5'-end (though this may come in the future).

For reads with contamination in *any* position, the program TagDust
(<http://genome.gsc.riken.jp/osc/english/dataresource/>) is
recommended. Scythe has the advantages of allowing fuzzier matching
and being base quality-aware, while TagDust has the advantages of very
fast matching (but allowing few mismatches, and not considering
quality) and FDR. Note that TagDust removes contaminated reads
*entirely*, while Scythe trims off contaminating sequence, leaving
valuable uncontaminated read sequence!

A possible pipeline would run FASTQ reads through Scythe, then
TagDust, then a quality-based trimmer, and finally through a read
quality statistics program such as qrqc
(<http://bioconductor.org/packages/devel/bioc/html/qrqc.html>) or FASTqc
(<http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/>).

## FAQ

### Does Scythe work with paired-end data?

Scythe does work with paired-end data. Each file must be run
separately, but Scythe will not remove reads entirely leaving
mismatched pairs.

In some cases, barcodes are ligated to both the 3'-end and 5'-end of
reads. 5'-end removal is trivial since base calling is near-perfect
there, but 3'-end removal can be trickier. Some users have created
Scythe adapter files that contain all possible barcodes concatenated
with possible adapters, so that both can be recognized and
removed. This has worked well and is recommended for cases when 3'-end
quality deteriorates and prevents barcode removal. Newer Illumina
chemistry (TruSeq) has the barcode separated from the fragment, so
that it appears as an entirely separate read that is used to
demultiplex sample reads by the Illumina pipeline.

### Does Scythe work on 5'-end or other contaminants?

No. Embracing the Unix tool philosophy that tools should do one thing
very well, Scythe just removes 3'-end contaminants where there could
be multiple base mismatches due to poor base quality. N-mismatch
algorithms (such as TagDust) don't consider base qualities. Scythe
will allow more mismatches in an alignment if the mismatched bases are
of low quality.

**Scythe only checks as far in as the entire adapter contaminant's
length.** However, some investigation has shown that Illumina
pipelines sometimes produce reads longer than the read length +
adapter length. The extra bases have always been observed to be
A's. Some testing has shown this can be addressed by appending A's to
the adapters in the adapters file. Since Scythe begins by checking for
contamination from the 5'-end of the adapter, this won't affect the
normal adapter contaminant cases.

### What does the numeric output from Scythe mean?

For each adapter in the file, the contaminants removed by position are
returned via standard error. For example:

    Adapter 1 'fake adapter' contamination occurences:
    [10, 2, 4, 5, 6]

indicates that "fake adapter" is 5 bases long (the length of the array
returned), and that there were 10 contaminants found of first base (-n
was set to 0 then), 2 of the first two bases, 4 contaminants of the
first 3 bases, 5 of the first 4 bases, etc.

### Does Scythe work on FASTA files?

No, as these have no quality information.

### How can I report a bug?

See the section below.

### How does Scythe compare to program "x"?

As far as I know, Scythe is the only program that employs a Bayesian
model that allows prior contaminant estimates to be used. This prior
is a more realistic approach than setting a fixed number of mismatches
because we can visually estimate it with the Unix tool `less`.

Scythe also looks at base-level qualities, *not* just a fixed level of
mismatches. A fixed number of mismatches is a bad approach with data
our group (the UC Davis Bioinformatics Core) has seen, as a small bad
quality run can quickly exhaust even a high numbers of fixed
mismatches and lead to higher false negatives.

### How do I cite Scythe?

Just link to this Github repository, until I finish the manuscript
(sorry).

## Reporting Bugs

Scythe is free software and is proved without a warranty. However, I
am proud of this software and I will do my best to provide updates,
bug fixes, and additional documentation as needed. Please report all
bugs and issues to Github's issue tracker
(http://github.com/vsbuffalo/scythe/issues). If you want to email me,
do so in addition to an issue request.

If you have a suggestion or comment on Scythe's methods, you can email
me directly.

## Is there a paper about Scythe?

I am currently writing a paper on Scythe's methods. In my preliminary
testing, Scythe has fewer false positives and false negatives than its
competitors.
