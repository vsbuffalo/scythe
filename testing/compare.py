## compare.py -- compare the output of a contaminant trimmer with the
## original simulated file to gauge accuracy.

from optparse import OptionParser
import re

parser = OptionParser()
parser.add_option("-r", "--readsfile", dest="reads_file",
                  help="original reads", metavar="FILE")
parser.add_option("-t", "--trimmedfile", dest="trimmed_file",
                  help="trimmed file output from scythe", metavar="FILE")

(options, args) = parser.parse_args()

def mean(ll):
    if len(ll) > 0:
        return sum(ll)/float(len(ll))
    return None

def FASTQIter(filename):
    """
    Return a dict with the elements of a FASTQ block.
    """
    fastq_file = open(filename, 'r')
    for line in fastq_file:
        header = line.strip()[1:]
        seq = fastq_file.next().strip()
        fastq_file.next()
        qual = fastq_file.next().strip()
        yield {'header':header, 'seq':seq, 'qual':qual}

def FASTQPairIter(filename1, filename2):
    """
    Iterates over a pair of FASTQ Files, returning them in tuples of
    dicts.
    """
    fastq_files = (open(filename1, 'r'), open(filename2, 'r'))
    complete = False
    while not complete:
        tmp = list()

        for fastq_file in fastq_files:
            line = fastq_file.next()

            if line is None:
                # check whether other file is done; if not error.
                if fastq_files.next().next() is not None:
                    raise Exception, "cannot continue: different lengthed files!"
                else:
                    # both are done; set breakout condition
                    compelete = True
                    break
                    
            header = line.strip()[1:]
            seq = fastq_file.next().strip()
            fastq_file.next()
            qual = fastq_file.next().strip()

            tmp.append({'seq': seq, 'header':header, 'qual':qual})
        yield tmp

contaminated = 0
uncontaminated = 0
trimmed = 0
untrimmed = 0
total = 0

false_positives = 0
true_positives = 0
true_negatives = 0
false_negatives = 0

for reads_block, trimmed_block in FASTQPairIter(options.reads_file, options.trimmed_file):
    total += 1
    # Check if this is a contaminated or uncontaminated entry using
    # headers from simulated reads.
    if re.search("-uncontaminated", reads_block['header']) is not None:
        uncontaminated += 1
        is_contaminated = False
    if re.search("-contaminated", reads_block['header']) is not None:
        # note we explicitly check both to check sums are correct
        is_contaminated = True
        contaminated += 1

    is_trimmed = len(reads_block['seq']) > len(trimmed_block['seq'])
    trimmed += int(is_trimmed)
    untrimmed += int(not is_trimmed)

    if is_trimmed and is_contaminated:
        true_positives += 1
    if is_trimmed and not is_contaminated:
        false_positives += 1
    if not is_trimmed and is_contaminated:
        false_negatives += 1
    if not is_trimmed and not is_contaminated:
        true_negatives += 1

## check that confusion matrix values add up as they should
assert(contaminated == true_positives + false_negatives)
assert(uncontaminated == false_positives + true_negatives)    
