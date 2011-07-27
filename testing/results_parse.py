## results-parse.py -- take output from scythe run on simulated reads
## and generate statistics.

from optparse import OptionParser
import re
import pdb

parser = OptionParser()
parser.add_option("-r", "--readsfile", dest="reads_file",
                  help="original reads", metavar="FILE")
parser.add_option("-m", "--matchesfile", dest="matches_file",
                  help="matches file output from scythe", metavar="FILE")
parser.add_option("-t", "--trimmedfile", dest="trimmed_file",
                  help="trimmed file output from scythe", metavar="FILE")

(options, args) = parser.parse_args()

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
    
def matchesIter(filename):
    """
    Return a dict with the elements of the matches file block.
    """
    matches_file = open(filename, 'r')
    for line in matches_file:
        p_contam, p_ncontam, adapter = [k.split(': ')[1] for k in line.split('; ')]
        header = matches_file.next().strip()
        read = matches_file.next().strip()
        matches = matches_file.next().strip()
        seq = matches_file.next().strip()
        qual = matches_file.next().strip()
        matches_file.next().strip()
        
        false_pos = re.search("-uncontaminated", header) is not None
        yield {'false_pos':false_pos, 'header':header, 'p_contam':p_contam, 'p_ncontam':p_ncontam}
        

## Gather initial statistics on the simulated reads actual contamination rate
uncontaminated = 0
contaminated = 0
total = 0
for block in FASTQIter(options.reads_file):
    if re.search("-uncontaminated", block['header']) is not None:
        uncontaminated += 1
    if re.search("-contaminated", block['header']) is not None:
        contaminated += 1
    total += 1
print "contaminated:", contaminated
print "uncontaminated:", uncontaminated
print "contamination rate:", contaminated/float(total)
print "total:", total

## Gather statistics on the number of false positives via the matches file
false_positives = 0
total_positives = 0
for block in matchesIter(options.matches_file):
    false_positives += 1 if block['false_pos'] else 0
    total_positives += 1

print "false positives:", false_positives
print "total positives:", total_positives
print "false positive rate:", false_positives/float(total)

print "false negatives:", contaminated-total_positives
print "false negative rate:", (contaminated-total_positives)/float(total)
print "sensitivity:", (total_positives - false_positives)/float(total_positives)
print "specificity:", 1-false_positives/float(total)
