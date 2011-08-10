## results-parse.py -- take output from scythe run on simulated reads
## and generate statistics.

from optparse import OptionParser
import re

parser = OptionParser()
parser.add_option("-r", "--readsfile", dest="reads_file",
                  help="original reads", metavar="FILE")
parser.add_option("-m", "--matchesfile", dest="matches_file",
                  help="matches file output from scythe", metavar="FILE")
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
        base_probs = [float(k) for k in matches_file.next().strip()[1:-1].split(', ')] # remove [, ], & split
        matches_file.next().strip()
        
        false_pos = re.search("-uncontaminated", header) is not None
        yield {'false_pos':false_pos, 'header':header, 'p_contam':p_contam, 'p_ncontam':p_ncontam, 'base_probs':base_probs}
        

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
print "contaminated\t", contaminated
print "uncontaminated\t", uncontaminated
print "contamination rate\t", contaminated/float(total) if total > 0 else 0
print "total\t", total

## Gather statistics on the number of false positives via the matches file
false_positives = 0
total_positives = 0
false_positives_qual = list()
true_positives_qual = list()
for block in matchesIter(options.matches_file):
    is_false_positive = block['false_pos']
    false_positives += 1 if is_false_positive else 0
    base_probs = mean(block['base_probs'])
    if is_false_positive:
        false_positives_qual.append(base_probs)
    else:
        true_positives_qual.append(base_probs)
    total_positives += 1

## Gather statistics on the number of false negatives from clipped file
false_negatives = 0
true_negatives = 0
for block in FASTQIter(options.trimmed_file):
    if re.search(r"-contaminated-[0-9]+$", block['header']) is not None:
        false_negatives += 1
    if re.search(r"-uncontaminated$", block['header']) is not None:
        true_negatives += 1

if total_positives > 0:
    true_positives = total_positives - false_positives
    total_negatives = true_negatives + false_negatives
    assert(contaminated == true_positives + false_negatives)
    assert(uncontaminated == false_positives + true_negatives)    

    print "false positives\t", false_positives
    print "true positives\t", true_positives
    print "total positives\t", total_positives
    print "false positive rate\t", false_positives/float(total)
    print "false positive mean qual\t", mean(false_positives_qual)
    print "true positive mean qual\t", mean(true_positives_qual)

    print "true negatives\t", true_negatives
    print "false negatives\t", false_negatives
    print "total negatives\t", true_negatives
    print "false negative rate\t", (false_negatives)/float(total)
    print "sensitivity\t", true_positives/float(contaminated) if contaminated > 0 else 0
    print "specificity\t", 1-false_positives/float(total)
