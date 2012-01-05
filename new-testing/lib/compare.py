## compare.py -- compare the output of a contaminant trimmer with the
## original simulated file to gauge accuracy.

from optparse import OptionParser
import re
import sys
import pdb

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
    Iterates over a pair of FASTQ Files, returning them in a list of
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

def FASTQStore(reads_file, trimmed_file):
    """
    Iterate over two FASTQ entries
    """
    fastq_files = {'reads': open(reads_file, 'r'), 'trimmed': open(trimmed_file, 'r')}
    store = dict()
    for key in fastq_files.keys():
        tmp_file = fastq_files[key]
        for line in tmp_file:
            header = line.strip()[1:]
            seq = tmp_file.next().strip()
            tmp_file.next()
            qual = tmp_file.next().strip()

            found = store.get(header)
            if not found:
                store[header] = {'reads': dict(), 'trimmed': dict()}
            store[header][key] = {'seq': seq, 'qual':qual}
    return store

contaminated = 0
uncontaminated = 0
trimmed = 0
untrimmed = 0
total = 0
left_out = 0 # the trimmer removed a read entirely from the dataset

total_positives = 0
total_negatives = 0
false_positives = 0
true_positives = 0
true_negatives = 0
false_negatives = 0
contam_diff = dict()


## New method in which pre-trimmed and trimmed reads are matched via
## dict keys. This allows for cases in which some trimmers remove
## entire FASTQ entries.
store = FASTQStore(options.reads_file, options.trimmed_file)

## We need to verify nothing funky is going on with the data.
## First, check that all headers (keys) have reads entries:

if not all([len(data['reads']) == 2 for header, data in store.items()]):
    raise Exception, "Simulated reads do not all have entries - are you using different reads file from trimmed file?"

for header in store:
    total += 1

    block = store[header]

    # Check if this is a contaminated or uncontaminated entry using
    # headers from simulated reads.
    if re.search("-uncontaminated", header) is not None:
        uncontaminated += 1
        is_contaminated = False
    if re.search("-contaminated", header) is not None:
        # note we explicitly check both to check sums are correct
        is_contaminated = True
        contaminated += 1

    ## Check if the trimmer removed this read entirely
    if len(block['trimmed']) == 0:
        left_out += 1
        is_trimmed = True # we consider this trimming the entire read

        offset = len(block['reads']['seq'])
        contam_diff[offset] = contam_diff.get(offset, 0) + 1
    elif len(block['trimmed']) == 2:
        # For simple accuracy
        is_trimmed = len(block['reads']['seq']) > len(block['trimmed']['seq'])

        # More complex (not only trimmed, but trimmed correctly)
        m = re.match(r".*-contaminated-([0-9]+)", header)
        if m is not None and is_trimmed and is_contaminated:
            m = int(m.group(1))
            offset = (len(block['reads']['seq']) - len(block['trimmed']['seq'])) - m
            contam_diff[offset] = contam_diff.get(offset, 0) + 1
    else:
        is_trimmed = False

    trimmed += int(is_trimmed)
    untrimmed += int(not is_trimmed)

    if is_trimmed:
        total_positives += 1
    if is_trimmed and is_contaminated:
        true_positives += 1
    if is_trimmed and not is_contaminated:
        false_positives += 1
    if not is_trimmed and is_contaminated:
        false_negatives += 1
    if not is_trimmed and not is_contaminated:
        true_negatives += 1

# for reads_block, trimmed_block in FASTQPairIter(options.reads_file, options.trimmed_file):
#     total += 1
#     # Check if this is a contaminated or uncontaminated entry using
#     # headers from simulated reads.
#     if re.search("-uncontaminated", reads_block['header']) is not None:
#         uncontaminated += 1
#         is_contaminated = False
#     if re.search("-contaminated", reads_block['header']) is not None:
#         # note we explicitly check both to check sums are correct
#         is_contaminated = True
#         contaminated += 1

#     # For simple accuracy
#     is_trimmed = len(reads_block['seq']) > len(trimmed_block['seq'])

#     # More complex (not only trimmed, but trimmed correctly)
#     m = re.match(r".*-contaminated-([0-9]+)", reads_block['header'])
#     if m is not None and is_trimmed and is_contaminated:
#         m = int(m.group(1))
#         offset = (len(reads_block['seq']) - len(trimmed_block['seq'])) - m
#         contam_diff[offset] = contam_diff.get(offset, 0) + 1

#     trimmed += int(is_trimmed)
#     untrimmed += int(not is_trimmed)

#     if is_trimmed:
#         total_positives += 1
#     if is_trimmed and is_contaminated:
#         true_positives += 1
#     if is_trimmed and not is_contaminated:
#         false_positives += 1
#     if not is_trimmed and is_contaminated:
#         false_negatives += 1
#     if not is_trimmed and not is_contaminated:
#         true_negatives += 1

## check that confusion matrix values add up as they should
assert(contaminated == true_positives + false_negatives)
assert(uncontaminated == false_positives + true_negatives)

# make string of offsets
offset_string = "".join(["%s\t%s\n" % (k, v) for k, v in contam_diff.items()])
sys.stderr.write(offset_string)


print "contaminated\t", contaminated
print "uncontaminated\t", uncontaminated
print "contamination rate\t", contaminated/float(total) if total > 0 else 0
print "total\t", total
print "trimmed\t", trimmed
print "untrimmed\t", untrimmed
print "left_out\t", left_out

true_positives = total_positives - false_positives
total_negatives = true_negatives + false_negatives

if total_positives > 0:
    ## check that confusion matrix values add up as they should
    assert(contaminated == true_positives + false_negatives)
    assert(uncontaminated == false_positives + true_negatives)    

print "false positives\t", false_positives
print "true positives\t", true_positives
print "total positives\t", total_positives

if total_positives > 0:
    print "false positive rate\t", false_positives/float(total)
else:
    print "false positive rate\tNA"
    
print "true negatives\t", true_negatives
print "false negatives\t", false_negatives
print "total negatives\t", total_negatives
print "false negative rate\t", (false_negatives)/float(total)

if true_negatives > 0:
    print "sensitivity\t", true_positives/float(contaminated) if contaminated > 0 else 0
    print "specificity\t", true_negatives/float(false_positives + true_negatives)
else:
    print "sensitivity\tNA"
    print "specificity\tNA"
