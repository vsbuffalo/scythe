#!/usr/bin/perl -w

# AUTHOR: Joseph Fass <joseph.fass@gmail.com>
# LAST REVISED: August 2012
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu

# profileTruSeqIDs.pl is the proprietary property of The Regents of
# the University of California (“The Regents.”) Copyright 2007-12 The
# Regents of the University of California, Davis campus. All Rights
# Reserved. Redistribution and use in source and binary forms, with
# or without modification, are permitted by nonprofit, research
# institutions for research use only, provided that the following
# conditions are met: Redistributions of source code must retain the
# above copyright notice, this list of conditions and the following
# disclaimer. Redistributions in binary form must reproduce the above
# copyright notice, this list of conditions and the following
# disclaimer in the documentation and/or other materials provided with
# the distribution. The name of The Regents may not be used to
# endorse or promote products derived from this software without
# specific prior written permission. The end-user understands that
# the program was developed for research purposes and is advised not
# to rely exclusively on the program for any reason. THE SOFTWARE
# PROVIDED IS ON AN "AS IS" BASIS, AND THE REGENTS HAVE NO OBLIGATION
# TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
# MODIFICATIONS. THE REGENTS SPECIFICALLY DISCLAIM ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE REGENTS BE LIABLE TO ANY
# PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, EXEMPLARY OR
# CONSEQUENTIAL DAMAGES, INCLUDING BUT NOT LIMITED TO PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES, LOSS OF USE, DATA OR PROFITS, OR
# BUSINESS INTERRUPTION, HOWEVER CAUSED AND UNDER ANY THEORY OF
# LIABILITY WHETHER IN CONTRACT, STRICT LIABILITY OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE AND ITS DOCUMENTATION, EVEN IF ADVISED OF THE POSSIBILITY
# OF SUCH DAMAGE. If you do not agree to these terms, do not download
# or use the software. This license may be modified only in a writing
# signed by authorized signatory of both parties.

# searches for TruSeq adapters in fastq, extracts and counts the 6 bp
# "barcodes," or multiplex identifiers

while (<>) {
	$seq = <>;
	<>;
	<>;
	if ($seq =~ /AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC([ACGT]{6})/) {
	    $count{$1}++;
	}
}
if (%count) {
    foreach $ID (reverse sort {$count{$a}<=>$count{$b}} keys %count) {
	print "$ID\t$count{$ID}\n";
    }
} else {
    print "\nNo ID's found! Reads too short? Try visual inspection, with a shorter search string, like \"AGATCGGAAGAG\".\n\n";
}
