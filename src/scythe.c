/*
  Scythe - A Simple Bayesian Adapter Contamination Trimmer
*/

#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <zlib.h>

#include "scythe.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)
static const float default_prior = 0.05;

#ifndef PROGRAM_NAME
#define PROGRAM_NAME "scythe"
#endif

#ifndef AUTHORS
#define AUTHORS "Vince Buffalo, UC Davis Bioinformatics Core\nEmail: <vsbuffaloAAAAA@ucdavis.edu> (poly-A tail removed)"
#endif

#ifndef VERSION
#define VERSION 0
#endif

/* Options drawn from GNU's coreutils/src/system.h */
/* These options are defined so as to avoid conflicting with option
   values used by commands */
enum {
  GETOPT_HELP_CHAR = (CHAR_MIN - 2),
  GETOPT_VERSION_CHAR = (CHAR_MIN - 3)
};
#define GETOPT_HELP_OPTION_DECL \
  "help", no_argument, NULL, GETOPT_HELP_CHAR
#define GETOPT_VERSION_OPTION_DECL \
  "version", no_argument, NULL, GETOPT_VERSION_CHAR
#define case_GETOPT_HELP_CHAR                   \
  case GETOPT_HELP_CHAR:                        \
    usage(EXIT_SUCCESS);                        \
    break;
#define case_GETOPT_VERSION_CHAR(Program_name, Version, Authors)        \
  case GETOPT_VERSION_CHAR:                                             \
  fprintf(stdout, "%s version %0.3f\nCopyright (c) 2011 The Regents "   \
          "of University of California, Davis Campus.\n"                \
          "%s is free software and comes with ABSOLUTELY NO WARRANTY.\n" \
          "Distributed under the MIT License.\n\nWritten by %s\n",      \
          Program_name, (double) Version, Program_name, Authors);       \
  exit(EXIT_SUCCESS);                                                   \
  break;
/* end code drawn from system.h */


static struct option long_options[] = {
  /* Options with an argument */
  {"adapter-file", required_argument, 0, 'a'},
  {"prior", required_argument, 0, 'p'},
  {"quality-type", required_argument, 0, 'q'},
  {"matches-file", required_argument, 0, 'm'},
  {"output-file", required_argument, 0, 'c'},
  {"min-match", required_argument, 0, 'n'},
  {"quiet", no_argument, 0, 'Q'},
  {"tag", no_argument, 0, 't'},
  {GETOPT_HELP_OPTION_DECL},
  {GETOPT_VERSION_OPTION_DECL},
  {NULL, 0, NULL, 0}
};

void usage(int status) {
  fputs("\nUsage: scythe -a adapter_file.fasta sequence_file.fastq\n\
Trim 3'-end adapter contaminants off sequence files. If no output file\n\
is specified, scythe will use stdout.\n\
\n\
Options:\n", stdout);
  printf("\
  -p, --prior		prior (default: %0.3f)\n", default_prior);
  fputs("\
  -q, --quality-type	quality type, either illumina, solexa, or sanger (default: illumina)\n\
  -m, --matches-file	matches file (default: no output)\n\
  -o, --output-file	output trimmed sequences file (default: stdout)\n\
  -t, --tag		add a tag to the header indicating Scythe cut a sequence (default: off)\n", stdout);
  fputs("\
  -n, --min-match	smallest contaminant to consider (default: 5)\n\
  --quiet		don't output statistics about trimming to stdout (default: off)\n\
  --help		display this help and exit\n\
  --version		output version information and exit\n", stdout);
  fputs("\n\
  Information on quality schemes:\n\
  phred			PHRED quality scores (e.g. from Roche 454). ASCII with no offset, range: [4, 60].\n\
  sanger		Sanger are PHRED ASCII qualities with an offset of 33, range: [0, 93]. From \n\
			NCBI SRA, or Illumina pipeline 1.8+.\n\
  solexa		Solexa (also very early Illumina - pipeline < 1.3). ASCII offset of\n", stdout);
  fputs("\
	 		64, range: [-5, 62]. Uses a different quality-to-probabilities conversion than other\n\
			schemes.\n\
  illumina		Illumina output from pipeline versions between 1.3 and 1.7. ASCII offset of 64,\n\
			range: [0, 62]\n", stdout);
  exit(status);
}

int main(int argc, char *argv[]) {
  kseq_t *seq;
  int i, l, index, min=5;
  int debug=0, verbose=1;
  int contaminated=0, total=0;
  quality_type qual_type=ILLUMINA;
  match *best_match;
  float *qprobs, prior=default_prior;
  adapter_array *aa;
  gzFile adapter_fp=NULL, output_fp=stdout, matches_fp=NULL, fp;
  int optc;
  char *match_string;
  char tag[14] = "";
  extern char *optarg;

  while (1) {
    int option_index = 0;
    optc = getopt_long(argc, argv, "dtp:a:o:q:m:o:n:", long_options, &option_index);

    if (optc == -1)
       break;
    switch (optc) {
      if (long_options[option_index].flag != 0)
        break;
      case 'a':
        adapter_fp = gzopen(optarg, "r");
        if (!adapter_fp) {
          fprintf(stderr, "Could not open adapter file '%s'.\n", optarg);
          return EXIT_FAILURE;
        }
        break;
      case 'd':
        debug = 1;
        break;
      case 't':
        strcpy(tag, ";;cut_scythe");
        printf("%s\n",tag);
        break;
      case 'Q':
        verbose = 0;
        break;
      case 'o':
        output_fp = fopen(optarg, "w");
        if (!output_fp) {
          fprintf(stderr, "Could not open output file '%s'.\n", optarg);
          return EXIT_FAILURE;
        }
        break;
      case 'm':
        matches_fp = fopen(optarg, "w");
        if (!matches_fp) {
          fprintf(stderr, "Could not open matches file '%s'.\n", optarg);
          return EXIT_FAILURE;
        }
        break;
      case 'n':
        min = atoi(optarg);
        break;
      case 'q':
        if (strcmp(optarg, "illumina") == 0)
          qual_type = ILLUMINA;
        else if (strcmp(optarg, "solexa") == 0)
          qual_type = SOLEXA;
        else if (strcmp(optarg, "sanger") == 0)
          qual_type = SANGER;
        else {
          fprintf(stderr, "Unknown quality type '%s'.\n", optarg);
          usage(EXIT_FAILURE);
        }          
        break;
      case 'p':
        /* This truncation is acceptable... priors need not be doubles */
        prior = (float) atof(optarg);
        if (prior > 1 || prior < 0) {
          fprintf(stderr, "Prior must be between 0 and 1\n");
          usage(EXIT_FAILURE);
        }
        break;
        
      case_GETOPT_HELP_CHAR;
      case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS);
      case '?':
        break;
      default:
        usage(EXIT_FAILURE);
      }
  }
  
  if (debug) {
    matches_fp = stdout;
    output_fp = stdout;
  }
  
  if ((index = optind) == argc) {
    fprintf(stderr, "No FASTQ file specified.\n");
    usage(EXIT_FAILURE);
  }

 
  /* load all adapter sequences into memory */
  if (!adapter_fp) {
    fprintf(stderr, "No adapter file specified.\n");
    usage(EXIT_FAILURE);
  }
  aa = load_adapters(adapter_fp);
  gzclose(adapter_fp);

  fp = gzopen(argv[index], "r");

  if (!fp) {
    fprintf(stderr, "FASTQ file '%s' not found.\n", argv[index]);
    return EXIT_FAILURE;
  }


  seq = kseq_init(fp);


  /* Loop through entire file, find best match against adapters, and
     score these best matches. Write trimmed sequences to file (or
     stdout), and record matches in a match file if specifed. 
  */
  while ((l = kseq_read(seq)) >= 0) {
    if (!seq->qual.s) {
      fputs("Sequence file missing or has malformed quality line.\n", stderr);
      usage(EXIT_FAILURE);
    }
      
    qprobs = qual_to_probs(seq->qual.s, qual_type);
    best_match = find_best_match(aa, seq->seq.s, qprobs, prior, 0.25, min);
    
    if (best_match && best_match->ps->is_contam) {
      /* Best match was classified as contaminated, print trimmed
         results and match entry (if specified). */
      contaminated++;
      fprintf(output_fp, 
              "@%s%s\n%.*s\n+%s%s\n%.*s\n", seq->name.s, tag,
              (int) seq->seq.l-best_match->n, seq->seq.s, seq->name.s, tag, 
              (int) seq->seq.l-best_match->n, seq->qual.s);
      
      /* print match for 5'-end; seq->seq.s is point to the
         approperiate place by taking the sequence length and
         subtracting the match length. */    
      
      if (matches_fp) {
        /* Make a string that indicates the position of the matches with "|"s. */
        match_string = fmt_matches((aa->adapters[best_match->adapter_index]).seq, 
                                   &(seq->seq.s)[seq->seq.l-best_match->n], 
                                   best_match->match, best_match->n);
        
        fprintf(matches_fp, "p(c|s): %f; p(!c|s): %f; adapter: %s\n%s\n%s\n%s\n", 
                best_match->ps->contam, best_match->ps->random,
                aa->adapters[best_match->adapter_index].name,
                seq->name.s, match_string, 
                &(seq->qual.s)[seq->qual.l-best_match->n]);
        fprint_float_array(matches_fp, qual_to_probs(&(seq->qual.s)[seq->qual.l-best_match->n], qual_type), best_match->n);
        fprintf(matches_fp, "\n\n");
        free(match_string);
      }

      /* record occurence in adapter position */ 
      aa->adapters[best_match->adapter_index].occurrences[best_match->n-1]++;

    } else {
      /* No trimming done; print sequence as is. */
      fprintf(output_fp, 
              "@%s\n%s\n+%s\n%s\n", seq->name.s, seq->seq.s, seq->name.s, seq->qual.s);
    }
    
    if (best_match)
      destroy_match(best_match);
    free(qprobs);
    total++;
  }
  
  /* If --quiet not specified, print out a summary at the end. */
  if (verbose) {
    printf("prior: %0.3f\n", prior);
    printf("\nAdapter Trimming Complete\ncontaminated: %d, uncontaminated: %d, total: %d\n", 
           contaminated, total-contaminated, total);
    printf("contamination rate: %f", contaminated/(float) total);
    for (i = 0; i < aa->n; i++) {
      printf("\nAdapter %d '%s' contamination occurences:\n", i+1, aa->adapters[i].name);
      print_uint_array(aa->adapters[i].occurrences, aa->adapters[i].length);
      printf("\n");
    }
  }

  kseq_destroy(seq);
  destroy_adapters(aa, MAX_ADAPTERS);
  gzclose(fp);
  return 0;
}
