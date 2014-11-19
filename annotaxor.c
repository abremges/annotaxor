/* The MIT License

   Copyright (c) 2014 Andreas Bremges <andreas@cebitec.uni-bielefeld.de>

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
*/

#include <stdint.h>   // uint64_t
#include <stdio.h>    // printf
#include <stdlib.h>   // calloc
#include <string.h>   // strtok
#include <sys/time.h> // gettimeofday (obsolescent)
#include <time.h>     // time
#include <unistd.h>   // getopt
#include <zlib.h>     // gzip


#include "kseq.h"
KSEQ_INIT(gzFile, gzread) // TODO Switch back
//KSTREAM_INIT(gzFile, gzread, 65536)

#include "khash.h"
KHASH_MAP_INIT_INT(TAX, uint32_t)
khash_t(TAX) *h;

KHASH_MAP_INIT_INT64(LCA, uint32_t)
khash_t(LCA) *hLca;

// typedef struct {
// 	uint32_t rank:4, taxid:28;
// } lca_t; // TODO ?

//static int kmerSize = 8;
//static int alphSize = 15;

static double start;
static int kmerSize = 8;

// Pack each amino acid into 4-bits by using a compressed alphabet
// Murphy et al. (2000), 15 letter alphabet, [0..14] = valid, 15 = N/A
// http://bio.math-inf.uni-greifswald.de/viscose/html/alphabets.html
unsigned char murphy15[128] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//  _, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O,
//  P, Q, R, S, T, U, V, W, X, Y, Z, _, _, _, _, _,
    0, 3, 0, 2,11,10, 8, 4,15, 1, 0,14, 1, 1,12, 0,
    7,13,14, 5, 6, 0, 1, 9, 0, 8, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 2,11,10, 8, 4,15, 1, 0,14, 1, 1,12, 0,
    7,13,14, 5, 6, 0, 1, 9, 0, 8, 0, 0, 0, 0, 0, 0
}; // TODO Go down to 15^8 memory, now 16^8, i..e change murphy15 code (?)


// 2bit-encoded triplets, A=00, C=01, G=10, T=11 (U=11)
unsigned char codonTable[64] = {
//   K, N, K, N, T, T, T, T, R, S, R, S, I, I, I, I,
    14,12,14,12, 6, 6, 6, 6,14, 5,14, 5, 1, 1, 1, 1,
//   Q, H, Q, H, P, P, P, P, R, R, R, R, L, L, L, L,
    13,15,13,15, 7, 7, 7, 7,14,14,14,14, 1, 1, 1, 1,
//   E, D, E, D, A, A, A, A, G, G, G, G, V, V, V, V,
    10,11,10,11, 3, 3, 3, 3, 4, 4, 4, 4, 1, 1, 1, 1,
//   X, Y, X, Y, S, S, S, S, X, C, W, C, L, F, L, F
     0, 8, 0, 8, 5, 5, 5, 5, 0, 2, 9, 2, 1, 8, 1, 8
}; // TODO double-check

// A = 00, C = 01, G = 10, T = 11, i.e. XOR with 3 -> compl.
unsigned char seq_fwd_table[128] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

//////////
// Taken from: https://github.com/lh3/seqtk/blob/master/seqtk.c
//////////
char comp_tab[] = {
    0, 1,   2,      3,      4, 5,   6,      7,      8, 9, 10,       11,     12, 13, 14,     15,
    16, 17, 18,     19,     20, 21, 22,     23,     24, 25, 26,     27,     28, 29, 30,     31,
    32, 33, 34,     35,     36, 37, 38,     39,     40, 41, 42,     43,     44, 45, 46,     47,
    48, 49, 50,     51,     52, 53, 54,     55,     56, 57, 58,     59,     60, 61, 62,     63,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,     92, 93, 94,     95,
    64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};
void seq_revcomp(kseq_t *seq) {
    int c0, c1;
    for (int i = 0; i < seq->seq.l>>1; ++i) { // reverse complement sequence
        c0 = comp_tab[(int)seq->seq.s[i]];
        c1 = comp_tab[(int)seq->seq.s[seq->seq.l - 1 - i]];
        seq->seq.s[i] = c1;
        seq->seq.s[seq->seq.l - 1 - i] = c0;
    }
    if (seq->seq.l & 1) // complement the remaining base
        seq->seq.s[seq->seq.l>>1] = comp_tab[(int)seq->seq.s[seq->seq.l>>1]];
    if (seq->qual.l) {
        for (int i = 0; i < seq->seq.l>>1; ++i) // reverse quality
            c0 = seq->qual.s[i], seq->qual.s[i] = seq->qual.s[seq->qual.l - 1 - i], seq->qual.s[seq->qual.l - 1 - i] = c0;
    }
}

double realtime() {
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

void addTaxon (const uint32_t key, const uint32_t val) {
	int ret;
	khiter_t k = kh_put(TAX, h, key, &ret); // TODO Check if present?
	if (ret) {
		kh_value(h, k) = val;
	} else {
		kh_del(TAX, h, k);
	}
}

void parseTaxonomy(const char *taxFile) {
	int dret;
	gzFile fp = strcmp(taxFile, "-") ? gzopen(taxFile, "r") : gzdopen(fileno(stdin), "r");
	kstream_t *ks = ks_init(fp);
	kstring_t str = {0,0,0};
	while (ks_getuntil2(ks, KS_SEP_LINE, &str, &dret, 0) >= 0) {
		char *key = strtok(str.s, "\t|");
		char *val = strtok(NULL, "\t|");
		addTaxon(atoi(key), atoi(val)); // TODO Sanity check, atoi returns zero upon failure
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str.s);
}

unsigned char as_to_int[128] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//  _, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O,
//  P, Q, R, S, T, U, V, W, X, Y, Z, _, _, _, _, _,
    0, 1, 0, 2, 3, 4, 5, 6, 7, 8, 0, 9,10,11,12, 0,
   13,14,15,16,17, 0,18,19, 0,20, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 2, 3, 4, 5, 6, 7, 8, 0, 9,10,11,12, 0,
   13,14,15,16,17, 0,18,19, 0,20, 0, 0, 0, 0, 0, 0
};

uint64_t getKmerInt(const char *sKmer) {
	uint64_t iKmer = 0;
	for (int i = 0; i < kmerSize; ++i) { // TODO variable k-mer length?
		int c = as_to_int[(int) sKmer[i]]; // TODO check string length!
		iKmer = ((iKmer << 5) | c);
	}
	return iKmer;
}

void parseKmerToLca(const char *lcaFile) {
	int dret;
	gzFile fp = strcmp(lcaFile, "-") ? gzopen(lcaFile, "r") : gzdopen(fileno(stdin), "r");
	kstream_t *ks = ks_init(fp);
	kstring_t str = {0,0,0};
	unsigned long n = 0;
	// TODO Declare some version and format in the first line of DB file, and check for correct one
	while (ks_getuntil2(ks, KS_SEP_LINE, &str, &dret, 0) >= 0) {
		char *kmer = strtok(str.s, "\t");
		char *taxId = strtok(NULL, "\t"); // TODO
		if (kmer && taxId) {
			int ret;
			khiter_t k = kh_put(LCA, hLca, getKmerInt(kmer), &ret);
			if (ret) {
				kh_value(hLca, k) = atoi(taxId);
			} else {
				kh_del(LCA, hLca, k);
			}
		} else {
			// TODO Something went terribly wrong...
			fprintf(stderr, "[%f] FUCK\n", (realtime()-start));
		}
		if ((++n % 1000000) == 0) fprintf(stderr, "[%f] Processed %lu kmers\n", (realtime()-start), n);
	}
	if ((n % 1000000) == 0) fprintf(stderr, "[%f] Processed %lu kmers\n", (realtime()-start), n);
	ks_destroy(ks);
	gzclose(fp);
	free(str.s);
}

KHASH_SET_INIT_INT(SPECIES)
khash_t(SPECIES) *khs;

void prepareTaxHashes() {
		int dret;
		gzFile fp = gzopen("valid_krona_taxa.txt", "r");
		kstream_t *ks = ks_init(fp);
		kstring_t str = {0,0,0};
		while (ks_getuntil2(ks, KS_SEP_LINE, &str, &dret, 0) >= 0) {
			int ret;
			kh_put(SPECIES, khs, atoi(str.s), &ret);
		}
		ks_destroy(ks);
		gzclose(fp);
		free(str.s);
}

uint32_t calcLca(uint32_t foo, uint32_t bar) {
    khash_t(TAX) *hFoo = kh_init(TAX);
	khiter_t iFoo;
	khiter_t k = kh_get(TAX, h, foo);
	while(k != kh_end(h) && kh_value(h, k) > 1) {
        int ret;
        iFoo = kh_put(TAX, hFoo, foo, &ret);
        if (ret) {
            kh_value(hFoo, iFoo) = 1; // TODO Switch to set instead of map?
        } else {
            kh_del(TAX, hFoo, iFoo);
        }
        foo = kh_value(h, k);
		k = kh_get(TAX, h, foo);
	}
    k = kh_get(TAX, h, bar);
	while(k != kh_end(h) && kh_value(h, k) > 1) {
        if ((iFoo = kh_get(TAX, hFoo, bar)) != kh_end(hFoo)) {
            kh_destroy(TAX, hFoo);
            return bar;
        }
		bar = kh_value(h, k);
		k = kh_get(TAX, h, bar);
	}
    kh_destroy(TAX, hFoo);
	return 1;
}

uint32_t getValidTaxId(uint32_t bar) {
    khiter_t k = kh_get(TAX, h, bar);
	while(k != kh_end(h) && kh_value(h, k) > 1) {
        if (kh_get(SPECIES, khs, bar) != kh_end(khs)) {
            return bar;
        }
		bar = kh_value(h, k);
		k = kh_get(TAX, h, bar);
	}
	return 1;
}

uint32_t calcDepth(int key) {
	key = getValidTaxId(key);
	int depth = 1;
	khiter_t k = kh_get(TAX, h, key);
	while(k != kh_end(h) && kh_value(h, k) != 1) {
		++depth;
		key = getValidTaxId(kh_value(h, k));
		k = kh_get(TAX, h, key);
	}
	return (1ULL<<(depth-1))<<(depth-1); // depth
}

uint32_t kmerTax(const uint64_t kmer) {
	khiter_t k = kh_get(LCA, hLca, kmer);
	if (k != kh_end(h)) {
		return kh_value(hLca, k);
	} else {
		return 0;
	}
}

void processSingleRead(const kseq_t *seq) {
	int index = 0;
	uint64_t forward = 0;
// 	fputs(seq->name.s, stdout);
// 	fputc('\t', stdout);
	khash_t(TAX) *th = kh_init(TAX); // TODO rename.
	khiter_t iter;
	for (int i = 0; i < seq->seq.l; ++i) {
		++index;
		int c = murphy15[(int) seq->seq.s[i]];
		uint64_t mask = (1ULL<<(kmerSize*5))-1;
		forward = ((forward << 5) | c) & mask;
		if (index >= kmerSize) { // TODO
			uint32_t taxon = kmerTax(forward);
			if (taxon > 0) {
				iter = kh_get(TAX, th, taxon);
				if (iter == kh_end(th)) {
					int ret;
					iter = kh_put(TAX, th, taxon, &ret);
					kh_value(th, iter) = calcDepth(taxon); // 1;
				} else {
					kh_value(th, iter) += calcDepth(taxon); // ++
				}
// 				fprintf(stdout, "%i;", taxon);
			}
		}
	}
// 	fputc('\t', stdout);
	int maxscore = 0;
	uint32_t maxtaxon = 1;
	for (iter = kh_begin(th); iter != kh_end(th); ++iter) {
		if (kh_exist(th, iter)) {
			uint32_t candidate = kh_key(th, iter);
			uint32_t key = candidate;
			int count = kh_val(th, iter);
				khiter_t k = kh_get(TAX, h, key);
				while(k != kh_end(h) && kh_value(h, k) > 1) {
					key = kh_value(h, k);
					k = kh_get(TAX, h, key);
					khiter_t iter2 = kh_get(TAX, th, key);
					if (iter2 != kh_end(th) && kh_exist(th, iter2)) {
						count += kh_val(th, iter2);
					}
				}
				if (count > maxscore) {
					maxscore = count;
					maxtaxon = candidate;
				} else if (count == maxscore) {
					// deal with co-optimal solutions: don't touch maxscore
					// set maxtaxon to lca of current maxtaxon and candidate
					maxtaxon = calcLca(maxtaxon, candidate);
				}
		}
	}
	fprintf(stdout, "%s\t%i\t%i/%lu\n", seq->name.s, getValidTaxId(maxtaxon), maxscore, seq->seq.l-7);
	kh_destroy(TAX, th);
}

void processNucleotideRead(kseq_t *seq) {
	int gmaxscore = 0;
	uint32_t gmaxtaxon = 1;
	for (int r = 0; r < 2; ++r) {
		for (int n = 0; n < 3; ++n) {
			khash_t(TAX) *th = kh_init(TAX);
			khiter_t iter;
			int index = 0;
			uint64_t forward = 0;
			for (int i = n; i < seq->seq.l-2; i+=3) {
				++index;
				int base, codon = 0;
				for (int p = 0; p < 3; ++p) {
					base = seq_fwd_table[(int) seq->seq.s[i+p]];
					if (base < 4) {
						codon = ((codon << 2) | base);
					} else {
						index = 0;
					}
				}
				int c = codonTable[codon];
				uint64_t mask = (1ULL<<(kmerSize*5))-1;
				forward = ((forward << 4) | c) & mask;
				if (index >= kmerSize) {
					uint32_t taxon = kmerTax(forward);
					if (taxon > 0) {
						iter = kh_get(TAX, th, taxon);
						if (iter == kh_end(th)) {
							int ret;
							iter = kh_put(TAX, th, taxon, &ret);
							kh_value(th, iter) = 1;
						} else {
							kh_value(th, iter)++;
						}
					}
				}
			}
			int maxscore = 0;
			uint32_t maxtaxon = 1;
			for (iter = kh_begin(th); iter != kh_end(th); ++iter) {
				if (kh_exist(th, iter)) {
					uint32_t candidate = kh_key(th, iter);
					
					uint32_t key = candidate;
					int count = kh_val(th, iter);
					
					khiter_t k = kh_get(TAX, h, key);
					while(k != kh_end(h) && kh_value(h, k) > 1) {
						key = kh_value(h, k);
						k = kh_get(TAX, h, key);
						
						khiter_t iter2 = kh_get(TAX, th, key);
						if (iter2 != kh_end(th) && kh_exist(th, iter2)) {
							count += kh_val(th, iter2);
						}
					}
					if (count > maxscore) {
						maxscore = count;
						maxtaxon = candidate;
					} else if (count == maxscore) {
						// deal with co-optimal solutions: don't touch maxscore
						// set maxtaxon to lca of current maxtaxon and candidate
						maxtaxon = calcLca(maxtaxon, candidate);
					}
				}
			}
			if (maxscore > gmaxscore) {
				gmaxscore = maxscore;
				gmaxtaxon = maxtaxon;
			} else if (maxscore == gmaxscore) {
				// deal with co-optimal solutions: don't touch gmaxscore
				// set gmaxtaxon to lca of current gmaxtaxon and maxtaxon
				gmaxtaxon = calcLca(gmaxtaxon, maxtaxon);
			}
			kh_destroy(TAX, th);
		}
		seq_revcomp(seq);
	}
	fprintf(stdout, "%s\t%i\t%i\n", seq->name.s, gmaxtaxon, gmaxscore);
}


void processReadFile(const char *file) {
	gzFile fp = strcmp(file, "-") ? gzopen(file, "r") : gzdopen(fileno(stdin), "r");
	kseq_t *seq = kseq_init(fp);
	int n = 0;
	while (kseq_read(seq) >= 0) {
		if ((++n % 100000) == 0) fprintf(stderr, "[%f] Processed %i protein sequences\n", (realtime()-start), n);
		processSingleRead(seq);
	}
	if ((n % 100000) != 0) fprintf(stderr, "[%f] Processed %i protein sequences\n", (realtime()-start), n);
	kseq_destroy(seq);
	gzclose(fp);
}

void processReadFile2(const char *file) {
	gzFile fp = strcmp(file, "-") ? gzopen(file, "r") : gzdopen(fileno(stdin), "r");
	kseq_t *seq = kseq_init(fp);
	int n = 0;
	while (kseq_read(seq) >= 0) {
		if ((++n % 100000) == 0) fprintf(stderr, "[%f] Processed %i nucleotide sequences\n", (realtime()-start), n);
		processNucleotideRead(seq);
	}
	if ((n % 100000) != 0) fprintf(stderr, "[%f] Processed %i nucleotide sequences\n", (realtime()-start), n);
	kseq_destroy(seq);
	gzclose(fp);
}

static int usage() {
	fprintf(stderr, "anotaxor -t nodes.dmp -l lcafile -p prot.faa -x read.fna\n");
	return 42;
}

int main(int argc, char *argv[]) {
	start = realtime();
	char *faaFile = 0, *fnaFile = 0, *lcaFile = 0, *taxFile = 0;
	int c;
    while((c = getopt(argc, argv, "k:p:x:l:t:")) != -1) {
        switch (c) {
			case 'k':
				kmerSize = atoi(optarg);
				break;
			case 'p':
				faaFile = optarg;
				break;
			case 'x':
				fnaFile = optarg;
				break;
			case 'l':
				lcaFile = optarg;
				break;
			case 't':
				taxFile = optarg;
				break;
            default:
				return usage();
        }
	}
	if (!((faaFile || fnaFile) && lcaFile && taxFile)) return usage();

	h = kh_init(TAX);
	hLca = kh_init(LCA);

    khs = kh_init(SPECIES);
    prepareTaxHashes();

	if (taxFile) {
		parseTaxonomy(taxFile);
	}

	fprintf(stderr, "[%f] Filling kmer to taxon map from %s\n", (realtime()-start), lcaFile);
	parseKmerToLca(lcaFile);

    if (faaFile != 0) {
        fprintf(stderr, "[%f] Processing nucleotide sequences in %s\n", (realtime()-start), faaFile);
        processReadFile(faaFile);
    } else {
        fprintf(stderr, "[%f] Processing protein sequences in %s\n", (realtime()-start), fnaFile);
        processReadFile2(fnaFile);
    }

	fprintf(stderr, "[%f] Done, thank you!\n", (realtime()-start));

	kh_destroy(TAX, h);
    kh_destroy(LCA, hLca);
    kh_destroy(SPECIES, khs);

	return 0;
}
