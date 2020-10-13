TnSeq data processing workflow
================================
Michael Jahn, Johannes Asplund-Samuelsson

## Description

The purpose of the following guide is to provide a simple, step-wise procedure for **Tn-Seq data analysis**. Tn-Seq (transposon sequencing) is the mapping of barcoded transposons to positions in a reference genome. This is the first of two major steps when analyzing barcoded transposon libraries. The second step is the extraction, PCR amplification, and sequencing of the already mapped barcodes from a transposon mutant library (see [Price et al., Nature, 2018](http://www.nature.com/articles/s41586-018-0124-0)). The initial steps of processing next generation sequencing data was directly adapted from [Morgan Price's Feba repository](https://bitbucket.org/berkeleylab/feba/src/master/). The user is referred to the Feba repository or the original Tn-BarSeq publications for further information.

## Prerequistes

- Linux environment with `bash`, `perl` installed
- a DNA alignment software, the default here is `blat`. It can be downloaded [here](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads)
- Perl scripts from M. Price/A. Arkin Lab, Berkeley
- `Fastq` sequencing data as obtained from Illumina runs (see `data/fastq`)
- Reference genome in `.fasta` format
- Model file containing the structure of the read, see `feba/primers/`

## Step 1: Read processing and mapping


All steps follow the description of the `Feba` workflow from Morgan Price. The first step is to run the read mapping script, `MapTnSeq.pl`. The following syntax and options are used with this script.

```
MapTnSeq.pl [ -debug ] [ -limit maxReads ] [ -minQuality $minQuality ]
            [-flanking $flanking]  [ -wobble $wobbleAllowed ]
            [ -minIdentity $minIdentity ] [ -minScore $minScore ] [ -delta $delta ]
            [ -tileSize $tileSize ] [ -stepSize $stepSize ]
            [-tmpdir $tmpdir ] [ -unmapped saveTo ] [ -trunc saveTo ]
            -genome fasta_file -model model_file -first fastq_file > output_file
```

**Important notes (directly taken from Feba repository)**

- The fastq file should have phred+33 ("sanger") encoding of quality scores (as in MiSeq or 2012+ HiSeq). If it is named .gz, it will be gunzipped before reading. If it contains paired-end reads, the second read (name matching " 2:") will be ignored.

- The model file contains 1-2 lines. The first line of the model file describes the beginning of the
read to the junction between the transposon and the genome. For `model_pKMW7` (the associated primer/PCR method), this is:

```
nnnnnnCGCCCTGCAGGGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCGGTTGAGATGTGTATAAGAGACAG
```

where the beginning nnnnnn means 6 random nucleotides and the 20 Ns are the random barcode. The second line of the model file is optional. It describes the vector's sequence after the inverted repeat (that is, after where the junction should be).  This allows `MapTnSeq.pl` to identify intact vector that did not integrate into the target genome. For instance, in `model_pKMW7`, this sequence is:

```
TCGACGGCTTGGTTTCATCAGCCATCCGCTTGCCCTCATCTGTTACGCCGGCGGTAGCCGGCCAGCCTCGCAGAGC
```

- All characters in the model read must be ACGT except for the optional block of n\'s at the front and a block of Ns which represents the barcode.

- This script does not handle sample multiplexing.

- The output file is tab-delimited and contains, for each usable read, the read name, the barcode, which scaffold the insertion lies in, the position of the insertion, the strand that the read matched, a boolean flag for if this mapping location is unique or not, the beginning and end of the hit to the genome in the read after trimming the transposon sequence, the bit score, and the % identity.

- `minQuality` specifies the minimum quality score for each character in the barcode. `flanking` specifies the minimum number of nucleotides on each side that must match exactly.

- `tileSize` and `stepSize` are parameters for BLAT.  If the portion of the read after the junction is short, try a lower stepSize. Also see BLAT's documentation.

- Use `-unmapped` or `-trunc` to save the portion of the reads past the junction. Only reads with barcodes and junctions are saved. `-unmapped` writes only the unmapped reads, in fasta format, with the barcode and the read name as the definition line. `-trunc` writes the remaining part of all of the reads with junctions, in fastq format, and it appends :barcode to the end of the read name (but before the space).

**Example**

- Download and copy `blat` executable to the `feba/bin/` directory so that the perl script finds it directly

- In a terminal, run the following command with customized file paths. The following example works with the data contained in this repository.

- Note: after testing of the different parameters for `MapTnSeq.pl`, none of them made a bigger difference except for reducing `-stepSize` and `-tileSize` from 11 to 7 (+10% reads mapped to genome).

```
cd feba/bin/

perl MapTnSeq.pl -genome ../../ref/GCF_000009285.1_ASM928v2_genomic.fna -model ../primers/model_pKMW7 -first ../../data/fastq/H16_S2_L001_R1_001.fastq.gz > ../../data/mapped/H16_S2_L001_R1_001.tsv
```


## Step 2: Filtering and quality control

Quoted from the Feba repository: `DesignRandomPool.pl` uses the output of `MapTnSeq.pl` to identify barcodes that consistently map to a unique location in the genome. These are the useful barcodes. Ideally, a mutant library has even insertions across the genome; has insertions in most of the protein-coding genes (except the essential ones); has a similar number of reads for insertions in most genes (i.e., no crazy positive selection for loss of a few genes); has insertions on both strands of genes (insertions on only one strand may indicate that the resistance marker's promoter is too weak); has tens of thousands or hundreds of thousands of useful barcodes; and the useful barcodes account for most of the reads.

The syntax for `DesignRandomPool.pl` is:
```
DesignRandomPool.pl -pool pool_file -genes genes_table [ -minN $minN ]
          [ -minFrac $minFrac ] [ -minRatio $minRatio ] [ -maxQBeg $maxQBeg ]
          MapTnSeq_file1 ... MapTnSeq_fileN
```

The script identifies the reliably mapped barcodes and writes a pool file, as well as auxilliary files pool.hit (strains per
gene), pool.unhit (proteins without insertions), pool.surprise.

The MapTnSeq files must be tab-delimited with fields `read, barcode, scaffold, pos, strand, uniq, qBeg, qEnd, score, identity` where qBeg and qEnd are the positions in the read, after the transposon sequence is removed, that match the genome. Gzipped (*.gz) input files are also supported.

The genes table must include the fields `scaffoldId, begin, end, strand, desc` and it should either include only protein-coding genes or it should include the field 'type' with type=1 for protein-coding genes.

**Optional arguments**

- `minN` is the minimum number of "good" reads for a barcode supporting its mapping.  (Good means a unique hit to the genome and qBeg=1.)
- `minFrac` is the minimum fraction of input reads for the barcode that agree with the preferred mapping.
- `minRatio` is the minimum ratio of reads for preferred mapping over the 2nd-most-frequent mapping.

**Example**

- Input files are the mapped reads in `data/mapped/*.tsv` and the gene list
- Note that the perl script calls `../lib/PoolStats.R` at the end which throws errors. Outcommented for now.
- Option `minN` was set to 1 (instead of 10) to assess distribution of all possible Tn insertions
- Output is a tab-separated table with one unique transposn/barcode per row, its frequency and other data

```
perl DesignRandomPool.pl -minN 1 -pool ../../data/pool/H16_barcode_pool.tsv -genes ../../ref/GCF_000009285.1_ASM928v2_genomic_trimmed.tsv ../../data/mapped/H16_S2_L001_R1_001.tsv
```

----------