# ntsm-bam

Run [ntsm](https://github.com/JustinChu/ntsm) sample-swap detection on aligned BAM/CRAM files, reading only the reads over its 96,287 human SNP sites instead of the whole file.

ntsm decides whether two sequencing runs come from the same individual by counting allele k-mers at a fixed SNP panel. `ntsmCount` normally consumes raw FASTQ. When the data are already aligned, nearly all reads are irrelevant to the panel. `ntsm_count.sh` uses the BAM/CRAM index to fetch only the reads overlapping the sites, streams them into `ntsmCount`, and writes a counts file for `ntsmEval`. Ready-made site BEDs for **hg19** (also b37 and hs37d5) and **hg38** are included.

## Requirements

- bash, awk and coreutils
- [samtools](https://www.htslib.org/) (tested with 1.21)
- ntsm (`ntsmCount`, `ntsmEval`) on `PATH`: `conda install bioconda::ntsm`
- ntsm's human site FASTA, `human_sites_n10.fa` (see Install)
- Coordinate-sorted BAM or CRAM input. CRAM also needs its reference FASTA, faidx-indexed.

## Install

```sh
git clone https://github.com/odinokov/ntsm-bam.git
cd ntsm-bam
curl -LO https://raw.githubusercontent.com/JustinChu/ntsm/663f9a548afa989f04ece08be22ae81a9cc1be6f/data/human_sites_n10.fa
md5sum human_sites_n10.fa    # 2b33daef838a468eda959b1f3be507de
```

The bundled BEDs were built from this exact site FASTA (ntsm commit `663f9a5`). A different site set needs its own BED.

## Quick start

```sh
./ntsm_count.sh -s human_sites_n10.fa -b ntsm_sites.hg38.bed -o counts sample1.bam sample2.bam
ntsmEval counts/sample1.counts.txt counts/sample2.counts.txt > summary.tsv
```

`summary.tsv` has one row per sample pair; `same` = 1 means ntsm calls them the same individual. For the faster PCA-based mode (`ntsmEval -a`) and the meaning of each column, see the [ntsm README](https://github.com/JustinChu/ntsm#readme).

## Usage

```text
ntsm_count.sh -s FASTA -b BED [-r REF] [-o DIR] [-t N] in.bam|in.cram ...
```

| Option | Meaning |
|---|---|
| `-s` | ntsm site FASTA (`human_sites_n10.fa`). Required. |
| `-b` | SNP-site BED4 matching the input reference assembly and the site FASTA. Required. |
| `-r` | Reference the reads were aligned to, faidx-indexed. Required for CRAM. |
| `-o` | Output directory (default `.`). |
| `-t` | samtools threads (default 4). |

The BED supplies the genomic coordinates; the script does not infer or verify the reference assembly. Use `ntsm_sites.hg38.bed` for hg38 alignments, `ntsm_sites.hg19.bed` for hg19/b37/hs37d5 alignments, or a BED for your assembly and site panel. Every input in one call must use the coordinate system of the supplied BED. Contig-name conversion (`chr1` ↔ `1`) does not convert coordinates between assemblies.

All inputs are checked before any is processed: each must exist, a CRAM needs `-r`, and basenames without `.bam`/`.cram` must be distinct.

For each input the script:

1. Converts the BED to the input's contig names (`chr1` or `1`) and drops contigs the header lacks. Exits with an error if no regions remain.
2. Indexes the input if it has no index newer than the file (writes `<input>.bai`/`.crai` next to it).
3. Runs `samtools view -u -M -F 0x900 -L sites.bed | samtools fastq | ntsmCount`.
4. Writes `DIR/<basename>.counts.txt` (only on success, so a failed run leaves no partial table) and `DIR/<basename>.counts.log` (regions used, reads extracted, sites with counts). It warns if under 10% of sites have counts.

### Caveats

- `errorRate` from `ntsmEval` is not meaningful on extracted reads. `score`, `cov` and genotypes are unaffected.
- Only reads aligned over the sites are counted. Unmapped reads, or reads carrying site k-mers that aligned elsewhere, are missed; a full-FASTQ `ntsmCount` run would see them.
- `samtools view -F 0x900` excludes secondary and supplementary alignments before piping to `samtools fastq`, which also excludes them by default. Duplicates are kept.
- `sample.bam` and `sample.cram` share an output basename, so passing both in one call is rejected.

## Site BEDs

Supply a tab-delimited BED4 with columns `contig`, `start`, `end`, and `site ID`. Coordinates are 0-based, half-open: a SNP at 1-based position 100 is `chr1\t99\t100\trsID` (with tab separators). Use one 1-bp interval per SNP, matching the panel in the site FASTA.

The bundled `ntsm_sites.<build>.bed` files use the rsID as the site ID and are naturally sorted, one line per site. Both files cover autosomes only, because ntsm removed the chrX sites from `human_sites_n10.fa` in commit `ab6ecaa`.

## Files

| File | Purpose |
|---|---|
| `ntsm_count.sh` | BAM/CRAM → ntsm counts |
| `ntsm_sites.hg19.bed`, `ntsm_sites.hg38.bed` | site BEDs |

## License

MIT, see [LICENSE](LICENSE). ntsm and its site panel are © Justin Chu, also under the MIT license; please cite [ntsm](https://github.com/JustinChu/ntsm) when you use this tool.
