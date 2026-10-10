#!/usr/bin/env bash
# ntsm_count.sh - ntsm k-mer counts from an aligned BAM/CRAM, reading only the
# site regions instead of streaming the whole file.
#
#   ntsm_count.sh -s FASTA -b BED [options] sample.bam [more.bam|.cram ...]
#
#   -s FASTA  ntsm site FASTA (human_sites_n10.fa)
#   -b BED    SNP-site BED4 matching the input assembly and the site FASTA
#   -r REF    reference FASTA the reads were aligned to, faidx-indexed; CRAM only
#   -o DIR    output directory (default: current directory)
#   -t N      samtools threads (default 4)
#
# The BED supplies the coordinates; the reference assembly is not checked.
# Contig-name conversion (chr1 vs 1) does not convert genomic coordinates.
#
# Output per input: DIR/{basename}.counts.txt (the ntsmCount table, ready for
# ntsmEval) and DIR/{basename}.counts.log.  Inputs must be coordinate-sorted;
# a missing or stale index is rebuilt.  Needs samtools and ntsmCount on PATH.
set -euo pipefail
usage() {
  awk 'NR==1{next} /^#/{print; next} {exit}' "$0" >&2
  exit 1
}
die() {
  echo "ntsm_count.sh: $*" >&2
  exit 1
}
sample() {
  local b
  b=$(basename "$1")
  b=${b%.bam}
  echo "${b%.cram}"
}

SITES=
REF=
BED=
OUT=.
T=4
while getopts "s:r:b:o:t:h" o; do case $o in
  s) SITES=$OPTARG ;; r) REF=$OPTARG ;; b) BED=$OPTARG ;;
  o) OUT=$OPTARG ;; t) T=$OPTARG ;; *) usage ;; esac done
shift $((OPTIND - 1))
if [ $# -lt 1 ] || [ -z "$SITES" ] || [ -z "$BED" ]; then usage; fi
[ -s "$SITES" ] || die "site FASTA not found: $SITES"
[ -s "$BED" ] || die "site BED not found or empty: $BED"
[ -z "$REF" ] || [ -s "$REF.fai" ] || die "$REF is not faidx-indexed (samtools faidx $REF)"
command -v samtools >/dev/null || die "samtools not on PATH"
command -v ntsmCount >/dev/null || die "ntsmCount not on PATH"
for IN in "$@"; do
  [ -s "$IN" ] || die "input not found: $IN"
  case $IN in *.cram) [ -n "$REF" ] || die "$IN: CRAM input needs -r REF" ;; esac
done
dup=$(for IN in "$@"; do sample "$IN"; done | sort | uniq -d)
[ -z "$dup" ] || die "inputs share an output name: $dup"
mkdir -p "$OUT"
# inside OUT so the final mv is a same-filesystem rename
tmp=$(mktemp -d "$OUT/.ntsm_count.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

nsites=$(cut -f4 "$BED" | sort -u | wc -l)
nbed=$(wc -l <"$BED")

# first non-empty index newer than the data; given to samtools via -X, not left to htslib's lookup
find_index() {
  local i
  for i in "$1.crai" "$1.bai" "${1%.bam}.bai" "$1.csi"; do
    if [ -s "$i" ] && [ "$i" -nt "$1" ]; then
      echo "$i"
      return
    fi
  done
}

for IN in "$@"; do
  base=$(sample "$IN")
  out="$OUT/$base.counts.txt"
  log="$OUT/$base.counts.log"
  refopt=()
  case $IN in *.cram) refopt=(--reference "$REF") ;; esac

  # BED contig names in the style of this file (chr1 vs 1), restricted to its contigs
  samtools view -H ${refopt[@]+"${refopt[@]}"} "$IN" |
    awk '$1=="@SQ"{for(i=2;i<=NF;i++) if($i~/^SN:/) print substr($i,4)}' >"$tmp/ctg.txt"
  awk -v OFS='\t' 'NR==FNR{c[$1]; next}
    { n=$1; sub(/^chr/,"",n); if(!(n in c)) n="chr" n; if(n in c) print n,$2,$3,$4 }' \
    "$tmp/ctg.txt" "$BED" >"$tmp/sites.bed"
  nreg=$(wc -l <"$tmp/sites.bed")
  [ "$nreg" -gt 0 ] || die "$base: no site region lies on a contig of $IN (contigs not named 1-22 or chr1-22?)"
  idx=$(find_index "$IN")
  if [ -z "$idx" ]; then
    case $IN in *.cram) idx=$IN.crai ;; *) idx=$IN.bai ;; esac
    echo "$base: no index newer than the input, writing $idx (input must be coordinate-sorted)" >&2
    samtools index -@ "$T" -o "$idx" "$IN"
  fi

  {
    echo "input: $IN"
    echo "index: $idx"
    echo "sites: $SITES"
    echo "bed: $BED ($nreg regions used of $nbed)"
  } >"$log"
  samtools view -@ "$T" -u -M -L "$tmp/sites.bed" -X ${refopt[@]+"${refopt[@]}"} "$IN" "$idx" 2>>"$log" |
    samtools fastq - 2>>"$log" |
    ntsmCount -s "$SITES" /dev/stdin >"$tmp/counts.txt" 2>>"$log" ||
    die "$base: extraction failed, see $log"
  mv "$tmp/counts.txt" "$out"

  hit=$(awk '!/^#/ && $2+$3>0{n++} END{print n+0}' "$out")
  echo "sites with k-mer hits: $hit of $nsites" >>"$log"
  echo "$base: $hit of $nsites sites have counts -> $out" >&2
  [ "$hit" -ge $((nsites / 10)) ] || echo "WARNING $base: under 10% of sites have counts (low depth, or a site BED from another build?)" >&2
done
