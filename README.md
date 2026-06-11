# ont-demultiplex

Custom demultiplexer for Oxford Nanopore Technologies (ONT) long-read FASTQ data.

Searches for barcodes at the 5' and 3' ends of each read in both forward and
reverse-complement orientations, using edit distance (via [edlib](https://github.com/Martinsos/edlib)).

This tool is described in a manuscript currently in submission (citation to be added upon publication).
It can be applied to data produced using protocols described in
[Karimian *et al.* (2024) *Science*](https://doi.org/10.1126/science.adh1937).

---

## Scripts

### `demultiplex_ont.py`
The main demultiplexer. Processes one or more FASTQ files and writes per-barcode
gzipped FASTQ output. Reads are assigned to the barcode with the single best
(lowest) edit distance match at either end of the read.

A read is marked **ambiguous** if the two best-matching barcodes are within
1 edit distance of each other, it is then written to a separate `*_ambiguous.fastq.gz`
file.

### `evaluate_error_rates.py`
Tests a range of edit-distance thresholds on a FASTQ file and reports how many
reads would be unambiguously assigned vs. ambiguous at each level. Can be used to evaluate 
`--max_errors` for your data.

---

## Installation

Requires Python 3 and the following libraries:

```bash
pip install biopython edlib
```

---

## Barcode file format

A plain CSV with a label and a sequence column and one barcode per line, no header:

```
SampleA.NB68,GAATCTAAGCAAACACGAAGGTGG
SampleB.NB88,TCTTCTACTACCGATCCGAAGCAG
```


---

## Usage

```bash
python demultiplex_ont.py \
    -i "reads/*.fastq.gz" \
    -b barcodes.csv \
    -o demux_output/ \
    -d 200 \
    -e 4
```

> **Note:** Quote glob patterns (`"*.fastq.gz"`) so the shell passes them
> to the script unexpanded.

| Argument | Description | Default |
|----------|-------------|---------|
| `-i` | Input FASTQ file(s); glob patterns accepted | required |
| `-b` | Barcode CSV file | required |
| `-o` | Output directory | required |
| `-d` | Bases from each end to search | 200 |
| `-e` | Maximum edit distance for a match | 4 |

### Output files

For each input file `reads.fastq.gz` and barcodes `SampleA`, `SampleB`, the script writes:

```
demux_output/
├── reads_SampleA_d200_e4.fastq.gz
├── reads_SampleB_d200_e4.fastq.gz
├── reads_ambiguous_d200_e4.fastq.gz
└── reads_unmatched_d200_e4.fastq.gz
```

The search distance and max errors are encoded in the output filenames to
prevent accidental overwriting when re-running with different parameters.

A per-file summary is printed to stdout, including a breakdown of matches by
end (head/tail) and orientation (forward/reverse-complement):

```
Barcode: SampleA (Total: 199)
  ├─ Found at Head: 0 (fwd) + 122 (rev) = 122
  └─ Found at Tail: 73 (fwd) + 4 (rev) = 77
```

---

## Example

A 300-read test FASTQ with reads from two barcodes and a
two-barcode CSV are provided in `example/`:

```bash
python demultiplex_ont.py \
    -i example/example_reads_mixed.fastq.gz \
    -b example/barcodes_example.csv \
    -o example/output/ \
    -d 200 \
    -e 4
```

Expected output:
```
Total reads processed: 300
Unambiguous matches: 295 (98.33%)
Ambiguous reads:     0 (0.00%)
Unmatched reads:     5 (1.67%)
-----------------------------------
Barcode: SampleA.NB68 (Total: 199)
  ├─ Found at Head: 0 (fwd) + 122 (rev) = 122
  └─ Found at Tail: 73 (fwd) + 4 (rev) = 77
Barcode: SampleB.NB88 (Total: 96)
  ├─ Found at Head: 0 (fwd) + 60 (rev) = 60
  └─ Found at Tail: 36 (fwd) + 0 (rev) = 36
```

---

## Notes on parameter choice

- **`--search_dist 200`** — Searches the first and last 200 bp of each read.
- **`--max_errors 4`** — Works well for 24 bp barcodes and chemistry used in Karimian 2024.
