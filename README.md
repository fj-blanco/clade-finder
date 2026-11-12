# clade-finder

Predict a Y-DNA clade from positive and negative SNP calls using a YFull tree.
The [YSEQ Clade Finder](https://lab.yseq.net/cladefinder/) provides a web interface.
The commands below run locally.

## Local installation

Use a virtual environment. The CLI is tested with Python 3.12.
`pytabix` builds a C extension, so a compiler and Python development headers are required.
On Debian or Ubuntu, install them with:

```bash
sudo apt-get install python3-venv python3-dev build-essential zlib1g-dev
```

Then, from the repository directory:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

The requirements cover the SNP-list CLI and database builder. The older VCF script
also needs `PyVCF3`. Web deployments must install dependencies in the Python
environment used by their server process.

## Build a database

Download and unzip a tree JSON file from [YFull's tree repository](https://github.com/YFullTeam/YTree).
Keep the original file and its version with your analysis so you can reproduce the result.

```bash
python createTreeInTabix.py /path/to/tree.json --output-dir tabix_output --index
```

This creates `cladeSNP.tsv`, `SNPclade.tsv`, their `.bgz` files and `.tbi` indexes.
Indexing uses pysam; separate `bgzip` and `tabix` executables are not needed for this command.
Existing generated files at these paths are replaced when rebuilding.
Without `--index`, the command writes only the TSV files.

Optional inputs:

- `--hg19-markers` and `--hg38-markers`: tab-separated position, SNP name and derived allele.
  Their output files are also compressed and indexed when `--index` is given.
- `--products`: tab-separated SNP name and YSEQ product ID.
- `--ignore`: SNP aliases to omit, one per line. This omits alias lookup entries;
  it does not remove the underlying mutation from the tree.

Run `python createTreeInTabix.py --help` for custom output paths. The old long option
names and nine positional arguments used by `createTabix.sh` are still accepted.

## Predict a clade

Create a text file containing SNP names followed by `+` or `-`, separated by commas
or newlines. Spaces around calls are allowed:

```text
Z6065+, L21-, DF13-, YP879+, M81-
```

```bash
python clade_finder.py calls.txt --database tabix_output
```

The default database directory is `tabix_output`. For an existing database with
different filenames, use `--clade-snps /path/to/cladeSNPs.bgz` and
`--snp-clades /path/to/SNPclades.bgz`. Each file needs its adjacent `.tbi` index.

Use `-` to read stdin and `--all` to return every candidate:

```bash
printf 'Z6065+, YP879+\n' | python clade_finder.py - --all
```

Successful commands write one JSON object, or an array with `--all`, to stdout and
exit with status 0. Input and database errors go to stderr and exit with status 2.
Diagnostics from the classifier also go to stderr, so stdout can be piped into
another program.

The score is the original classifier's ranking score, not a probability. Candidates
are sorted by descending score, then clade name for ties. A tied best score produces
a warning; use `--all` to inspect the alternatives. SNPs absent from the chosen tree
are ignored by the classifier. An input with no usable positive evidence fails.
The CLI accepts signed SNP lists only, not VCF files or consumer genotype CSV files.

## Check the installation

A small synthetic tree is included for testing. It is not a reference for real samples.

```bash
python createTreeInTabix.py test/fixtures/tree.json --output-dir /tmp/clade-example --index
python clade_finder.py test/sample_for_clade_finder.txt --database /tmp/clade-example
```

Expected output:

```json
{"clade": "J-YP879", "score": 2.0, "nextPrediction": {"clade": "J-Z6065", "score": 1.0}}
```

Run the integration tests with:

```bash
python -m unittest discover -s test -v
```
