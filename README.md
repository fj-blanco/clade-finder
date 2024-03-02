# clade-finder

This is a fork of Provyn's [clade-finder](https://github.com/hprovyn/clade-finder). The CladeFinder is running on YSEQ's website ([YSEQ Clade Finder](https://cladefinder.yseq.net/)) and has been proven to be a useful tool for Y-DNA research.

This fork is intended to be used in my batch processing pipelines to process .vcf and .csv with variant calls. I will document the usage of the repository for such cases as I go along. I will remove the php web app. I will focus on Y-DNA as well, at least in principle.

## Installation

Clone the repository and install the Python requirements into your Python environment:

```bash
pip install -r requirements.txt
```

On Ubuntu/Debian-based Linux Distributions:

You need to install `tabix` from the official package repositories by running:

### On Ubuntu/Debian-based Linux Distributions

```bash
sudo apt-get update
sudo apt-get install -y tabix
```

## Creating the tree in Tabix format

Get the [YFull tree](https://github.com/YFullTeam/YTree), for example (you may want to check for the latest version):

```bash
mkdir -p yfull_tree
cd yfull_tree
wget -O ytree.zip https://github.com/YFullTeam/YTree/blob/master/ytree/tree_11.05.0.zip?raw=true
unzip ytree.zip
rm ytree.zip
cd ..
```

Now create the tree in Tabix format:

```bash
python createTreeInTabix.py yfull_tree/tree_11.05.0.json --create_tabix_files
```
