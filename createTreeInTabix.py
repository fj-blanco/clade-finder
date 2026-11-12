# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 13:22:52 2019

@author: hunte
Additional contributions by: fj-blanco
"""

import os
import json
import argparse
from pathlib import Path
import sys
import tempfile

hierarchy = {}
childMap = {}
snps = {}

TABIX_OUTPUT_DIR = "./tabix_output"

def parseTreeJSON(fil):
    hierarchy.clear()
    childMap.clear()
    snps.clear()
    thefile = open(fil)
    root = json.load(thefile)
    thefile.close()
    recurseTreeJson(root)
    return (root["id"])

def replaceAsNecessary(snp):
    return snp.replace("(","_L_PAREN_").replace(")","_R_PAREN_").replace("+","_PLUS_").replace("-","_MINUS_").replace(" ","").replace(".","_DOT_")

def getToIgnore(file):
    ignore = []
    with open(file, "r") as r:
        for line in r.readlines():
            snp = line.replace("\n","")
            if snp != "":
                ignore.append(snp)
    return ignore

def parseSNPsString(snpsString):
    thesnps = set([])
    for snp in snpsString.split(", "):        
        replaced = replaceAsNecessary(snp)
        if replaced != "":
            thesnps.add(replaced)
    return thesnps
            
def recurseTreeJson(node):
    global hierarchy
    global snps
    global childMap
    if "children" in node:
        childMap[node["id"]] = []
        for child in node["children"]:
            if child["id"][-1] != "*":
                childMap[node["id"]].append(child["id"])
                hierarchy[child["id"]] = node["id"]
                snps[child["id"]] = parseSNPsString(child["snps"])
                recurseTreeJson(child)

def getMappingOfSamenameSNPtoUniq():
    uniqsnps = set([])
    for clade in snps:
        for snp in snps[clade]:
            uniqsnps.add(snp)
    samenameSNPToUniqSNP = {}
    for snp in uniqsnps:
        if "/" in snp:
            samenamesnps = snp.split("/")
            for samenamesnp in samenamesnps:
                samenameSNPToUniqSNP[samenamesnp] = snp
    
    return samenameSNPToUniqSNP

def getUniqSNPtoProducts(productsFilePath):
    snpToProducts = {}
    with open(productsFilePath, "r") as r:
        for line in r.readlines():
            splt = line.replace("\n","").split("\t")
            if len(splt) == 2:
                snp = replaceAsNecessary(splt[0])
                snpToProducts[snp] = splt[1]
    r.close()
    uniqSNPtoProducts = {}
    samenameSNPtoUniqSNP = getMappingOfSamenameSNPtoUniq()
    for snp in snpToProducts:
        if snp in samenameSNPtoUniqSNP:
            uniqSNPtoProducts[samenameSNPtoUniqSNP[snp]] = snpToProducts[snp]
        else:
            uniqSNPtoProducts[snp] = snpToProducts[snp]
    return uniqSNPtoProducts

def processPositionMarkers(inputTSVPath, outputFilePath):
    with open(outputFilePath, "w") as w:
        with open(inputTSVPath, "r") as r:
            for line in r.readlines():
                splt = line.strip().split("\t")
                if len(splt) == 3 and splt[0] != "":
                    marker_safe = replaceAsNecessary(splt[1])
                    w.write("\t".join([splt[0], "1", "1", marker_safe, splt[2]]) + "\n")
                else:
                    print("ignored: " + ",".join(splt))

def compressed_path(file_path):
    path = Path(file_path)
    return path.with_suffix('.bgz') if path.suffix == '.tsv' else Path(str(path) + '.bgz')

def create_tabix_file(file_path):
    import pysam
    path = Path(file_path)
    output = compressed_path(path)
    with path.open() as source:
        rows = sorted(source, key=lambda row: (row.split('\t')[0], int(row.split('\t')[1]), row))
    with tempfile.TemporaryDirectory(dir=path.parent) as temporary:
        sorted_path = Path(temporary) / 'sorted.tsv'
        sorted_path.write_text(''.join(rows))
        compressed = Path(temporary) / 'sorted.bgz'
        pysam.tabix_compress(str(sorted_path), str(compressed), force=True)
        pysam.tabix_index(str(compressed), seq_col=0, start_col=1, end_col=2, zerobased=False)
        os.replace(compressed, output)
        os.replace(str(compressed) + '.tbi', str(output) + '.tbi')
    return output

def createTextFile(cladeSNPFilePath, SNPcladeFilePath, uniqSnpToProducts, toIgnore=(),
                   hg19positionMarkersTSV=None, hg38positionMarkersTSV=None,
                   hg19positionMarkersFilePath=None, hg38positionMarkersFilePath=None):
    with open(cladeSNPFilePath, "w") as w:
        for clade in snps:
            for snp in snps[clade]:
                w.write("\t".join([clade, "1", "1", snp, "."]) + "\n")
            if clade in hierarchy:
                w.write("\t".join([clade, "2", "2", hierarchy[clade], "."]) + "\n")
            if clade in childMap:
                for child in childMap[clade]:
                    w.write("\t".join([clade, "3", "3", child, "."]) + "\n")

    with open(SNPcladeFilePath, "w") as w:
        for clade in snps:
            for snp in snps[clade]:
                w.write("\t".join([snp, "1", "1", clade, "."]) + "\n")
                if "/" in snp:
                    for same_name_snp in snp.split("/"):
                        if same_name_snp not in toIgnore:
                            w.write("\t".join([same_name_snp, "2", "2", snp, "."]) + "\n")
        for uniqSNP in uniqSnpToProducts:
            w.write("\t".join([uniqSNP, "3", "3", uniqSnpToProducts[uniqSNP], "."]) + "\n")

    # Handle hg19 position markers if file paths are provided
    if hg19positionMarkersTSV and hg19positionMarkersFilePath:
        processPositionMarkers(hg19positionMarkersTSV, hg19positionMarkersFilePath)

    # Handle hg38 position markers if file paths are provided
    if hg38positionMarkersTSV and hg38positionMarkersFilePath:
        processPositionMarkers(hg38positionMarkersTSV, hg38positionMarkersFilePath)

def main(args):
    toIgnore = getToIgnore(args.toIgnoreFilePath) if args.toIgnoreFilePath else []
    parseTreeJSON(args.treeFile)
    if not snps:
        raise ValueError('The tree has no usable clades.')
    uniqSnpToProducts = getUniqSNPtoProducts(args.productsFilePath) if args.productsFilePath else {}
    output_dir = Path(args.output_dir)
    clade_path = Path(args.cladeSNPFilePath) if args.cladeSNPFilePath else output_dir / 'cladeSNP.tsv'
    snp_path = Path(args.SNPcladeFilePath) if args.SNPcladeFilePath else output_dir / 'SNPclade.tsv'
    outputs = [clade_path, snp_path]
    inputs = [args.treeFile, args.toIgnoreFilePath, args.productsFilePath]
    for build in ('hg19', 'hg38'):
        source = getattr(args, build + 'positionMarkersTSV')
        target = getattr(args, build + 'positionMarkersFilePath')
        if target and not source:
            raise ValueError('--{}-output requires --{}-markers.'.format(build, build))
        if source:
            target = target or str(output_dir / (build + 'positionMarkers.tsv'))
            setattr(args, build + 'positionMarkersFilePath', target)
            inputs.append(source)
            outputs.append(Path(target))
    # Check every destination before opening any output file.
    destinations = list(outputs)
    if args.create_tabix_files:
        import pysam
        for path in outputs:
            compressed = compressed_path(path)
            destinations.extend([compressed, Path(str(compressed) + '.tbi')])
    resolved = [path.resolve() for path in destinations]
    if len(set(resolved)) != len(resolved):
        raise ValueError('Output paths must be distinct.')
    for source in filter(None, inputs):
        if not Path(source).is_file():
            raise ValueError('Missing input file: {}'.format(source))
        if Path(source).resolve() in resolved:
            raise ValueError('An output path would overwrite an input file: {}'.format(source))
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
    createTextFile(str(clade_path), str(snp_path), uniqSnpToProducts, toIgnore,
                   args.hg19positionMarkersTSV, args.hg38positionMarkersTSV,
                   args.hg19positionMarkersFilePath, args.hg38positionMarkersFilePath)
    if args.create_tabix_files:
        for path in outputs:
            print(create_tabix_file(path))

def parse_args(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    # Keep the positional interface used by createTabix.sh.
    if len(argv) == 9 and not any(arg.startswith('--') for arg in argv):
        names = ['--hg19-markers', '--hg38-markers', '--clade-snps', '--snp-clades',
                 '--hg19-output', '--hg38-output', '--products', '--ignore']
        argv = [argv[0]] + [value for pair in zip(names, argv[1:]) for value in pair]
    parser = argparse.ArgumentParser(description='Build a local SNP database from a YFull tree JSON file.')
    parser.add_argument('treeFile', help='YFull tree JSON file.')
    parser.add_argument('--output-dir', default=TABIX_OUTPUT_DIR, help='Output directory (default: tabix_output).')
    parser.add_argument('--clade-snps', '--cladeSNPFilePath', dest='cladeSNPFilePath', default='',
                        help='Override the clade-to-SNP TSV output path.')
    parser.add_argument('--snp-clades', '--SNPcladeFilePath', dest='SNPcladeFilePath', default='',
                        help='Override the SNP-to-clade TSV output path.')
    for build in ('hg19', 'hg38'):
        parser.add_argument('--' + build + '-markers', '--' + build + 'positionMarkersTSV',
                            dest=build + 'positionMarkersTSV', default='',
                            help=build + ' marker TSV: position, SNP name, derived allele.')
        parser.add_argument('--' + build + '-output', '--' + build + 'positionMarkersFilePath',
                            dest=build + 'positionMarkersFilePath', default='',
                            help='Override the ' + build + ' marker output path.')
    parser.add_argument('--products', '--productsFilePath', dest='productsFilePath', default='',
                        help='Optional TSV mapping SNP names to YSEQ product IDs.')
    parser.add_argument('--ignore', '--toIgnoreFilePath', dest='toIgnoreFilePath', default='',
                        help='Optional file of SNP aliases to omit, one per line.')
    parser.add_argument('--index', '--create_tabix_files', dest='create_tabix_files', action='store_true',
                        help='Also create .bgz files and Tabix indexes alongside each TSV.')
    return parser.parse_args(argv)

if __name__ == '__main__':
    try:
        main(parse_args())
    except (OSError, ValueError, KeyError) as error:
        sys.exit('{}: error: {}'.format(Path(sys.argv[0]).name, error))
