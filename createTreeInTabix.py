# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 13:22:52 2019

@author: hunte
"""

import os
import json
import argparse
import subprocess

hierarchy = {}
childMap = {}
snps = {}

TABIX_OUTPUT_DIR = "./tabix_output"

def parseTreeJSON(fil):
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

def create_tabix_file(file_path):
    # Determine the base name for the file without the .tsv extension
    base_name = os.path.basename(file_path).replace('.tsv', '')
    output_path = os.path.join(TABIX_OUTPUT_DIR, base_name)

    # Sort, compress, and create the bgz file
    sort_command = f"sort {file_path} -k1,1 -k2n | bgzip > {output_path}.bgz"
    subprocess.run(sort_command, shell=True, check=True)

    # Index the compressed file
    tabix_command = f"tabix -s 1 -b 2 -e 3 {output_path}.bgz"
    subprocess.run(tabix_command, shell=True, check=True)

def createTextFile(cladeSNPFilePath, SNPcladeFilePath, uniqSnpToProducts, toIgnore=[], 
                   hg19positionMarkersTSV=None, hg38positionMarkersTSV=None, 
                   hg19positionMarkersFilePath=None, hg38positionMarkersFilePath=None):
    # Handle Clade to SNP and SNP to Clade mappings
    # We create the cladeSNP file in the TABIX_OUTPUT_DIR folder if it doesn't exist

    with open(cladeSNPFilePath, "w") as w:
        for clade in snps:
            for snp in snps[clade]:
                w.write("\t".join([clade, "1", "1", snp, "."]) + "\n")
            if clade in hierarchy:
                w.write("\t".join([clade, "2", "2", hierarchy[clade], "."]) + "\n")
            if clade in childMap:
                for child in childMap[clade]:
                    w.write("\t".join([clade, "3", "3", child, "."]) + "\n")

    # We create the SNPclade file in the TABIX_OUTPUT_DIR folder if it doesn't exist
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
    treeFile = args.treeFile
    parseTreeJSON(treeFile)

    # Handle optional arguments
    if args.productsFilePath:
        uniqSnpToProducts = getUniqSNPtoProducts(args.productsFilePath)
    else:
        uniqSnpToProducts = {}
    if not args.cladeSNPFilePath:
        cladeSNPFilePath = f"{TABIX_OUTPUT_DIR}/cladeSNP.tsv"
    if not args.SNPcladeFilePath:
        SNPcladeFilePath = f"{TABIX_OUTPUT_DIR}/SNPclade.tsv"

    # Create the cladeSNP and SNPclade files
    createTextFile(cladeSNPFilePath, SNPcladeFilePath, uniqSnpToProducts, toIgnore)
    # Create the tabix files if the flag is set
    if args.create_tabix_files:
        print("Creating tabix files...")
        create_tabix_file(cladeSNPFilePath)
        create_tabix_file(SNPcladeFilePath)

def parse_args():
    parser = argparse.ArgumentParser(description='Process a phylogenetic tree and generate SNP and clade mappings.')
    parser.add_argument('treeFile', type=str, help='Path to the YFull tree JSON file.')
    parser.add_argument('--cladeSNPFilePath', type=str, default='', help='Path to the output file for clade to SNP mappings. Optional')
    parser.add_argument('--SNPcladeFilePath', type=str, default='', help='Path to the output file for SNP to clade mappings. Optional')
    parser.add_argument('--hg19positionMarkersTSV', type=str, default='', help='Path to the TSV file with hg19 position markers. Optional.')
    parser.add_argument('--hg38positionMarkersTSV', type=str, default='', help='Path to the TSV file with hg38 position markers. Optional.')
    parser.add_argument('--hg19positionMarkersFilePath', type=str, default='', help='Output path for the processed hg19 position markers. Optional.')
    parser.add_argument('--hg38positionMarkersFilePath', type=str, default='', help='Output path for the processed hg38 position markers. Optional.')
    parser.add_argument('--productsFilePath', type=str, default='', help='Path to the products file mapping unique SNPs to YSEQ products. Optional.')
    parser.add_argument('--toIgnoreFilePath', type=str, default='', help='Path to the file listing SNPs to ignore. Optional.')
    parser.add_argument('--create_tabix_files', action='store_true', help='Create tabix files for cladeSNP and SNPclade files.')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    # if cladeSNPFilePath or SNPcladeFilePath are not provided, we create them in the TABIX_OUTPUT_DIR folder
    if not args.cladeSNPFilePath or not args.SNPcladeFilePath:
        # Ensure the TABIX_OUTPUT_DIR exists
        if not os.path.exists(TABIX_OUTPUT_DIR):
            os.makedirs(TABIX_OUTPUT_DIR)
    main(args)