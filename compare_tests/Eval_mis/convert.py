#!/usr/bin/env python3
import argparse
import fileinput
import re
import sys
'''
Create a molecule extents file in BED format from contig reports output after running QUAST

@author: Yaman Malkoc
'''



class Molecule:
    "A long-reads molecule representation"

    def __init__(
            self, rname, source, my_type, start,
            end, score, orientation, phase, attributes):
        self.rname = rname
        self.source = source
        self.my_type = my_type
        self.start = start
        self.end = end
        self.score = score
        self.orientation = orientation
        self.phase = phase 
        self.attributes = attributes

    def print_bed(self, file):
        "Print this molecule to a BED file"
        print(
            self.rname, self.source, self.my_type, self.start, self.end,
            self.score, self.orientation, self.phase, self.attributes,
            sep="\t", file=file)

class ConnectedMolecule:
    def __init__(
            self, chr1, start1, end1,
            chr2, start2, end2, unique_id):
        self.chr1 = chr1
        self.start1 = start1
        self.end1 = end1
        self.chr2 = chr2
        self.start2 = start2
        self.end2 = end2
        self.unique_id = unique_id

    def print_bedpe(self, file):
        "Prints this misassembly in a BEDPE format"
        print(self.chr1, self.start1, self.end1,
                self.chr2, self.start2, self.end2,
                self.unique_id, sep="\t", file=file)

def parse_arguments():
    "Parsing the arguments given to the command"
    parser = argparse.ArgumentParser(
            description="Convert a mis_contigs.info representation of misassemblies from a QUAST run to a BED format. ")
    parser.add_argument( 
            metavar="mis_contigs", dest="input_misassemblies_file", 
            help="Input mis_contigs.info file to be converted")
    parser.add_argument("-o", "--bed",
            metavar="FILE", dest="bed",
            help="Output file in BED format that is to be generated after running this command. ")
    parser.add_argument("-b", "--bedpe", metavar="FILE", dest="bedpe",
            help="When set, output file in BEDPE format to show the connection between the misassembled regions")
    
    return parser.parse_args()

def convert_file(input_file, bed, bedpe):
    "Splitting columns of the paf file and filtering out the unnecessary data(for bed file) and rearrange columns."
    with open(input_file, 'r') as infofile:
        bedfile = open(bed, "w+") if bed is not None else sys.stdout
        "Set the tag to add GFF3-style attributes (more information when hovered over the misassemblies)"
        print("##gff-version 3.2.1", file=bedfile)
        bedpefile = open(bedpe, "w+") if bedpe is not None else None
        lines = infofile.readlines()
        current_contig = lines[0]
        first_piece_inverted = False
        second_piece_inverted = False
        misassembly_pattern = re.compile("Extensive misassembly \((.*)\).*between ([0-9]*) ([0-9]*) and ([0-9]*) ([0-9]*)")
        line_count = 1
        for line in lines:
            if line is not None:
                match = misassembly_pattern.match(line)
                
                if(match is None):
                    contig_name = line.split("\n")
                    current_contig = contig_name[0]
                else:
                    misassembly_type = match.group(1)
                    #Assign coordinate values and swap them if they are in an ascending order
                    mis_contig_1a = match.group(2)
                    mis_contig_1b = match.group(3)
                    if (int(mis_contig_1a) > int(mis_contig_1b)):
                        mis_contig_1a, mis_contig_1b = mis_contig_1b, mis_contig_1a
                        first_piece_inverted = True
                    mis_contig_2a = match.group(4)
                    mis_contig_2b = match.group(5)
                    if (int(mis_contig_2a) > int(mis_contig_2b)):
                        mis_contig_2a, mis_contig_2b = mis_contig_2b, mis_contig_2a
                        second_piece_inverted = True
                    if(bedpefile is not None): 
                        newConnectedMolecule = ConnectedMolecule(current_contig, str(mis_contig_1a), str(mis_contig_1b),
                                                                 current_contig, str(mis_contig_2a), str(mis_contig_2b), str(line_count))
                        newConnectedMolecule.print_bedpe(bedpefile)
                        line_count = line_count + 1
                    newMolecule = Molecule(current_contig, ".", str(misassembly_type), str(mis_contig_1a), str(mis_contig_1b), 
                            ".", "+" if first_piece_inverted is False else "-", ".", 
                            "ID="+str(line_count))
                    newMolecule.print_bed(bedfile)
                    line_count = line_count + 1
                    newMolecule = Molecule(current_contig, ".", str(misassembly_type), str(mis_contig_2a), str(mis_contig_2b), 
                            ".", "+" if second_piece_inverted is False else "-", ".", 
                            "ID="+str(line_count) + ";Parent=" + str(line_count - 1))
                    newMolecule.print_bed(bedfile)
                    line_count = line_count + 1
                    first_piece_inverted = False
                    second_piece_inverted = False

        if bedfile is not sys.stdout:
            bedfile.close()

if __name__ == '__main__':
    arguments = parse_arguments()
    input_paf = "/dev/stdin" if arguments.input_misassemblies_file == "-" else arguments.input_misassemblies_file
    convert_file(input_paf, arguments.bed, arguments.bedpe)
