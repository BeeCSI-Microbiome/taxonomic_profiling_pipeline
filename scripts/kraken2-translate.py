#!/usr/bin/env python
"""
# =============================================================================
GLOBALS & IMPORTS
# =============================================================================
"""
import os
import re
from collections import deque

# ranks that krakefaction needs for lineage
TAXA_RANKS = ['d','p','c','o','f','g','s', 's1']

"""
# =============================================================================
# =============================================================================
"""

"""
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Class - TaxonNode
------
PURPOSE
-------
This class represents a node of a taxonomic level, with links to subtaxa and 
supertaxon.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""
class TaxonNode:

    def __init__(self, clade_count, taxa_count, rank, taxid, name, supertaxon=None):
        self.clade_count = int(clade_count) # num of db entries of this + subtaxa
        self.taxa_count = int(taxa_count)   # num of db entries of this exact taxa
        self.rank = rank.lower()
        self.taxid = taxid
        self.name = name.strip().replace(' ', '_')
        self.supertaxon = supertaxon        # link to supertaxon node
        self.subtaxa = []                   # list of links to subtaxa nodes
        self.subtaxa_sum = int(taxa_count)
        self.lineage = ''

    def create_lineage(self):
        """Returns appropriate lineage string, based on ancestors"""
        # first case (domain)
        if self.rank == 'd':
            return 'd__' + self.name
        elif self.rank in TAXA_RANKS:
            return self.supertaxon.lineage + '|{}__{}'.format(self.rank, self.name)
        elif self.supertaxon:
            return self.supertaxon.lineage
        else:
            return ''

    def add_child(self, child):
        """Add a subtaxon child node to this node
           Add the new node's taxa_count to the subtaxa_sum of all supertaxa"""
        # add to parent's subtaxa list
        self.subtaxa.append(child)
        # add parent as child's supertaxa
        child.supertaxon = self

        if child.taxa_count != 0:
            # add tax_count to parent's sum
            self.subtaxa_sum += child.taxa_count

            curr_parent = self
            # add to any further ancestor's sums
            while curr_parent.supertaxon is not None:
                curr_parent = curr_parent.supertaxon
                curr_parent.subtaxa_sum += child.taxa_count

        # now that child has parent, create lineage
        child.lineage = child.create_lineage()

        return

    def has_full_subtaxa(self):
        """Return boolean whether all subtaxa have been added to this node"""
        return self.subtaxa_sum == self.clade_count

    def __str__(self):
        s = 'Taxon name: {} - Rank {}\n'.format(self.name, self.rank)
        s += '\tLineage: {}\n'.format(self.lineage)
        s += '\tTaxa ID: {}\n'.format(self.taxid)
        s += '\tClade count: {}\n'.format(self.clade_count)
        s += '\tTaxa count: {}\n'.format(self.taxa_count)
        s += '\tSubtaxa sum: {}\n'.format(self.subtaxa_sum)
        if self.supertaxon is not None:
            s += '\tSupertaxa: {}\n'.format(self.supertaxon.name)
        else:
            s += '\tSupertaxa: None\n'
        s += '\tSubtaxa: {}\n'.format(' '.join([subtaxon.name for subtaxon in self.subtaxa]))
        return s

# def get_superlineage(node):
#     """Recursively finds the most recent lineage string from ancestors"""
#     # if supertaxon's lineage is not an empty string, return it
#     if node.supertaxon.lineage:
#         return node.supertaxon.lineage
#     # Else, recurse
#     else:
#         return get_superlineage(node.supertaxon)

"""
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""

"""
# =============================================================================
MAIN
# =============================================================================
"""
def main():
    """Verify inputs, then run script"""
    # get inspection filename and filtered kraken output file paths from snakemake
    db_inspection, infile_paths = get_inputs()
    run(db_inspection, infile_paths)
    
"""
# =============================================================================
# =============================================================================
"""

def get_inputs():
    """Get the inputs from snakemake.
       Database inspection file, and paths to filtered kraken output files"""
    # retrieve the name of the inspection file from snakemake
    db_inspection = snakemake.input[0]

    # retrieve the paths of filtered files from snakemake
    infile_paths = snakemake.input[1:]

    # verifiy files
    check_file_existence(db_inspection, infile_paths)

    return(db_inspection, infile_paths)


def check_file_existence(db_ins, in_paths):
    """Validate file existence, raising error upon failure"""
    # check database inspection file
    if not os.path.isfile(db_ins):
        raise RuntimeError(
            "ERROR: could not open database inspection file {}\n".format(db_ins))
    
    # check each kraken input file
    for f in in_paths:
        if not os.path.isfile(f):
            raise RuntimeError(
                "ERROR: could not open input file {}\n".format(f))

    return
        
"""
# =============================================================================
"""
def run(db_inspection, infile_paths):
    """Main logical control of the script occurs within"""
    taxon_dict = get_taxa_dict(db_inspection)

    for infile in infile_paths:
        f = open(infile, 'r')
        reads = f.readlines()
        f.close()
        create_outfile(reads, taxon_dict, infile)
"""
# =============================================================================
"""    

def get_taxa_dict(db_ins):
    """Produce a taxonomy dictionary with taxids as keys and lineage strings as values"""
    # get database inspection file as list of lines
    with open(db_ins, 'r') as f:
        inspection_lines = f.readlines()

    # create the taxon tree
    root = create_tree(inspection_lines)
    
    # create the taxid/lineage dictionary
    taxa_dict = {}
    build_taxa_dict(taxa_dict, root)

    return taxa_dict


def create_tree(inspection_lines):
    """Create the taxonomy tree and return the root node"""
    # make deque for leftpopping of inspection lines
    inspection_lines = deque(inspection_lines)
    
    # create root node
    root_line = inspection_lines.popleft().split('\t')
    root_node = TaxonNode(*root_line[1:])

    # stack for tracking nodes with more subtaxa to add
    stack = []
    stack.append(root_node)

    # while there are still lines to make nodes from
    while inspection_lines:
        # get the node at the top of the stack
        curr_node = stack[-1]
        new_node = TaxonNode(*inspection_lines.popleft().split('\t')[1:])

        # if current node has full subtaxa list, pop stack until otherwise
        if curr_node.has_full_subtaxa():
            curr_node = find_next_parent(stack)

        # add the new node to the parent node
        curr_node.add_child(new_node)
        # put new node on top of stack if it still needs subtaxa
        if not new_node.has_full_subtaxa():
            stack.append(new_node)
        
    
    return root_node
    

def find_next_parent(stack):
    """Recursively pops stack until a node that still needs children is found"""
    node = stack[-1]
    if not node.has_full_subtaxa():
        return stack[-1]
    else:
        stack.pop()
        return find_next_parent(stack)


def build_taxa_dict(t_dict, node):
    """Recursively add taxid/lineage from tree to taxa dictionary"""
    # add current node if it has lineage
    if node.lineage:
        t_dict[str(node.taxid)] = node.lineage
    else:
        # if it has no lineage, it is root or subroot
        t_dict[str(node.taxid)] = 'root'
    for sub in node.subtaxa:
        build_taxa_dict(t_dict, sub)


def create_outfile(reads, t_dict, filename):
    """Creates the formatted output table file"""
    with open(filename.replace('filtered', 'translated'), 'w') as f:
        for line in reads:
            columns = line.split('\t')
            # data to be first entry in output file
            datum = columns[1]
            # find the taxid
            tax_search = re.search('(?<=taxid )(\d*)', columns[2])
            taxid = tax_search.group(1)
            if taxid in t_dict.keys():
                #print(datum + '\t' + t_dict[taxid])
                if line == reads[-1]:
                    f.write(datum + '\t' + t_dict[taxid])
                else:
                    f.write(datum + '\t' + t_dict[taxid] + '\n')
            #else:
               # l.debug('Unable to translate: {}'.format(columns[2]))


def print_tree(root):
    """Prints out tree top down"""
    print(root)
    for sub in root.subtaxa:
        print_tree(sub)
        
    return

if __name__ == '__main__':
    main()
