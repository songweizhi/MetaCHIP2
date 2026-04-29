import os
import glob
import random
import dendropy
import argparse
from ete3 import Tree
from distutils.spawn import find_executable


def concat_trees(tree_file_1, tree_file_2, concatenated_tree):

    tree1 = Tree(tree_file_1, quoted_node_names=True, format=1)
    tree2 = Tree(tree_file_2, quoted_node_names=True, format=1)
    tree1_with_root = Tree()
    tree2_with_root = Tree()
    tree1_with_root.add_child(tree1)
    tree2_with_root.add_child(tree2)
    concat_tree = Tree()
    concat_tree.add_child(tree1_with_root)
    concat_tree.add_child(tree2_with_root)
    concat_tree.write(outfile=concatenated_tree, quoted_node_names=True, format=1)

rooted_tree_bac         = 'rooted_species_tree_bac.tree'
rooted_tree_ar          = 'rooted_species_tree_ar.tree'
rooted_tree_combined    = 'rooted_species_tree_combined.tree'

concat_trees(rooted_tree_bac, rooted_tree_ar, rooted_tree_combined)
