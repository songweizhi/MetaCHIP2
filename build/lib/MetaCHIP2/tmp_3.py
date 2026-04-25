import os
import glob
import random
import dendropy
import argparse
from ete3 import Tree
from distutils.spawn import find_executable


def root_tree_at_middle_point(input_tree, tree_file_rooted):

    tree = Tree(input_tree, quoted_node_names=True, format=1)
    midpoint = tree.get_midpoint_outgroup()
    tree.set_outgroup(midpoint)
    tree_str = tree.write(format=1)
    with open(tree_file_rooted, 'w') as f:
        f.write(tree_str)

root_tree_at_middle_point('/Users/songweizhi/Desktop/666/ar53.unrooted.tree', '/Users/songweizhi/Desktop/666/ar53.rooted.tree')
print(os.path.isfile('/Users/songweizhi/Desktop/666/ar53.rooted.tree'))
