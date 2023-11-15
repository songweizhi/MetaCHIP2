import random
import dendropy
from ete3 import Tree


def get_smallest_outgroup(tree_object):

    min_outgroup_leaf_num = 99999
    for each_root_child in tree_object.children:
        leaf_list = each_root_child.get_leaf_names()
        if len(leaf_list) < min_outgroup_leaf_num:
            min_outgroup_leaf_num = len(leaf_list)

    out_group_leaf_list = []
    for each_root_child in tree_object.children:
        leaf_list = each_root_child.get_leaf_names()
        if len(leaf_list) == min_outgroup_leaf_num:
            out_group_leaf_list = leaf_list

    return out_group_leaf_list


def sep_taxon_str(taxon_string):

    taxon_string_split = taxon_string.strip().split(';')
    taxon_p = taxon_string_split[1]
    taxon_c = taxon_string_split[2]
    taxon_o = taxon_string_split[3]
    taxon_f = taxon_string_split[4]
    taxon_g = taxon_string_split[5]

    return taxon_p, taxon_c, taxon_o, taxon_f, taxon_g


def subset_and_rename_tree(tree_file_in, to_keep_leaf_list, rename_dict):

    input_tree = Tree(tree_file_in, quoted_node_names=True, format=1)

    # subset tree
    subset_tree = input_tree.copy()
    subset_tree.prune(to_keep_leaf_list, preserve_branch_length=True)

    # rename leaf
    for each_leaf in subset_tree:
        leaf_name_new = rename_dict.get(each_leaf.name, each_leaf.name)
        each_leaf.name = leaf_name_new

    return subset_tree


def root_with_outgroup(input_tree, out_group_list, tree_file_rooted):

    """
    Reroot the tree using the given outgroup.
    modified based on: https://github.com/Ecogenomics/GTDBTk/blob/master/gtdbtk/reroot_tree.py

    input_tree:  File containing Newick tree to rerooted.
    output_tree: Name of file for rerooted tree.
    outgroup:    Labels of taxa in outgroup.
    """

    tree = dendropy.Tree.get_from_path(input_tree, schema='newick', rooting='force-rooted', preserve_underscores=True)

    outgroup_in_tree = set()
    ingroup_leaves = set()
    for n in tree.leaf_node_iter():
        if n.taxon.label in out_group_list:
            outgroup_in_tree.add(n.taxon)
        else:
            ingroup_leaves.add(n)

    # Since finding the MRCA is a rooted tree operation, the tree is first rerooted on an ingroup taxa. This
    # ensures the MRCA of the outgroup can be identified so long as the outgroup is monophyletic. If the
    # outgroup is polyphyletic trying to root on it is ill-defined. To try and pick a "good" root for
    # polyphyletic outgroups, random ingroup taxa are selected until two of them give the same size
    # lineage. This will, likely, be the smallest bipartition possible for the given outgroup though
    # this is not guaranteed.

    mrca = tree.mrca(taxa=outgroup_in_tree)
    mrca_leaves = len(mrca.leaf_nodes())
    while True:
        rnd_ingroup = random.sample(list(ingroup_leaves), 1)[0]
        tree.reroot_at_edge(rnd_ingroup.edge, length1=0.5 * rnd_ingroup.edge_length, length2=0.5 * rnd_ingroup.edge_length)
        mrca = tree.mrca(taxa=outgroup_in_tree)
        if len(mrca.leaf_nodes()) == mrca_leaves:
            break

        mrca_leaves = len(mrca.leaf_nodes())

    if mrca.edge_length is not None:
        tree.reroot_at_edge(mrca.edge, length1=0.5 * mrca.edge_length, length2=0.5 * mrca.edge_length)
        tree.write_to_path(tree_file_rooted, schema='newick', suppress_rooting=True, unquoted_underscores=True)


def root_tree_by_gtdb(gtdb_ref_tree, gtdb_gnm_metadata, user_gnm_taxon, user_gnm_tree, user_tree_rooted, gnm_domain):

    tree = Tree(gtdb_ref_tree, quoted_node_names=True, format=1)
    bac120_ref_tree_gnm_list = tree.get_leaf_names()
    bac120_ref_tree_gnm_set = {i for i in bac120_ref_tree_gnm_list}

    # read in user_gnm_taxon
    user_gnm_taxon_dict_p = dict()
    user_gnm_taxon_dict_c = dict()
    user_gnm_taxon_dict_o = dict()
    user_gnm_taxon_dict_f = dict()
    user_gnm_taxon_dict_g = dict()
    for each_gnm in open(user_gnm_taxon):
        if not each_gnm.startswith('user_genome\t'):
            each_gnm_split = each_gnm.strip().split('\t')
            gnm_id = each_gnm_split[0]
            gnm_taxon = each_gnm_split[1]

            count_current_gnm = False
            if gnm_domain == 'bac':
                if 'd__Bacteria' in gnm_taxon:
                    count_current_gnm = True
            elif gnm_domain == 'ar':
                if 'd__Archaea' in gnm_taxon:
                    count_current_gnm = True

            if count_current_gnm is True:

                gnm_p, gnm_c, gnm_o, gnm_f, gnm_g = sep_taxon_str(gnm_taxon)

                if gnm_p not in user_gnm_taxon_dict_p:
                    user_gnm_taxon_dict_p[gnm_p] = set()
                if gnm_c not in user_gnm_taxon_dict_c:
                    user_gnm_taxon_dict_c[gnm_c] = set()
                if gnm_o not in user_gnm_taxon_dict_o:
                    user_gnm_taxon_dict_o[gnm_o] = set()
                if gnm_f not in user_gnm_taxon_dict_f:
                    user_gnm_taxon_dict_f[gnm_f] = set()
                if gnm_g not in user_gnm_taxon_dict_g:
                    user_gnm_taxon_dict_g[gnm_g] = set()

                user_gnm_taxon_dict_p[gnm_p].add(gnm_id)
                user_gnm_taxon_dict_c[gnm_c].add(gnm_id)
                user_gnm_taxon_dict_o[gnm_o].add(gnm_id)
                user_gnm_taxon_dict_f[gnm_f].add(gnm_id)
                user_gnm_taxon_dict_g[gnm_g].add(gnm_id)

    # determine rooting rank, start from phylum
    rooting_rank = ''
    rooting_rank_taxon_dict = dict()
    if len(user_gnm_taxon_dict_p) > 1:
        rooting_rank = 'p'
        rooting_rank_taxon_dict = user_gnm_taxon_dict_p
    elif len(user_gnm_taxon_dict_c) > 1:
        rooting_rank = 'c'
        rooting_rank_taxon_dict = user_gnm_taxon_dict_c
    elif len(user_gnm_taxon_dict_o) > 1:
        rooting_rank = 'o'
        rooting_rank_taxon_dict = user_gnm_taxon_dict_o
    elif len(user_gnm_taxon_dict_f) > 1:
        rooting_rank = 'f'
        rooting_rank_taxon_dict = user_gnm_taxon_dict_f
    elif len(user_gnm_taxon_dict_g) > 1:
        rooting_rank = 'g'
        rooting_rank_taxon_dict = user_gnm_taxon_dict_g

    if rooting_rank == '':
        print('All user genomes are from the same genus, program exited!')
        exit()

    col_index = {}
    canditate_gnms_rooting_rank = dict()
    counted_taxons_rooting_rank = set()
    for each_ref in open(gtdb_gnm_metadata):
            each_ref_split = each_ref.strip().split('\t')
            if each_ref.startswith('accession	ambiguous_bases'):
                col_index = {key: i for i, key in enumerate(each_ref_split)}
            else:
                ref_accession = each_ref_split[0]
                gtdb_taxonomy = each_ref_split[col_index['gtdb_taxonomy']]
                if ref_accession in bac120_ref_tree_gnm_set:
                    gnm_p, gnm_c, gnm_o, gnm_f, gnm_g = sep_taxon_str(gtdb_taxonomy)

                    gnm_rooting_rank = ''
                    if rooting_rank == 'p':
                        gnm_rooting_rank = gnm_p
                    elif rooting_rank == 'c':
                        gnm_rooting_rank = gnm_c
                    elif rooting_rank == 'o':
                        gnm_rooting_rank = gnm_o
                    elif rooting_rank == 'f':
                        gnm_rooting_rank = gnm_f
                    elif rooting_rank == 'g':
                        gnm_rooting_rank = gnm_g

                    # rooting_rank
                    if gnm_rooting_rank in rooting_rank_taxon_dict:
                        if gnm_rooting_rank not in counted_taxons_rooting_rank:
                            counted_taxons_rooting_rank.add(gnm_rooting_rank)
                            canditate_gnms_rooting_rank[ref_accession] = gnm_rooting_rank

    ref_tree_rooting_rank = subset_and_rename_tree(gtdb_ref_tree, canditate_gnms_rooting_rank, canditate_gnms_rooting_rank)

    # get the smallest out group taxon set
    smallest_outgroup_taxon_list = get_smallest_outgroup(ref_tree_rooting_rank)

    user_gnm_taxon_dict_rooting_rank = dict()
    if rooting_rank == 'p':
        user_gnm_taxon_dict_rooting_rank = user_gnm_taxon_dict_p
    elif rooting_rank == 'c':
        user_gnm_taxon_dict_rooting_rank = user_gnm_taxon_dict_c
    elif rooting_rank == 'o':
        user_gnm_taxon_dict_rooting_rank = user_gnm_taxon_dict_o
    elif rooting_rank == 'f':
        user_gnm_taxon_dict_rooting_rank = user_gnm_taxon_dict_f
    elif rooting_rank == 'g':
        user_gnm_taxon_dict_rooting_rank = user_gnm_taxon_dict_g

    # get the smallest out group genome set
    out_group_gnm_set_1 = set()
    out_group_gnm_set_2 = set()
    for each_rooting_rank_taxon in user_gnm_taxon_dict_rooting_rank:
        gnm_member_set = user_gnm_taxon_dict_rooting_rank[each_rooting_rank_taxon]
        if each_rooting_rank_taxon in smallest_outgroup_taxon_list:
            out_group_gnm_set_1.update(gnm_member_set)
        else:
            out_group_gnm_set_2.update(gnm_member_set)

    # select the smaller one as outgroup
    if len(out_group_gnm_set_1) < len(out_group_gnm_set_2):
        out_group_gnm_set = out_group_gnm_set_1
    else:
        out_group_gnm_set = out_group_gnm_set_2

    # root user tree with identified out group genomes
    root_with_outgroup(user_gnm_tree, out_group_gnm_set, user_tree_rooted)


def concat_trees(tree_file_1, tree_file_2, concatenated_tree):

    tree1 = Tree(tree_file_1)
    tree2 = Tree(tree_file_2)
    tree1_with_root = Tree()
    tree2_with_root = Tree()
    tree1_with_root.add_child(tree1)
    tree2_with_root.add_child(tree2)
    concat_tree = Tree()
    concat_tree.add_child(tree1_with_root)
    concat_tree.add_child(tree2_with_root)
    concat_tree.write(outfile=concatenated_tree, format=1)


def RootTree(db_dir, user_gnm_taxon, user_gnm_tree_bac, user_gnm_tree_ar, rooted_tree_bac, rooted_tree_ar, rooted_tree_combined):

    # define file name
    gtdb_ref_tree_ar    = '%s/ar53_r214.tree'           % db_dir
    gtdb_ref_tree_bac   = '%s/bac120_r214.tree'         % db_dir
    gtdb_gnm_meta_ar    = '%s/ar53_metadata_r214.tsv'   % db_dir
    gtdb_gnm_meta_bac   = '%s/bac120_metadata_r214.tsv' % db_dir

    contain_ar_gnm = False
    contain_bac_gnm = False
    for each_genome in open(user_gnm_taxon):
        if not each_genome.startswith('user_genome'):
            each_split = each_genome.strip().split('\t')
            taxon_str = each_split[1]
            if 'd__Bacteria' in taxon_str:
                contain_bac_gnm = True
            if 'd__Archaea' in taxon_str:
                contain_ar_gnm = True

    # root bacterial tree
    if (contain_bac_gnm is True) and (contain_ar_gnm is False):
        root_tree_by_gtdb(gtdb_ref_tree_bac, gtdb_gnm_meta_bac, user_gnm_taxon, user_gnm_tree_bac, rooted_tree_combined,'bac')

    # root archaeal tree
    elif (contain_bac_gnm is False) and (contain_ar_gnm is True):
        root_tree_by_gtdb(gtdb_ref_tree_ar,  gtdb_gnm_meta_ar,  user_gnm_taxon, user_gnm_tree_ar,  rooted_tree_combined,'ar')

    # root both trees
    elif (contain_bac_gnm is True) and (contain_ar_gnm is True):
        root_tree_by_gtdb(gtdb_ref_tree_bac, gtdb_gnm_meta_bac, user_gnm_taxon, user_gnm_tree_bac, rooted_tree_bac,'bac')
        root_tree_by_gtdb(gtdb_ref_tree_ar,  gtdb_gnm_meta_ar,  user_gnm_taxon, user_gnm_tree_ar,  rooted_tree_ar, 'ar')
        concat_trees(rooted_tree_bac, rooted_tree_ar, rooted_tree_combined)



db_dir                  = '/Users/songweizhi/DB/GTDB_r214'
tax_file_to_use         = '/Users/songweizhi/Desktop/999/NorthSea_69_MAGs.bac120.summary.tsv'
inferred_bac120_tree    = '/Users/songweizhi/Desktop/999/bac120.unrooted.tree'
inferred_ar53_tree      = ''
user_tree_rooted_bac    = '/Users/songweizhi/Desktop/999/rooted_species_tree_bac.tree'
user_tree_rooted_ar     = '/Users/songweizhi/Desktop/999/rooted_species_tree_ar.tree'
user_tree_rooted        = '/Users/songweizhi/Desktop/999/rooted_species_tree.tree'


db_dir                  = '/scratch/PI/ocessongwz/DB/gtdbtk_db'
tax_file_to_use         = '/scratch/PI/boqianpy/shanzhang/01_MetaCHIP2/02_NorthSea/NorthSea_69_MAGs.bac120.summary.tsv'
inferred_bac120_tree    = 'bac120.unrooted.tree'
inferred_ar53_tree      = ''
user_tree_rooted_bac    = 'rooted_species_tree_bac.tree'
user_tree_rooted_ar     = 'rooted_species_tree_ar.tree'
user_tree_rooted        = 'rooted_species_tree.tree'


RootTree(db_dir, tax_file_to_use, inferred_bac120_tree, inferred_ar53_tree, user_tree_rooted_bac, user_tree_rooted_ar, user_tree_rooted)
