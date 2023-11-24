import os
import argparse
import pandas as pd
from pycirclize import Circos
import matplotlib as mpl
mpl.use('Agg')


circos_usage = '''
======================= circos example commands =======================

MetaCHIP2 circos -o circos_1.pdf -i cir_plot_matrix.csv
MetaCHIP2 circos -o circos_2.pdf -l detected_HGTs.txt -g grouping.txt

# grouping.txt (tab separated)
S2_sal_fde_bin1	c__Clostridia
S2_sal_fde_bin2	c__Bacteroidia
S2_sal_fde_bin6	c__Bacilli

=======================================================================
'''


def get_circos_matrix(grouping_file, detected_hgts_txt, hgt_matrix):

    # get genome to group dict
    gnm_to_group_dict = dict()
    for genome in open(grouping_file):
        gnm_id = genome.strip().split('\t')[0]
        group_id = genome.strip().split('\t')[1]
        gnm_to_group_dict[gnm_id] = group_id

    grp_set = set()
    col_index = {}
    d2r_hgt_num_dict = dict()
    for each in open(detected_hgts_txt):
        each_split = each.strip().split('\t')
        if each.startswith('Gene_1\tGene_2\tIdentity'):
            col_index = {key: i for i, key in enumerate(each_split)}
        else:
            direction = each_split[col_index['direction']]
            if '%)' in direction:
                direction = direction.split('(')[0]

            direction_split = direction.split('-->')
            gnm_d = direction_split[0]
            gnm_r = direction_split[1]
            grp_d = gnm_to_group_dict[gnm_d]
            grp_r = gnm_to_group_dict[gnm_r]
            grp_set.add(grp_d)
            grp_set.add(grp_r)
            key_grp_d_to_r = '%s_d2r_%s' % (grp_d, grp_r)
            if key_grp_d_to_r not in d2r_hgt_num_dict:
                d2r_hgt_num_dict[key_grp_d_to_r] = 1
            else:
                d2r_hgt_num_dict[key_grp_d_to_r] += 1

    grp_list_sorted = sorted([i for i in grp_set])

    hgt_matrix_handle = open(hgt_matrix, 'w')
    hgt_matrix_handle.write('\t' + '\t'.join(grp_list_sorted) + '\n')
    for each_d_grp in grp_list_sorted:
        num_list = [each_d_grp]
        for each_r_grp in grp_list_sorted:
            key_d_to_r = '%s_d2r_%s' % (each_d_grp, each_r_grp)
            hgt_num = d2r_hgt_num_dict.get(key_d_to_r, 0)
            num_list.append(str(hgt_num))
        hgt_matrix_handle.write('\t'.join(num_list) + '\n')
    hgt_matrix_handle.close()


def pycircos(data_matrix, sep_symbol, plot_out):

    matrix_df = pd.read_csv(data_matrix, sep=sep_symbol, header=0, index_col=0)
    interval = round((matrix_df.max()).max()/10)*5                  # get tick interval
    circos = Circos.initialize_from_matrix(matrix_df,
                                           start=-265,              # Plot start degree (-360 <= start < end <= 360)
                                           end=95,                  # Plot end degree (-360 <= start < end <= 360)
                                           space=1,                 # Space degree(s) between sector
                                           r_lim=(90, 95),          # Outer track radius limit region (0 - 100)
                                           cmap="tab10",            # Colormap assigned to each outer track and link.
                                           # order='desc',          # asc, desc; sort in ascending(or descending) order by node size.
                                           ticks_interval=interval, # Ticks interval. If None, ticks are not plotted.
                                           ticks_kws=dict(label_size=3.5, label_orientation="vertical"),    # font size of tick labels
                                           label_kws=dict(size=6, orientation="vertical"),
                                           link_kws=dict(direction=1, color='white', ec="black", lw=0))
    fig = circos.plotfig()
    fig.savefig(plot_out, dpi=100)


def circos(args):

    input_matrix        = args['i']
    detected_hgts_txt   = args['l']
    grouping_file       = args['g']
    output_plot         = args['o']

    # plot with data matrix
    if (input_matrix is not None) and (detected_hgts_txt is None) and (grouping_file is None):
        pycircos(input_matrix, '\t', output_plot)

    # plot with detected_hgts_txt and grouping_file
    elif (input_matrix is None) and (detected_hgts_txt is not None) and (grouping_file is not None):
        hgt_matrix = '%s.matrix.txt' % output_plot
        get_circos_matrix(grouping_file, detected_hgts_txt, hgt_matrix)
        pycircos(hgt_matrix, '\t', output_plot)

    else:
        print('Please refers to the example commands for usage, program exited!')
        exit()

    print('Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',   required=False, default=None,         help='input matrix')
    parser.add_argument('-l',   required=False, default=None,         help='MetaCHIP produced detected_HGTs.txt')
    parser.add_argument('-g',   required=False, default=None,         help='grouping file')
    parser.add_argument('-o',   required=True,                        help='output plot')
    args = vars(parser.parse_args())
    circos(args)
