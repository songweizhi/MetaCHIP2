import os
import argparse


circos_usage = '''
======================= circos example commands =======================

MetaCHIP2 circos -o circos_1.pdf -i cir_plot_matrix.csv
MetaCHIP2 circos -o circos_2.pdf -l detected_HGTs.txt -g grouping.txt

=======================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def circos_plot_from_hgt_and_grouping(detected_hgts_txt, grouping_file, circos_HGT_R, pwd_plot_circos):

    # get genome to group dict
    genome_to_group_dict = {}
    for genome in open(grouping_file):
        group_id2 = genome.strip().split(',')[0]
        genome_name = genome.strip().split(',')[1]
        genome_to_group_dict[genome_name] = group_id2

    # define file name
    f_path, f_base, f_ext = sep_path_basename_ext(pwd_plot_circos)
    pwd_cir_plot_t1                 = '%s/%s_matrix_t1.txt'                % (f_path, f_base)
    pwd_cir_plot_t1_sorted          = '%s/%s_matrix_t1_sorted.txt'         % (f_path, f_base)
    pwd_cir_plot_t1_sorted_count    = '%s/%s_matrix_t1_sorted_count.txt'   % (f_path, f_base)
    pwd_cir_plot_matrix_filename    = '%s/%s_matrix.txt'                   % (f_path, f_base)

    transfers = []
    col_index = {}
    for each in open(detected_hgts_txt):
        each_split = each.strip().split('\t')
        if each.startswith('Gene_1\tGene_2\tIdentity'):
            col_index = {key: i for i, key in enumerate(each_split)}
        else:
            Direction = each_split[col_index['direction']]
            if '%)' in Direction:
                Direction = Direction.split('(')[0]

            transfers.append(Direction)

    tmp1 = open(pwd_cir_plot_t1, 'w')
    all_group_id = []
    for each_t in transfers:
        each_t_split    = each_t.split('-->')
        donor           = each_t_split[0]
        recipient       = each_t_split[1]
        donor_group     = genome_to_group_dict[donor]
        recipient_group = genome_to_group_dict[recipient]
        if donor_group not in all_group_id:
            all_group_id.append(donor_group)
        if recipient_group not in all_group_id:
            all_group_id.append(recipient_group)
        tmp1.write('%s,%s\n' % (donor_group, recipient_group))
    tmp1.close()

    os.system('cat %s | sort > %s' % (pwd_cir_plot_t1, pwd_cir_plot_t1_sorted))

    current_t = ''
    count = 0
    tmp2 = open(pwd_cir_plot_t1_sorted_count, 'w')
    for each_t2 in open(pwd_cir_plot_t1_sorted):
        each_t2 = each_t2.strip()
        if current_t == '':
            current_t = each_t2
            count += 1
        elif current_t == each_t2:
            count += 1
        elif current_t != each_t2:
            tmp2.write('%s,%s\n' % (current_t, count))
            current_t = each_t2
            count = 1
    tmp2.write('%s,%s\n' % (current_t, count))
    tmp2.close()

    # read in count as dict
    transfer_count = {}
    for each_3 in open(pwd_cir_plot_t1_sorted_count):
        each_3_split = each_3.strip().split(',')
        key = '%s,%s' % (each_3_split[0], each_3_split[1])
        value = each_3_split[2]
        transfer_count[key] = value

    all_group_id = sorted(all_group_id)

    matrix_file = open(pwd_cir_plot_matrix_filename, 'w')
    matrix_file.write('\t' + '\t'.join(all_group_id) + '\n')
    for each_1 in all_group_id:
        row = [each_1]
        for each_2 in all_group_id:
            current_key = '%s,%s' % (each_2, each_1)
            if current_key not in transfer_count:
                row.append('0')
            else:
                row.append(transfer_count[current_key])
        matrix_file.write('\t'.join(row) + '\n')
    matrix_file.close()

    # get plot with R
    if len(all_group_id) > 1:
        r_cmd = 'Rscript %s -m %s -s 11 -p %s' % (circos_HGT_R, pwd_cir_plot_matrix_filename, pwd_plot_circos)
        print(r_cmd)
        os.system(r_cmd)

    # rm tmp files
    os.system('rm %s' % pwd_cir_plot_t1)
    os.system('rm %s' % pwd_cir_plot_t1_sorted)
    os.system('rm %s' % pwd_cir_plot_t1_sorted_count)


def circos(args):

    input_matrix        = args['i']
    detected_hgts_txt   = args['l']
    grouping_file       = args['g']
    font_size           = args['s']
    output_plot         = args['o']

    # get path to circos_HGT.R
    current_file_path   = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    circos_HGT_R        = '%s/circos_HGT.R' % current_file_path

    # plot with data matrix
    if (input_matrix is not None) and (detected_hgts_txt is None) and (grouping_file is None):
        os.system('Rscript %s -m %s -s %s -p %s' % (circos_HGT_R, input_matrix, font_size, output_plot))
        print('Done, plot exported to %s' % output_plot)

    # plot with detected_hgts_txt and grouping_file
    elif (input_matrix is None) and (detected_hgts_txt is not None) and (grouping_file is not None):
        circos_plot_from_hgt_and_grouping(detected_hgts_txt, grouping_file, circos_HGT_R, output_plot)
    else:
        print('Please refers to the example commands for usage, program exited!')
        exit()

    print('Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',   required=False, default=None,         help='input matrix')
    parser.add_argument('-l',   required=False, default=None,         help='MetaCHIP produced detected_HGTs.txt')
    parser.add_argument('-g',   required=False, default=None,         help='grouping file')
    parser.add_argument('-s',   required=False, type=int, default=12, help='font size, default: 12')
    parser.add_argument('-o',   required=True,                        help='output plot')
    args = vars(parser.parse_args())
    circos(args)


'''

cd /Users/songweizhi/Desktop
MetaCHIP2 circos -o circos_f.pdf -s 11 -l detected_HGTs.txt -g Bacilli_plus_78_clade_f.txt
MetaCHIP2 circos -o circos_g.pdf -s 11 -l detected_HGTs.txt -g Bacilli_plus_78_clade.txt
MetaCHIP2 circos -o circos_f2.pdf -s 11 -i circos_f_matrix.txt
MetaCHIP2 circos -o circos_g2.pdf -s 11 -i circos_g_matrix.txt

cd /Users/songweizhi/Desktop
Rscript /Users/songweizhi/PycharmProjects/MetaCHIP2/MetaCHIP2/circos_HGT.R -m cir_plot_matrix.csv -p cir_plot_matrix.pdf -s 11

'''
