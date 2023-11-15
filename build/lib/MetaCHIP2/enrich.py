import os
import glob
import pandas
import argparse
import numpy as np
from Bio import SeqIO
import multiprocessing as mp
from datetime import datetime
from Bio.SeqRecord import SeqRecord


enrich_usage = '''
============================ enrich example commands ============================

# This module was prepared to produce a plot as Fig. 9 in the MetaCHIP paper
MetaCHIP2 enrich -faa faa_files -o demo_1HGT -db_dir /Users/songweizhi/DB/COG2020 -t 10 -diamond -f -hgt1 detected_HGTs.faa -desc
MetaCHIP2 enrich -faa faa_files -o demo_2HGT -db_dir /Users/songweizhi/DB/COG2020 -t 10 -diamond -f -hgt1 HGTs_setting1.faa -hgt2 HGTs_setting2.faa

# Prepare DB files (version 2020):
cd path/to/your/COG_db_dir
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.fa.gz
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/Readme.2020-11-25.txt
gunzip cog-20.fa.gz
makeblastdb -in cog-20.fa -dbtype prot -parse_seqids -logfile cog-20.fa.log
diamond makedb --in cog-20.fa --db cog-20.fa.dmnd --quiet

# The complete descriptions of COG categories can be found here:
https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab

=================================================================================
'''

time_format = '[%Y-%m-%d %H:%M:%S] '

def sep_path_basename_ext(file_in):
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


def dna2aa(dna_file, aa_file):
    query_aa_handle = open(aa_file, 'w')
    for each in SeqIO.parse(dna_file, 'fasta'):
        each_aa = each.seq.translate()
        each_aa_record = SeqRecord(each_aa)
        each_aa_record.id = each.id
        each_aa_record.description = each.description
        SeqIO.write(each_aa_record, query_aa_handle, 'fasta')
    query_aa_handle.close()


def best_hit(args):

    file_in = args['i']
    file_out = args['o']
    file_out_handle = open(file_out, 'w')
    best_hit_line = ''
    best_hit_query_id = ''
    best_hit_score = 0
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = float(blast_hit_split[11])

        if best_hit_query_id == '':
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score
        elif (query_id == best_hit_query_id) and (bit_score > best_hit_score):
            best_hit_score = bit_score
            best_hit_line = blast_hit
        elif query_id != best_hit_query_id:
            file_out_handle.write(best_hit_line)
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score
    file_out_handle.write(best_hit_line)
    file_out_handle.close()


def COG2020_worker(argument_list):

    pwd_input_file =                    argument_list[0]
    pwd_prot2003_2014 =                 argument_list[1]
    protein_id_to_cog_id_dict =         argument_list[2]
    cog_id_to_category_dict =           argument_list[3]
    cog_id_to_description_dict =        argument_list[4]
    cog_category_list =                 argument_list[5]
    cog_category_to_description_dict =  argument_list[6]
    output_folder =                     argument_list[7]
    thread_num =                        argument_list[8]
    run_diamond =                       argument_list[9]
    evalue_cutoff =                     argument_list[10]

    input_seq_no_path, input_seq_no_ext, input_seq_ext = sep_path_basename_ext(pwd_input_file)
    current_output_folder = '%s/%s_COG2020_wd' % (output_folder, input_seq_no_ext)

    pwd_blastp_output =             '%s/%s_blastp.tab'                  % (current_output_folder, input_seq_no_ext)
    pwd_blastp_output_besthits =    '%s/%s_blastp_besthits.tab'         % (current_output_folder, input_seq_no_ext)
    pwd_query_to_cog_txt =          '%s/%s_query_to_cog.txt'            % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_GeneNumber =     '%s/%s_func_stats_GeneNumber.txt'   % (current_output_folder, input_seq_no_ext)
    os.mkdir(current_output_folder)

    # run blastp
    if run_diamond is False:
        os.system('blastp -query %s -db %s -out %s -evalue %s -outfmt 6 -show_gis -num_threads %s' % (pwd_input_file, pwd_prot2003_2014, pwd_blastp_output, evalue_cutoff, thread_num))
    else:
        os.system('diamond blastp -q %s --db %s.dmnd --out %s --evalue %s --outfmt 6 --threads %s --quiet' % (pwd_input_file, pwd_prot2003_2014, pwd_blastp_output, evalue_cutoff, thread_num))

    # keep only best hits
    best_hit({'i': pwd_blastp_output, 'o': pwd_blastp_output_besthits})

    # get query_to_ref_protein_dict
    query_to_ref_protein_dict = {}
    for each_hit in open(pwd_blastp_output_besthits):
        each_hit_split = each_hit.strip().split('\t')
        each_hit_query = each_hit_split[0]
        each_hit_subject = each_hit_split[1]
        each_hit_subject_no_dot = '_'.join(each_hit_subject.split('.'))
        query_to_ref_protein_dict[each_hit_query] = each_hit_subject_no_dot

    # get query sequences list
    query_seq_list = []
    for query_seq in SeqIO.parse(pwd_input_file, 'fasta'):
        query_seq_list.append(query_seq.id)

    # export annotation
    genes_with_cog = set()
    cog_id_num_dict = {}
    cog_cate_num_dict = {}
    cog_id_to_gene_member_dict = {}
    cog_cate_to_gene_member_dict = {}
    pwd_query_to_cog_txt_handle = open(pwd_query_to_cog_txt, 'w')
    pwd_query_to_cog_txt_handle.write('Query\tCOG\tCategory\tDescription\n')
    for query_gene in sorted(query_seq_list):
        if query_gene not in query_to_ref_protein_dict:
            pwd_query_to_cog_txt_handle.write('%s\n' % (query_gene))
        else:
            db_protein_id = query_to_ref_protein_dict[query_gene]
            if db_protein_id not in protein_id_to_cog_id_dict:
                pwd_query_to_cog_txt_handle.write('%s\n' % (query_gene))
            else:
                cog_id_list = protein_id_to_cog_id_dict[db_protein_id]
                for cog_id in cog_id_list:
                    cog_cate = cog_id_to_category_dict[cog_id]
                    cog_des = cog_id_to_description_dict[cog_id]
                    pwd_query_to_cog_txt_handle.write('%s\t%s\t%s\t%s\n' % (query_gene, cog_id, cog_cate, cog_des))
                    genes_with_cog.add(query_gene)

                    # update cog_id_num_dict
                    if cog_id not in cog_id_num_dict:
                        cog_id_num_dict[cog_id] = 1
                        cog_id_to_gene_member_dict[cog_id] = [query_gene]
                    else:
                        cog_id_num_dict[cog_id] += 1
                        cog_id_to_gene_member_dict[cog_id].append(query_gene)

                    # update cog_cate_num_dict
                    for each_cog_cate in cog_cate:
                        if each_cog_cate not in cog_cate_num_dict:
                            cog_cate_num_dict[each_cog_cate] = 1
                            cog_cate_to_gene_member_dict[each_cog_cate] = [query_gene]
                        else:
                            cog_cate_num_dict[each_cog_cate] += 1
                            cog_cate_to_gene_member_dict[each_cog_cate].append(query_gene)
    pwd_query_to_cog_txt_handle.close()

    #################### export func_stats_GeneNumber ####################

    pwd_func_stats_GeneNumber_handle = open(pwd_func_stats_GeneNumber, 'w')
    pwd_func_stats_GeneNumber_handle.write('Category\tGeneNumber\tDescription\n')
    for each_cog_cate in cog_category_list:
        each_cog_cate_GeneNumber = 0
        if each_cog_cate in cog_cate_num_dict:
            each_cog_cate_GeneNumber = cog_cate_num_dict[each_cog_cate]
        pwd_func_stats_GeneNumber_handle.write('%s\t%s\t%s\n' % (each_cog_cate, each_cog_cate_GeneNumber, cog_category_to_description_dict[each_cog_cate]))
    pwd_func_stats_GeneNumber_handle.close()


def COG2020(file_in, file_extension, depth_file, DB_dir, num_threads, run_diamond, evalue_cutoff, output_folder):

    pwd_cog_20_fa           = '%s/cog-20.fa'         % DB_dir
    pwd_cog_20_fa_diamond   = '%s/cog-20.fa.dmnd'    % DB_dir
    pwd_cog_20_cog_csv      = '%s/cog-20.cog.csv'    % DB_dir
    pwd_cog_20_def_tab      = '%s/cog-20.def.tab'    % DB_dir
    pwd_fun_20_tab          = '%s/fun-20.tab'        % DB_dir

    ############################################ check whether db file exist ###########################################

    # check whether db file exist
    unfound_inputs = []
    for each_input in [pwd_cog_20_fa, pwd_cog_20_def_tab, pwd_fun_20_tab]:
        if (not os.path.isfile(each_input)) and (not os.path.isdir(each_input)):
            unfound_inputs.append(each_input)
    if len(unfound_inputs) > 0:
        for each_unfound in unfound_inputs:
            print('%s not found' % each_unfound)
        exit()

    if run_diamond is True:
        if os.path.isfile(pwd_cog_20_fa_diamond) is False:
            print('DB file for diamond not found, please refers to the help info for diamond db preparation')
            print('Program exited!')
            exit()

    ################################################# read db into dict ################################################

    # get protein_to_cog_dict (cog-20.cog.csv)
    protein_to_cog_dict = {}
    for each_line in open(pwd_cog_20_cog_csv):
        each_line_split = each_line.strip().split(',')
        protein_id = each_line_split[2]
        protein_id_no_dot = '_'.join(protein_id.split('.'))
        cog_id = each_line_split[6]
        if protein_id_no_dot not in protein_to_cog_dict:
            protein_to_cog_dict[protein_id_no_dot] = {cog_id}
        else:
            protein_to_cog_dict[protein_id_no_dot].add(cog_id)

    # get cog_id_to_category_dict and cog_id_to_description_dict (cognames2003-2014.tab)
    cog_id_to_category_dict = {}
    cog_id_to_description_dict = {}
    for cog_id_to_cate_des in open(pwd_cog_20_def_tab, encoding='windows-1252'):
        if not cog_id_to_cate_des.startswith('#'):
            cog_id_to_cate_des_split = cog_id_to_cate_des.strip().split('\t')
            cog_id = cog_id_to_cate_des_split[0]
            cog_cate = cog_id_to_cate_des_split[1]
            cog_des = cog_id_to_cate_des_split[2]
            cog_id_to_category_dict[cog_id] = cog_cate
            cog_id_to_description_dict[cog_id] = cog_des

    # get cog_category_to_description_dict (fun2003-2014.tab)
    cog_category_list = []
    cog_category_to_description_dict = {}
    for cog_category in open(pwd_fun_20_tab):
        if not cog_category.startswith('#'):
            cog_category_split = cog_category.strip().split('\t')
            cog_category_list.append(cog_category_split[0])
            cog_category_to_description_dict[cog_category_split[0]] = cog_category_split[2][:15]

    ################################################## if input is file ################################################

    # if input is file
    if os.path.isfile(file_in) is True:
        if depth_file is not None:
            if os.path.isfile(depth_file) is False:
                print(datetime.now().strftime(time_format) + 'specified depth file not found, program exited!')
                exit()

        COG2020_worker([file_in, pwd_cog_20_fa, protein_to_cog_dict, cog_id_to_category_dict, cog_id_to_description_dict,
                        cog_category_list, cog_category_to_description_dict, output_folder, num_threads, run_diamond, evalue_cutoff])

    ################################################ if input is folder ################################################

    # if input is folder
    else:
        # check whether input folder exist
        if os.path.isdir(file_in) is False:
            print(datetime.now().strftime(time_format) + 'input folder not found, program exited')
            exit()
        else:
            # check whether input genome exist
            input_file_re = '%s/*.%s' % (file_in, file_extension)
            input_file_name_list = [os.path.basename(file_name) for file_name in glob.glob(input_file_re)]

            if len(input_file_name_list) == 0:
                print(datetime.now().strftime(time_format) + 'input file not found, program exited')
                exit()

            ######################################################### main #########################################################

            print(datetime.now().strftime(time_format) + 'Running COG annotation for %s files with %s cores' % (len(input_file_name_list), num_threads))

            list_for_multiple_arguments_COG = []
            for input_file in input_file_name_list:
                pwd_input_file = '%s/%s' % (file_in, input_file)
                list_for_multiple_arguments_COG.append([pwd_input_file, pwd_cog_20_fa, protein_to_cog_dict, cog_id_to_category_dict, cog_id_to_description_dict,
                                                        cog_category_list, cog_category_to_description_dict, output_folder, 1, run_diamond, evalue_cutoff])

            # run COG annotaion files with multiprocessing
            pool = mp.Pool(processes=num_threads)
            pool.map(COG2020_worker, list_for_multiple_arguments_COG)
            pool.close()
            pool.join()


def turn_to_percentage(number_list):
    number_list_percent = []
    for each_element in number_list:
        each_element_percent = float("{0:.2f}".format(each_element / sum(number_list)))
        number_list_percent.append(each_element_percent)
    return number_list_percent


def boxplot_matrix_COG(input_folder, output_csv, in_percent, skip_1st_row, with_functional_description):

    column_order_list = ['J', 'A', 'K', 'L', 'B', 'D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O', 'X', 'C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q', 'R', 'S']
    input_files = '%s/*.txt' % input_folder
    file_list = [os.path.basename(file_name) for file_name in glob.glob(input_files)]
    file_list = sorted(file_list)

    ############################################ get category_num_lol with dict ############################################

    genome_name_list = []
    detected_category_id_list = set()
    annotation_results_dict = {}
    cate_to_description_dict = {}
    for each_file in file_list:
        genome_name = '_'.join(each_file.split('_')[:-2])
        genome_name_list.append(genome_name)
        current_annotation_results_dict = {}
        n = 0
        for each_category in open('%s/%s' % (input_folder, each_file)):
            each_category_split = each_category.strip().split('\t')
            category_id = each_category_split[0]
            category_num_str = each_category_split[1]
            category_description = each_category_split[2]

            if skip_1st_row is True:
                if n > 0:
                    current_annotation_results_dict[category_id] = int(category_num_str)
                    detected_category_id_list.add(category_id)
            else:
                current_annotation_results_dict[category_id] = int(category_num_str)
                detected_category_id_list.add(category_id)

            if category_id not in cate_to_description_dict:
                cate_to_description_dict[category_id] = category_description

            n += 1

        annotation_results_dict[genome_name] = current_annotation_results_dict

    # reorder columns
    column_order_list_detected = []
    for cate_id in column_order_list:
        if cate_id in detected_category_id_list:
            column_order_list_detected.append(cate_id)

    column_order_list_detected_with_description = []
    for detected_cate in column_order_list_detected:
        column_order_list_detected_with_description.append('%s_%s' % (detected_cate, cate_to_description_dict[detected_cate]))

    # get category_num_lol
    category_num_lol = []
    for genome in genome_name_list:

        current_genome_annotation_results_dict = annotation_results_dict[genome]
        current_genome_annotation_results_list = []
        for cog_cate in column_order_list_detected:
            cog_cate_num = 0
            if cog_cate in current_genome_annotation_results_dict:
                cog_cate_num = current_genome_annotation_results_dict[cog_cate]
            current_genome_annotation_results_list.append(cog_cate_num)

        # turn absolute number to percentage if specified
        if in_percent is True:
            current_genome_annotation_results_list_in_percent = turn_to_percentage(current_genome_annotation_results_list)
            category_num_lol.append(current_genome_annotation_results_list_in_percent)
        else:
            category_num_lol.append(current_genome_annotation_results_list)

    ########################################################################################################################

    # turn list into arrary
    category_num_arrary = np.array(category_num_lol)

    # add row and column name to dataframe
    if with_functional_description is True:
        category_num_df = pandas.DataFrame(category_num_arrary, index=genome_name_list, columns=column_order_list_detected_with_description)
    else:
        category_num_df = pandas.DataFrame(category_num_arrary, index=genome_name_list, columns=column_order_list_detected)

    # write out
    category_num_df.to_csv(output_csv)


def enrich(args):

    op_dir          = args['o']
    file_in         = args['faa']
    file_ext        = args['x']
    hgt1_faa        = args['hgt1']
    hgt2_faa        = args['hgt2']
    db_dir          = args['db_dir']
    num_threads     = args['t']
    run_diamond     = args['diamond']
    evalue_cutoff   = args['e']
    force_overwrite = args['f']
    include_desc    = args['desc']

    current_file_path = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    COG_boxplot_last1row_R = '%s/COG_boxplot_last1row.R' % current_file_path
    COG_boxplot_last2row_R = '%s/COG_boxplot_last2row.R' % current_file_path

    # define file name
    fun_stats_dir           = '%s/fun_stats_files'          % op_dir
    df_file                 = '%s/df_COG2020_fun_stats.txt' % op_dir
    plot_file               = '%s/df_COG2020_fun_stats.png' % op_dir
    cog_annotation_dir_faa  = '%s/COG_gnm'                  % op_dir
    cog_annotation_dir_hgt1 = '%s/COG_hgt1'                 % op_dir
    cog_annotation_dir_hgt2 = '%s/COG_hgt2'                 % op_dir

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('%s exist, program exited!' % op_dir)
            exit()
    os.mkdir(op_dir)

    # annotate gnm
    os.mkdir(cog_annotation_dir_faa)
    COG2020(file_in, file_ext, None, db_dir, num_threads, run_diamond, evalue_cutoff, cog_annotation_dir_faa)

    # annotate hgt1
    os.mkdir(cog_annotation_dir_hgt1)
    COG2020(hgt1_faa, '', None, db_dir, num_threads, run_diamond, evalue_cutoff, cog_annotation_dir_hgt1)

    # annotate hgt2
    if hgt2_faa is not None:
        os.mkdir(cog_annotation_dir_hgt2)
        COG2020(hgt2_faa, '', None, db_dir, num_threads, run_diamond, evalue_cutoff, cog_annotation_dir_hgt2)

    # copy fun_stats files into fun_stats_dir
    os.mkdir(fun_stats_dir)
    cp_cmd_faa  = 'cp %s/*_COG2020_wd/*_func_stats_GeneNumber.txt %s/'                                      % (cog_annotation_dir_faa, fun_stats_dir)
    cp_hgt1_faa = 'cp %s/*_COG2020_wd/*_func_stats_GeneNumber.txt %s/zzzzz_HGT1_func_stats_GeneNumber.txt'  % (cog_annotation_dir_hgt1, fun_stats_dir)
    cp_hgt2_faa = 'cp %s/*_COG2020_wd/*_func_stats_GeneNumber.txt %s/zzzzz_HGT2_func_stats_GeneNumber.txt'  % (cog_annotation_dir_hgt2, fun_stats_dir)
    os.system(cp_cmd_faa)
    os.system(cp_hgt1_faa)
    if hgt2_faa is not None:
        os.system(cp_hgt2_faa)

    # get matrix command
    boxplot_matrix_COG(fun_stats_dir, df_file, True, True, include_desc)

    # get plot with R
    if hgt2_faa is None:
        get_plot_cmd = 'Rscript %s -i %s -o %s' % (COG_boxplot_last1row_R, df_file, plot_file)
    else:
        get_plot_cmd = 'Rscript %s -i %s -o %s' % (COG_boxplot_last2row_R, df_file, plot_file)
    print(get_plot_cmd)
    os.system(get_plot_cmd)


if __name__ == '__main__':

    enrich_parser = argparse.ArgumentParser()
    enrich_parser.add_argument('-o',       required=True,                              help='output plot')
    enrich_parser.add_argument('-faa',     required=True,                              help='faa files of MAGs, produced by Prokka')
    enrich_parser.add_argument('-x',       required=False, default='faa',              help='file extension, default: faa')
    enrich_parser.add_argument('-hgt1',    required=True,                              help='amino acid sequences of HGTs, required')
    enrich_parser.add_argument('-hgt2',    required=False, default=None,               help='amino acid sequences of HGTs, e.g., predicted with a different approach, optional')
    enrich_parser.add_argument('-db_dir',  required=True,                              help='COG_db_dir')
    enrich_parser.add_argument('-diamond', required=False, action='store_true',        help='run diamond (for big dataset), default is NCBI blastp')
    enrich_parser.add_argument('-t',       required=False, type=int, default=1,        help='number of threads')
    enrich_parser.add_argument('-e',       required=False, default=0.001, type=float,  help='evalue cutoff, default: 0.001')
    enrich_parser.add_argument('-desc',    required=False, action='store_true',        help='include functional description in x-axis labels')
    enrich_parser.add_argument('-f',       required=False, action="store_true",        help='force overwrite')
    args = vars(enrich_parser.parse_args())
    enrich(args)


'''

cd /Users/songweizhi/Desktop/enrichHGT
/usr/local/bin/python3.7 /Users/songweizhi/PycharmProjects/MetaCHIP2/MetaCHIP2/enrich.py -faa faa_files -x faa -o demo_1HGT_desc -db_dir /Users/songweizhi/DB/COG2020 -t 10 -diamond -f -desc -hgt1 zHGTs.faa 

cd /Users/songweizhi/Desktop/enrichHGT
/usr/local/bin/python3.7 /Users/songweizhi/PycharmProjects/MetaCHIP2/MetaCHIP2/enrich.py -faa faa_files -x faa -o demo_2HGT_desc -db_dir /Users/songweizhi/DB/COG2020 -t 10 -diamond -f -desc -hgt1 zHGTs.faa -hgt2 zHGTs_mmseqs.faa

'''
