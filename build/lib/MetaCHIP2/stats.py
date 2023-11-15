# import os
# import glob
# from Bio import SeqIO
# import matplotlib.pyplot as plt
# 
#
# def sep_path_basename_ext(file_in):
#
#     f_path, file_name = os.path.split(file_in)
#     if f_path == '':
#         f_path = '.'
#     f_base, f_ext = os.path.splitext(file_name)
#     return f_path, f_base, f_ext

'''
should we indicate whether the gene is essential or not when we discuss the gene been affect by SNV in the discussion section?




e.g., at the begining of the
'''

# def Average(lst):
#     return sum(lst) / len(lst)
#
#
# # grouping_file    = '/Users/songweizhi/Desktop/MetaCHIP2/HGT_Shen/gnm_clade_v1.txt'
# # detected_hgt_txt = '/Users/songweizhi/Desktop/MetaCHIP2/HGT_Shen/demo_MetaCHIP_wd/detected_HGTs.txt'
# # recipient_ffn    = '/Users/songweizhi/Desktop/MetaCHIP2/HGT_Shen/demo_MetaCHIP_wd/recipient.ffn'
# #
# #
# grouping_file    = '/Users/songweizhi/Desktop/MetaCHIP2/HGT_Shen/gnm_clade_v1.txt'
# # detected_hgt_txt = '/Users/songweizhi/Desktop/MetaCHIP2/HGT_Shen/demo_v2_MetaCHIP_wd/detected_HGTs.txt'
# # recipient_ffn    = '/Users/songweizhi/Desktop/MetaCHIP2/HGT_Shen/demo_v2_MetaCHIP_wd/recipient.ffn'
#
# gnm_dir = '/Users/songweizhi/Desktop/MetaCHIP2/HGT_Shen/gnms'
# gnm_ext = 'fna'
#
# ffn_dir = '/Users/songweizhi/Desktop/MetaCHIP2/HGT_Shen/ffn_files'
# ffn_ext = 'ffn'
#
# gnm_to_group_dict = dict()
# for each_genome in open(grouping_file):
#     each_genome_split = each_genome.strip().split(',')
#     group_id = each_genome_split[0]
#     genome_name = each_genome_split[1]
#     gnm_to_group_dict[genome_name] = group_id
#
# # gnm_file_re = '%s/*.%s' % (gnm_dir, gnm_ext)
# # gnm_file_list = glob.glob(gnm_file_re)
# # print(gnm_file_list)
# #
# # grp_to_gc_dict = dict()
# # for each_gnm in gnm_file_list:
# #     _, f_base, _ = sep_path_basename_ext(each_gnm)
# #     gnm_grp = gnm_to_group_dict[f_base]
# #     total_gc_num = 0
# #     total_len = 0
# #     for each_ctg in SeqIO.parse(each_gnm, 'fasta'):
# #         ctg_seq = str(each_ctg.seq).upper()
# #         ctg_len = len(ctg_seq)
# #         ctg_gc_num = ctg_seq.count('G') + ctg_seq.count('C')
# #         total_gc_num += ctg_gc_num
# #         total_len += ctg_len
# #     gnm_gc = total_gc_num*100/total_len
# #
# #
# #     if gnm_grp not in grp_to_gc_dict:
# #         grp_to_gc_dict[gnm_grp] = [gnm_gc]
# #     else:
# #         grp_to_gc_dict[gnm_grp].append(gnm_gc)
# #
# # for each_grp in sorted(list(grp_to_gc_dict.keys())):
# #     value_list = grp_to_gc_dict[each_grp]
# #     mean_gc = Average(value_list)
# #     print('%s\t%s' % (each_grp, mean_gc))
#
# '''
# A	46.94859274256098
# B	47.278857138701355
# C1	52.55834210601024
# C2	51.847601440016284
# C3	52.0317749186422
# D	48.03136944446467
# E	54.98494711794365
# '''
#
#
# ffn_file_re = '%s/*.%s' % (ffn_dir, ffn_ext)
# ffn_file_list = glob.glob(ffn_file_re)
# print(ffn_file_list)
#
#
# grp_to_gene_gc_dict = dict()
# for each_ffn in ffn_file_list:
#     _, gnm_id, _ = sep_path_basename_ext(each_ffn)
#     grp_id = gnm_to_group_dict[gnm_id]
#
#     if grp_id not in grp_to_gene_gc_dict:
#         grp_to_gene_gc_dict[grp_id] = []
#
#     for each_gene in SeqIO.parse(each_ffn, 'fasta'):
#         gene_id = each_gene.id
#         gene_seq = str(each_gene.seq).upper()
#         gene_len = len(gene_seq)
#         gene_gc_num = gene_seq.count('G') + gene_seq.count('C')
#         gene_gc = gene_gc_num*100/gene_len
#         gene_gc = float("{0:.2f}".format(gene_gc))
#         #print('%s\t%s\t%s\t%s' % (grp_id, gnm_id, gene_id, gene_gc))
#         grp_to_gene_gc_dict[grp_id].append(gene_gc)
#
# print(len(grp_to_gene_gc_dict))
# print(grp_to_gene_gc_dict.keys())
#
# fig, axs = plt.subplots(4, 1, sharey=True, tight_layout=True)
#
# # We can set the number of bins with the *bins* keyword argument.
# axs[0].hist(grp_to_gene_gc_dict['A'], bins=100)
# axs[1].hist(grp_to_gene_gc_dict['B'], bins=100)
# axs[2].hist(grp_to_gene_gc_dict['C'], bins=100)
# axs[3].hist(grp_to_gene_gc_dict['D'], bins=100)
# plt.show()
#
# # gnm_to_hgt_num_dict = dict()
# # grp_to_hgt_num_dict = dict()
# # for each_hgt in open(detected_hgt_txt):
# #     if not each_hgt.startswith('Gene_1\tGene_2\tIdentity'):
# #
# #         each_hgt_split = each_hgt.strip().split('\t')
# #         recipient_gnm = each_hgt_split[-1].split('-->')[-1]
# #         recipient_grp = gnm_to_group_dict[recipient_gnm]
# #
# #         if recipient_gnm not in gnm_to_hgt_num_dict:
# #             gnm_to_hgt_num_dict[recipient_gnm] = 1
# #         else:
# #             gnm_to_hgt_num_dict[recipient_gnm] += 1
# #
# #         if recipient_grp not in grp_to_hgt_num_dict:
# #             grp_to_hgt_num_dict[recipient_grp] = 1
# #         else:
# #             grp_to_hgt_num_dict[recipient_grp] += 1
# #
# # print(gnm_to_hgt_num_dict)
# # print(grp_to_hgt_num_dict)
# #
# # grp_to_gc_content_dict = dict()
# # grp_to_gc_num_dict = dict()
# # grp_to_total_len_dict = dict()
# # for each_seq in SeqIO.parse(recipient_ffn, 'fasta'):
# #     recipient_gnm = '_'.join(each_seq.id.split('_')[:-1])
# #     recipient_grp = gnm_to_group_dict[recipient_gnm]
# #     gc_num = str(each_seq.seq).count('G') + str(each_seq.seq).count('C')
# #     gc_content = gc_num*100/len(each_seq.seq)
# #     gc_content = float("{0:.2f}".format(gc_content))
# #     if recipient_grp not in grp_to_gc_content_dict:
# #         grp_to_gc_content_dict[recipient_grp] = [gc_content]
# #         grp_to_gc_num_dict[recipient_grp] = gc_num
# #         grp_to_total_len_dict[recipient_grp] = len(each_seq.seq)
# #     else:
# #         grp_to_gc_content_dict[recipient_grp].append(gc_content)
# #         grp_to_gc_num_dict[recipient_grp] += gc_num
# #         grp_to_total_len_dict[recipient_grp] += len(each_seq.seq)
# #
# #
# # print(grp_to_gc_content_dict)
# #
# #
# # for each_grp in sorted(list(grp_to_gc_content_dict.keys())):
# #     gc_list = grp_to_gc_content_dict[each_grp]
# #     mean_gc = Average(gc_list)
# #     print('%s\t%s\t%s' % (each_grp, mean_gc, gc_list))
# #
# # print()
# # for each_grp in sorted(list(grp_to_gc_content_dict.keys())):
# #     total_gc = grp_to_gc_num_dict[each_grp]/grp_to_total_len_dict[each_grp]
# #     print('%s\t%s' % (each_grp,total_gc))
