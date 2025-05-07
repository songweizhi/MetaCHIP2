import os
import glob
import warnings
warnings.filterwarnings("ignore")


pwd_blast_result_filtered_folder = '/Users/songweizhi/Desktop/MetaCHIP2/MetaBAT_1634_138/blastn_results_filtered'
blast_result_filtered_file_re = '%s/*_filtered.tab' % pwd_blast_result_filtered_folder
blast_result_filtered_file_list = [os.path.basename(file_name) for file_name in glob.glob(blast_result_filtered_file_re)]

for filtered_blast_result in blast_result_filtered_file_list:
    genome_id = filtered_blast_result.split('_blastn_filtered.tab')[0]
    print(genome_id)
    pwd_file = '%s/%s' % (pwd_blast_result_filtered_folder, filtered_blast_result)
    for each_identity in open(pwd_file):
        each_identity_split = each_identity.strip().split('\t')
        query = each_identity_split[0]
        subject = each_identity_split[1]
        identity = float(each_identity_split[2])
        query_genome_name = '_'.join(query.split('_')[:-1])
        subject_genome_name = '_'.join(subject.split('_')[:-1])
        print(each_identity_split)
        print(query_genome_name)
        print(subject_genome_name)


















