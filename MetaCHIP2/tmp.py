import os


def sep_path_basename_ext(file_in):
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


grouping_file    = '/Users/songweizhi/Desktop/MetaCHIP2/HGT_Shen/Bacilli_plus_78_clade.txt'
detected_hgt_txt = '/Users/songweizhi/Desktop/MetaCHIP2/HGT_Shen/ms/detected_HGTs.txt'


gnm_to_group_dict = dict()
for each_genome in open(grouping_file):
    each_genome_split = each_genome.strip().split(',')
    group_id = each_genome_split[0]
    genome_name = each_genome_split[1]
    gnm_to_group_dict[genome_name] = group_id


a = 0
b = 0
c = 0
d = 0
e = 0
f = 0
m = 0
n = 0
for each_hgt in open(detected_hgt_txt):
    if not each_hgt.startswith('Gene_1\tGene_2\tIdentity'):
        each_hgt_split = each_hgt.strip().split('\t')
        d_gnm = each_hgt_split[-1].split('-->')[0]
        r_gnm = each_hgt_split[-1].split('-->')[-1]
        d_grp = gnm_to_group_dict[d_gnm]
        r_grp = gnm_to_group_dict[r_gnm]

        if (d_grp in ['A', 'B', 'C', 'D']) and (r_grp in ['A', 'B', 'C', 'D']) :
            a += 1
        elif (d_grp in ['A', 'B', 'C', 'D']) or (r_grp in ['A', 'B', 'C', 'D']) :
            b += 1
        elif (d_grp not in ['A', 'B', 'C', 'D']) and (r_grp not in ['A', 'B', 'C', 'D']) :
            c += 1
        else:
            d += 1

        if (d_grp not in ['A', 'B', 'C', 'D']) and (r_grp in ['A']) :
            e += 1

        if (d_grp not in ['A', 'B', 'C', 'D']) and (r_grp in ['B']) :
            f += 1

        if (d_grp not in ['A', 'B', 'C', 'D']) and (r_grp in ['C']) :
            m += 1

        if (d_grp not in ['A', 'B', 'C', 'D']) and (r_grp in ['D']) :
            n += 1

print(a)
print(b)
print(c)
print(d)
print(e)
print(f)
print(m)
print(n)
