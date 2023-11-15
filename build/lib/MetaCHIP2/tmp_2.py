

kegg_cnt_match_PB2_C2_tsv                   = '/Users/songweizhi/Desktop/000/kegg_cnt_match_PB2_C2.tsv'
gnm_id                                      = 'Bradyrhizobium_sp_ORS_285'
gnm_id                                      = 'Bradyrhizobium_sp_ORS_278_ORS278'
kegg_cnt_match_PB2_C2_tsv_with_locus_tag    = '/Users/songweizhi/Desktop/000/kegg_cnt_match_PB2_C2_%s.tsv'  % gnm_id


ko_list = []
for each_ko in open(kegg_cnt_match_PB2_C2_tsv):
    ko_list.append(each_ko.split('\t')[0])

ko_to_locus_tag_dict = dict()
for each_ko in open('/Users/songweizhi/Desktop/000/seq_info.tbl'):
    each_ko_split = each_ko.strip().split('\t')
    ko_id = each_ko_split[0]
    if ko_id in ko_list:
        for each_gnm in each_ko_split[1:]:
            if gnm_id in each_gnm:
                each_gnm_split = each_gnm.split(':')
                for each_g2 in each_gnm_split[1:]:
                    ko_to_locus_tag_dict[ko_id] = each_g2

kegg_cnt_match_PB2_C2_tsv_with_locus_tag_handle = open(kegg_cnt_match_PB2_C2_tsv_with_locus_tag, 'w')
for each_ko in open(kegg_cnt_match_PB2_C2_tsv):
    ko_id = each_ko.split('\t')[0]
    if ko_id == 'KEGGid':
        locus_tag_str = 'locus_tag'
    else:
        locus_tag_str = ko_to_locus_tag_dict.get(ko_id, 'na')
    kegg_cnt_match_PB2_C2_tsv_with_locus_tag_handle.write('%s\t%s\n' % (each_ko.strip(), locus_tag_str))
kegg_cnt_match_PB2_C2_tsv_with_locus_tag_handle.close()
