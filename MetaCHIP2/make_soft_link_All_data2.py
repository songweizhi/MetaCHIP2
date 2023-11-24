#!/bin/python

import os

dst = f'/home-user/lling/4_Population/222PB_strains/222Protein/'

file = "/home-user/lling/4_Population/222PB_strains/All_PB_fullname.txt"
file = '/Users/songweizhi/Desktop/All_PB_fullname.txt'
f1 = open(file, 'r')

for i in f1.readlines():
    i_split = i.strip().split('_')
    print(i_split)
    gnm_id = i_split[-1]

    src = f'/home-user/lling/4_Population/222PB_strains/protein2/{i.strip()}.faa'
    print(src)
    # si = i.strip("\n").split("_", )[-1]
    # gi = si[-1]
    #
    # src = f'/home-user/lling/4_Population/222PB_strains/protein2/{i}.faa'
    #
    # if "HKCCYL" in si:
    #     print(i)
    #     print(si)
    #     print(gi)
    #     ln_cmd = 'ln -s %s %s' % (src, '')
    #     print(ln_cmd)
    #     # print(src, dst + gi.strip('\n') + '.faa')
        #os.symlink(src, dst + gi.strip('\n') + '.faa')
    #
    # elif "SZCCHN" in gi:
    #     print(src, dst + gi.strip('\n') + '.faa')
    #     os.symlink(src, dst + gi.strip('\n') + '.faa')
    #
    # else:
    #     print(src, dst + i + '.faa')
    #     os.symlink(src, dst + i + '.faa')



'''

ln -s source_file symbolic_link

'''