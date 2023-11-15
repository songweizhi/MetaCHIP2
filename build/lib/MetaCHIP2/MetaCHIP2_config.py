import os


# extract path to this config file
pwd_config_file = os.path.realpath(__file__)
config_file_path = '/'.join(pwd_config_file.split('/')[:-1])

# specify full path to corresponding executables at the right side of colon
config_dict = {'config_file_path'       : config_file_path,
               'mafft'           : 'mafft',
               'blastp'          : 'blastp',
               'blastn'          : 'blastn',
               'makeblastdb'     : 'makeblastdb',
               'fasttree'        : 'FastTree',
               'ranger_mac'      : '%s/Ranger-DTL-Dated.mac'    % config_file_path,  # do not edit this line
               'ranger_linux'    : '%s/Ranger-DTL-Dated.linux'  % config_file_path,  # do not edit this line
               'circos_HGT_R'    : '%s/circos_HGT.R'            % config_file_path   # do not edit this line
               }
