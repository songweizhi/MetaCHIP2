import os
import argparse


infer_usage = '''
====================== infer example command ======================

MetaCHIP2 infer -o Demo -i gnm_folder -x fa -t 12

Note: To run the infer module, you need to have gtdbtk installed 
on your system. 

===================================================================
'''


def infer(args):

    input_gnm_dir       = args['i']
    file_extension      = args['x']
    output_dir          = args['o']
    num_threads         = args['t']
    force_overwrite     = args['f']

    # create output folder
    if (os.path.isdir(output_dir) is True) and (force_overwrite is False):
        print('Output folder detected, program exited: %s' % output_dir)
        exit()
    else:
        if os.path.isdir(output_dir) is True:
            os.system('rm -r %s' % output_dir)
        os.mkdir(output_dir)

    # define file name
    msa_bac120_gz       = '%s/align/gtdbtk.bac120.user_msa.fasta.gz'                            % output_dir
    msa_bac120          = '%s/align/gtdbtk.bac120.user_msa.fasta'                               % output_dir
    msa_ar53_gz         = '%s/align/gtdbtk.ar53.user_msa.fasta.gz'                              % output_dir
    msa_ar53            = '%s/align/gtdbtk.ar53.user_msa.fasta'                                 % output_dir

    # prepare commands
    cmd_identify        = 'gtdbtk identify --genome_dir %s -x %s --out_dir %s --cpus %s'        % (input_gnm_dir, file_extension, output_dir, num_threads)
    cmd_align           = 'gtdbtk align --identify_dir %s --out_dir %s --cpus %s'               % (output_dir, output_dir, num_threads)
    cmd_gunzip_bac120   = 'gunzip %s'                                                           % msa_bac120_gz
    cmd_gunzip_ar53     = 'gunzip %s'                                                           % msa_ar53_gz
    cmd_infer_bac120    = 'gtdbtk infer --msa_file %s --out_dir %s --cpus %s --prefix bac120'   % (msa_bac120, output_dir, num_threads)
    cmd_infer_ar53      = 'gtdbtk infer --msa_file %s --out_dir %s --cpus %s --prefix ar53'     % (msa_ar53, output_dir, num_threads)

    print(cmd_identify)
    os.system(cmd_identify)
    print(cmd_align)
    os.system(cmd_align)

    if os.path.isfile(msa_bac120_gz):
        print(cmd_gunzip_bac120)
        os.system(cmd_gunzip_bac120)
        print(cmd_infer_bac120)
        os.system(cmd_infer_bac120)

    if os.path.isfile(msa_ar53_gz):
        print(cmd_gunzip_ar53)
        os.system(cmd_gunzip_ar53)
        print(cmd_infer_ar53)
        os.system(cmd_infer_ar53)

    inferred_bac120_tree = '%s/bac120.unrooted.tree' % output_dir
    inferred_ar53_tree   = '%s/ar53.unrooted.tree'   % output_dir

    if os.path.isfile(inferred_bac120_tree):
        print('Inferred bacterial tree:\t%s' % inferred_bac120_tree)
    if os.path.isfile(inferred_ar53_tree):
        print('Inferred archaeal tree:\t%s'  % inferred_ar53_tree)

    print('Done!')


if __name__ == '__main__':

    infer_parser = argparse.ArgumentParser(usage=infer_usage)
    infer_parser.add_argument('-i', required=True,                       help='genome folder')
    infer_parser.add_argument('-x', required=False, default='fna',       help='genome file extension, default: fna')
    infer_parser.add_argument('-o', required=True,                       help='output folder')
    infer_parser.add_argument('-t', required=False, type=int, default=1, help='number of threads')
    infer_parser.add_argument('-f', required=False, action="store_true", help='force overwrite existing results')
    args = vars(infer_parser.parse_args())
    infer(args)
