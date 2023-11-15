
Commands to get a plot like Fig. 9 in the [MetaCHIP paper](https://doi.org/10.1186/s40168-019-0649-y)

    module load python/3.7.4
    source ~/mypython3env/bin/activate
    module load diamond/2.1.6
    cd /srv/scratch/z5265700/NASA
    BioSAK COG2020 -m P -t 12 -db_dir /srv/scratch/z5265700/DB/COG2020 -i faa_files -x faa -diamond
    # BioSAK COG2020 -m P -t 30 -db_dir /srv/scratch/z5265700/DB/COG2020 -i faa_files -x faa -diamond
    
    
    # plot
    cd /srv/scratch/z5265700/Shan_z5095298/z5095298/Weizhi/NASA
    python3 /srv/scratch/z5265700/Shan_z5095298/z5095298/Weizhi/Scripts/boxplot_matrix_COG.py -in faa_files_COG2020_fun_stats -out faa_files_COG2020_fun_stats.txt -in_percent -skip_1st_row
    
    cd /Users/songweizhi/Desktop
    #Rscript /Users/songweizhi/PycharmProjects/NASA/COG_boxplot_last1row.R -i faa_files_COG2020_fun_stats.txt -o faa_files_COG2020_fun_stats.png
    Rscript /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/COG_boxplot_last1row.R -i faa_files_COG2020_fun_stats_blastn.txt -o faa_files_COG2020_fun_stats_blastn.png
    Rscript /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/COG_boxplot_last1row.R -i faa_files_COG2020_fun_stats_mmseqs.txt -o faa_files_COG2020_fun_stats_mmseqs.png
