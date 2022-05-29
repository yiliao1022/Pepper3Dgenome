library(CALDER)

contact_mat_file_LJTZ2 = "/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/100k/Chr12/matrixTZ2.Chr12.csv.txt"
CALDER_main(contact_mat_file_LJTZ2, chr=12, bin_size=10E4, out_dir='./Chr12_LJTZ2', sub_domains=TRUE, save_intermediate_data=FALSE)
