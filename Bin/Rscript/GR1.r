library(CALDER)

contact_mat_file_LJGR1 = "/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/100k/Chr12/matrixGR1.Chr12.csv.txt"
CALDER_main(contact_mat_file_LJGR1, chr=12, bin_size=10E4, out_dir='./Chr12_LJGR1', sub_domains=TRUE, save_intermediate_data=FALSE)
