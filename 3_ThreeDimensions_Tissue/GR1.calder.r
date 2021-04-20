library(CALDER)

contact_mat_file_LJGR1 = "/home/yiliao/OLDISK/genome_assembly/hic_explorer/7_hic2cool/40k/cooltools_raw/ajust/Final/Chr01/matrixGR1.chr01.csv.txt"
CALDER_main(contact_mat_file_LJGR1, chr=1, bin_size=40E3, out_dir='./Chr01_LJGR1', sub_domains=TRUE, save_intermediate_data=FALSE)
