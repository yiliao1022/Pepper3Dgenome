library(CALDER)

contact_mat_file_LJHL1 = "/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/500k/Chr01/matrixHL1.chr01.csv.txt"
contact_mat_file_LJHL2 = "/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/500k/Chr01/matrixHL2.chr01.csv.txt"
contact_mat_file_LJGR1 = "/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/500k/Chr01/matrixGR1.chr01.csv.txt"
contact_mat_file_LJGR1 = "/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/500k/Chr01/matrixGR2.chr01.csv.txt"
contact_mat_file_LJTZ1 = "/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/500k/Chr01/matrixTZ1.chr01.csv.txt"
contact_mat_file_LJTZ2 = "/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/500k/Chr01/matrixTZ2.chr01.csv.txt"
contact_mat_file_LJYP1 = "/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/500k/Chr01/matrixYP1.chr01.csv.txt"
contact_mat_file_LJYP2 = "/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/500k/Chr01/matrixYP2.chr01.csv.txt"
CALDER_main(contact_mat_file_LJHL1, chr=1, bin_size=50E4, out_dir='./Chr01_LJHL1', sub_domains=TRUE, save_intermediate_data=FALSE)
CALDER_main(contact_mat_file_LJHL2, chr=1, bin_size=50E4, out_dir='./Chr01_LJHL2', sub_domains=TRUE, save_intermediate_data=FALSE)
CALDER_main(contact_mat_file_LJGR1, chr=1, bin_size=50E4, out_dir='./Chr01_LJGR1', sub_domains=TRUE, save_intermediate_data=FALSE)
CALDER_main(contact_mat_file_LJGR2, chr=1, bin_size=50E4, out_dir='./Chr01_LJGR2', sub_domains=TRUE, save_intermediate_data=FALSE)
CALDER_main(contact_mat_file_LJTZ1, chr=1, bin_size=50E4, out_dir='./Chr01_LJTZ1', sub_domains=TRUE, save_intermediate_data=FALSE)
CALDER_main(contact_mat_file_LJTZ2, chr=1, bin_size=50E4, out_dir='./Chr01_LJTZ2', sub_domains=TRUE, save_intermediate_data=FALSE)
CALDER_main(contact_mat_file_LJYP1, chr=1, bin_size=50E4, out_dir='./Chr01_LJYP1', sub_domains=TRUE, save_intermediate_data=FALSE)
CALDER_main(contact_mat_file_LJYP2, chr=1, bin_size=50E4, out_dir='./Chr01_LJYP2', sub_domains=TRUE, save_intermediate_data=FALSE)
