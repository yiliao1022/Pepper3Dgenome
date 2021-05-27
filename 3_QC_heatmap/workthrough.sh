
# Calculate the map resolution (Using 10k matrix as an example

Step1:
hicConvertFormat -m ../LJGRCF_mapping/matrix/hic_matrix_10000_3.53.h5 --outFileName LJGR2.10k.ginteractions --inputFormat h5 --outputFormat ginteractions
hicConvertFormat -m ../LJGR_mapping/matrix/hic_matrix_10000_3.53.h5 --outFileName LJGR1.10k.ginteractions --inputFormat h5 --outputFormat ginteractions
hicConvertFormat -m ../LJHL_mapping/matrix/hic_matrix_10000_3.53.h5 --outFileName LJHL1.10k.ginteractions --inputFormat h5 --outputFormat ginteractions
hicConvertFormat -m ../LYHL_mapping/matrix/hic_matrix_10000_3.53.h5 --outFileName LJHL2.10k.ginteractions --inputFormat h5 --outputFormat ginteractions
hicConvertFormat -m ../LJTZ_mapping/matrix/hic_matrix_10000_3.53.h5 --outFileName LJTZ1.10k.ginteractions --inputFormat h5 --outputFormat ginteractions 
hicConvertFormat -m ../LJTZCF_mapping/matrix/hic_matrix_10000_3.53.h5 --outFileName LJTZ2.10k.ginteractions --inputFormat h5 --outputFormat ginteractions
hicConvertFormat -m ../LJYPCF_mapping/matrix/hic_matrix_10000_3.53.h5 --outFileName LJYP2.10k.ginteractions --inputFormat h5 --outputFormat ginteractions
hicConvertFormat -m ../LYJP_mapping/matrix/hic_matrix_10000_3.53.h5 --outFileName LJYP1.10k.ginteractions --inputFormat h5 --outputFormat ginteractions

Step2:
perl Getresolution.pl LJGR2.10k.ginteractions ### This script will give how much percent the loci have at least 1000 contacts
