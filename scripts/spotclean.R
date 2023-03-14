# 2.12.2022
# Author: Antti Kiviaho
#
# Script for running spotclean on a set of prostate ST samples

library(SpotClean)

data_path = '/lustre/scratch/kiviaho/prostate_spatial/results'
res_path = '/lustre/scratch/kiviaho/prostate_spatial/results/after_spotclean'


samples = c('BPH_651','BPH_665','BPH_688',
'CRPC-278','CRPC-489','CRPC-697','CRPC_531',
'PC-03-6712','PC_00_16338_II','PC_01_06342_VAS',
'PC_01_14451_OIK','PC_02_05601_OIK','PC_02_10136_VAS',
'PC_03_01669_TUTKV','PC_15420OIK','PC_4980','PC_7875OIK')

for(i in 1:length(samples)){
	sample = samples[i]
	print(sample)
	# load count matrix and slide metadata
	data_dir <- file.path(data_path,sample,'outs')
	raw_data = read10xRawH5(h5_file = file.path(data_dir,'raw_feature_bc_matrix.h5'))


	spatial_dir <- file.path(data_path,sample,'outs/spatial')

	slide_info <- read10xSlide(tissue_csv_file=file.path(spatial_dir,
"tissue_positions_list.csv"),
	tissue_img_file = file.path(spatial_dir,
"tissue_lowres_image.png"),
	scale_factor_file = file.path(spatial_dir,
"scalefactors_json.json"))

	# Create slide object
	slide_obj <- createSlide(raw_data,slide_info)

	spotclean_res = spotclean(slide_obj)
	spotclean_res$metadata
	spotclean_res
	cleaned_data = SummarizedExperiment::assays(spotclean_res)$decont

	# Save data into a matrix
	write.table(cleaned_data,file=file.path(res_path,paste0(sample,'_counts_after_spotclean.csv')),sep=',')

	print('')
	print('##################################')
	print('')
}

