# Set local variables
REFS = data/references
FIGS = results/figures
TABLES = data/process/tables
PROC = data/process
FINAL = submission
CODE = code
METADATA = data/raw/metadata


# utility function to print various variables. For example, running the
# following at the command line:
#
#	make print-BAM
#
# will generate:
#	BAM=data/raw_june/V1V3_0001.bam data/raw_june/V1V3_0002.bam ...
print-%:
	@echo '$*=$($*)'

#################################################################################
#																				#
# Part 1: Get the references 													#
#																				#
# We will need several reference files to complete the analyses including the   #
# SILVA reference alignment and RDP reference taxonomy.                         #
#																				#
#################################################################################

# We want the latest greatest reference alignment and the SILVA reference
# alignment is the best reference alignment on the market. This version is from
# v123 and described at http://blog.mothur.org/2015/12/03/SILVA-v123-reference-files/
# We will use the SEED v. 123, which contain 12,083 bacterial sequences. This
# also contains the reference taxonomy. We will limit the databases to only
# include bacterial sequences.

$(REFS)/silva.seed.align :
	wget -N http://mothur.org/w/images/1/15/Silva.seed_v123.tgz
	tar xvzf Silva.seed_v123.tgz silva.seed_v123.align silva.seed_v123.tax
	mothur "#get.lineage(fasta=silva.seed_v123.align, taxonomy=silva.seed_v123.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v123.pick.align, processors=8)"
	mv silva.seed_v123.pick.align $(REFS)/silva.seed.align
	rm Silva.seed_v123.tgz silva.seed_v123.*

$(REFS)/silva.v4.align : $(REFS)/silva.seed.align
	mothur "#pcr.seqs(fasta=$(REFS)/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
	mv $(REFS)/silva.seed.pcr.align $(REFS)/silva.v4.align

# Next, we want the RDP reference taxonomy. The current version is v10 and we
# use a "special" pds version of the database files, which are described at
# http://blog.mothur.org/2014/10/28/RDP-v10-reference-files/

$(REFS)/trainset14_032015.% :
	wget -N http://mothur.org/w/images/8/88/Trainset14_032015.pds.tgz
	tar xvzf Trainset14_032015.pds.tgz trainset14_032015.pds/trainset14_032015.pds.*
	mv trainset14_032015.pds/* $(REFS)/
	rmdir trainset14_032015.pds
	rm Trainset14_032015.pds.tgz

##################################################################################
#																				 #
# Part 2: Run data through mothur 												 #
#																				 #
#	Process fastq data through the generation of files that will be used in the  #
# overall analysis.																 #
#																				 #
##################################################################################

# This runs the download of needed fastq files and runs sequence processing
# using the mothur program.  It runs up to the cluster.split step.  The
# majority of the processing can be found on www.mothur.org under the MiSeq
# SOP.

$(PROC)/unmatched.% :
	bash $(CODE)/mothur.batch

# This runs the final clustering and generation of OTUs.  A de novo clustering
# using the average neighbor algorithm was used to create OTUs.  The alpha and
# beta diversity was calculated as well. A subsampling fo the shared file is
# also performed within this call.

$(PROC)/final.% :
	bash $(CODE)/mothurCluster.batch


#################################################################################
#																				#
# Metadata Processing and General Analysis 										#
#																				#
#																				#
#################################################################################

# This modifies the meta data files by adding necessary categories (e.g. lesion)
# for files that will be used for all downstream analysis.

$(PROC)/mod_metadata/metaI_final.csv\
$(PROC)/mod_metadata/metaF_final.csv\
$(PROC)/mod_metadata/good_metaf_final.csv : $(PROC)/metaI.txt\
$(PROC)/metaF.txt $(PROC)/followup_samples.csv\
code/make_metadata_tables.R
	R -e "source('code/make_metadata_tables.R')"


# This analyzes and compares all alpha diversity metrics for lesion, adenoma,
# and carcinoma for initial and follow up samples.

$(TABLES)/alpha_table_summary.csv : $(PROC)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.groups.ave-std.summary code/Run_Alpha_Diversity_tests.R
	R -e "source('code/Run_Alpha_Diversity_tests.R')"

# This code runs the comparison of initial and follow up for the
# adenoma and carcinoma with respect to FIT and thetayc distances.

$(TABLES)/change_theta_fit_summary.csv\
$(TABLES)/difference_table.csv : $(PROC)/final.thetayc.0.03.lt.ave.dist\
$(PROC)/mod_metadata/metaF_final.csv code/Run_change_theta_Fit.R
	R -e "source('code/Run_change_theta_Fit.R')"

# This code creates the NMDS needed and also runs the PERMANOVA
# analysis for initial and follow up for either adenoma or
# carcinoma.

$(TABLES)/beta_diver_summary.csv\
$(TABLES)/thetayc_adn_IF.csv\
$(TABLES)/thetayc_srn_IF.csv\
$(TABLES)/thetayc_crc_IF.csv : $(PROC)/final.thetayc.0.03.lt.ave.dist\
$(PROC)/mod_metadata/metaF_final.csv code/Run_Beta_Diversity_tests.R
	R -e "source('code/Run_Beta_Diversity_tests.R')"


# This code runs comparisons checking for differences in time between
# initial and follow up samples for adenoma or carcinoma.

$(TABLES)/time_pvalues.csv\
$(TABLES)/time_summary_data.csv\
$(TABLES)/time_datatable.csv : $(PROC)/final.thetayc.0.03.lt.ave.dist\
$(PROC)/mod_metadata/metaI_final.csv $(PROC)/mod_metadata/metaF_final.csv\
$(PROC)/mod_metadata/good_metaf_final.csv code/Run_Supplemental_time_table.R
	R -e "source('code/Run_Supplemental_time_table.R')"


####################################################################################
#																				   #
# Model building  Adenoma Treatment												   #
#																				   #
#																				   #
####################################################################################

# Set up and run treatment models
exploratory/adn_treatment_model.RData : $(PROC)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/Adns_treatment_model.R
	R -e "source('code/Adns_treatment_model.R')"

$(TABLES)/adn_treatment_test_tune_data.csv\
$(TABLES)/adn_treatment_ROC_model_summary.csv\
$(TABLES)/adn_treatment_test_data_roc.csv\
$(TABLES)/adn_treatment_raw_mda_values.csv\
$(TABLES)/adn_treatment_MDA_Summary.csv : exploratory/adn_treatment_model.RData\
code/Adns_combine_agg_treat_model.R
	R -e "source('code/Adns_combine_agg_treat_model.R')"



####################################################################################
#																				   #
# Model building  Adenoma Treatment Random Labels								   #
#																				   #
#																				   #
####################################################################################

# Set up and run treatment models
exploratory/adn_randomized_treatment_model.RData : $(PROC)/mod_metadata/metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/Adns_randomization_treatment_model.R
	R -e "source('code/Adns_randomization_treatment_model.R')"

$(TABLES)/adn_randomization_treatment_test_tune_data.csv\
$(TABLES)/adn_randomization_treatment_ROC_model_summary.csv\
$(TABLES)/adn_randomization_treatment_test_data_roc.csv\
$(TABLES)/adn_randomization_treatment_imp_vars_summary.csv\
$(TABLES)/adn_randomization_treatment_top_vars_MDA_Summary.csv\
$(TABLES)/adn_randomization_treatment_top_vars_MDA_full_data.csv : exploratory/adn_randomized_treatment_model.RData\
source('code/Adns_randomization_combine_agg_treat_model.R
	R -e "source('code/Adns_randomization_combine_agg_treat_model.R')"



####################################################################################
#																				   #
# Model building  Adenoma Towards Normal										   #
#																				   #
#																				   #
####################################################################################


# Sets up the target file names
MODEL_NUMBER=$(shell seq 1 100)
ADN_TITLE=$(addprefix exploratory/adn_RF_model_,$(MODEL_NUMBER))
ADN_MODELS=$(addsuffix .RData,$(ADN_TITLE))

$(ADN_MODELS) : $(PROC)/final.0.03.subsample.shared\
$(PROC)/mod_metadata/metaI_final.csv $(PROC)/mod_metadata/good_metaf_final.csv\
code/adn_reference_run_RF.R code/adn_RF_reference.pbs code/setup_adn_RF_test.R\
$(CODE)/adn_createDuplicates.sh $(CODE)/adn_create_pbs.sh $(CODE)/adn_qsubmission.sh
	mkdir $(CODE)/adn
	R -e "source('code/setup_adn_RF_test.R')"
	bash $(CODE)/adn_createDuplicates.sh
	bash $(CODE)/adn_create_pbs.sh
	bash $(CODE)/adn_qsubmission.sh

# This code combines all the runs from the adenoma models and aggregates them
# together.  I also collects relevant information (e.g. AUCs).
# It also grabs the most important OTUs based on MDA and
# frequency they've occured in the 100 different runs.

$(TABLES)/adn_test_data_splits.csv\
$(TABLES)/adn_test_tune_data.csv\
$(TABLES)/adn_ROC_model_summary.csv\
$(TABLES)/adn_test_data_roc.csv\
$(TABLES)/adn_raw_mda_values.csv\
$(TABLES)/adn_MDA_Summary.csv : code/Run_adn_Combine_Testing_pull_imp_OTUs.R
	R -e "source('code/Run_adn_Combine_Testing_pull_imp_OTUs.R')"

# This code creates a 100 different 80/20 splits but with only the most
# important OTUs.  Each of the reduced adenoma models are stored as .RData
# files in the exploratory directory.

# Sets up the target file names
RED_ADN_TITLE=$(addprefix exploratory/adn_Reducedfeatures_RF_model_,$(MODEL_NUMBER))
RED_ADN_MODELS=$(addsuffix .RData,$(RED_ADN_TITLE))

$(RED_ADN_MODELS) : $(TABLES)/adn_full_test_data.csv\
$(TABLES)/adn_rf_wCV_imp_vars_summary.csv code/adn_RF_reduced_vars_reference.pbs\
code/adn_reference_run_reduced_feature_RF.R code/Run_adn_reduce_feature_lesion_model.R\
$(CODE)/adn_createDuplicates_reducedVars.sh $(CODE)/adn_create_reducedVars_pbs.sh\
$(CODE)/adn_qsubmission_reducedVars.sh
	mkdir $(CODE)/reduced_adn
	R -e "source('code/Run_adn_reduce_feature_lesion_model.R')"
	bash $(CODE)/adn_createDuplicates_reducedVars.sh
	bash $(CODE)/adn_create_reducedVars_pbs.sh
	bash $(CODE)/adn_qsubmission_reducedVars.sh

# This code gathers all the data together from the 100 different reduced adenoma model runs.
# It also stores the MDA infomration for the OTUs used in this reduced model.

$(TABLES)/adn_reduced_test_data_splits.csv\
$(TABLES)/adn_reduced_test_tune_data.csv\
$(TABLES)/adn_Reduced_ROC_model_summary.csv\
$(TABLES)/adn_reduced_test_data_roc.csv\
$(TABLES)/adn_reduced_auc_summary.csv\
$(TABLES)/adn_reduced_model_top_vars_MDA_Summary.csv\
$(TABLES)/adn_reduced_lesion_model_top_vars_MDA.csv : code/Run_combine_aggregate_reduced_adn_model.R
	R -e "source('code/Run_combine_aggregate_reduced_adn_model.R')"

# This code uses the entire normal vs adenoma cohort to generate the best model for the
# reduced ladenoma model.

$(TABLES)/adn_reduced_test_data_roc.csv\
$(TABLES)/adn_reduced_follow_up_probability_summary.csv : $(TABLES)/adn_reduced_test_tune_data.csv\
$(TABLES)/adn_Reduced_ROC_model_summary.csv $(TABLES)/adn_reduced_test_data_roc.csv\
$(TABLES)/adn_reduced_auc_summary.csv $(PROC)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/Run_adn_reduced_best_model.R
	R -e "source('code/Run_adn_reduced_best_model.R')"


####################################################################################
#																				   #
# Model building SRN Treatment   												   #
#																				   #
#																				   #
####################################################################################

# Set up and run treatment models
exploratory/srn_treatment_model.RData : $(PROC)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/SRN_treatment_model.R
	R -e "source('code/SRN_treatment_model.R')"

$(TABLES)/srn_treatment_test_tune_data.csv\
$(TABLES)/srn_treatment_ROC_model_summary.csv\
$(TABLES)/srn_treatment_test_data_roc.csv\
$(TABLES)/srn_treatment_raw_mda_values.csv\
$(TABLES)/srn_treatment_MDA_Summary.csv : exploratory/srn_treatment_model.RData\
code/SRN_combine_agg_treat_model.R
	R -e "source('code/SRN_combine_agg_treat_model.R')"



####################################################################################
#																				   #
# Model building SRN Treatment Random Labels	 								   #
#																				   #
#																				   #
####################################################################################

# Set up and run treatment models
exploratory/srn_randomized_treatment_model.RData : $(PROC)/mod_metadata/metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/SRN_randomization_treatment_model.R
	R -e "source('code/SRN_randomization_treatment_model.R')"

$(TABLES)/srn_randomization_treatment_test_tune_data.csv\
$(TABLES)/srn_randomization_treatment_ROC_model_summary.csv\
$(TABLES)/srn_randomization_treatment_test_data_roc.csv\
$(TABLES)/srn_randomization_treatment_imp_vars_summary.csv\
$(TABLES)/srn_randomization_treatment_top_vars_MDA_Summary.csv\
$(TABLES)/srn_randomization_treatment_top_vars_MDA_full_data.csv : exploratory/srn_randomized_treatment_model.RData\
source('code/SRN_randomization_combine_agg_treat_model.R
	R -e "source('code/SRN_randomization_combine_agg_treat_model.R')"



######################################################################################
#																					 #
# Model building  SRN towards normal 												 #
#																					 #
#																					 #
######################################################################################

# Sets up the target file names
SRN_TITLE=$(addprefix exploratory/srn_RF_model_,$(MODEL_NUMBER))
SRN_MODELS=$(addsuffix .RData,$(SRN_TITLE))

$(SRN_MODELS) : $(PROC)/final.0.03.subsample.shared\
$(PROC)/mod_metadata/metaI_final.csv $(PROC)/mod_metadata/good_metaf_final.csv\
code/srn_reference_run_RF.R code/srn_RF_reference.pbs code/setup_srn_RF_test.R\
$(CODE)/srn_createDuplicates.sh $(CODE)/srn_create_pbs.sh $(CODE)/srn_qsubmission.sh
	mkdir $(CODE)/srn
	R -e "source('code/setup_srn_RF_test.R')"
	bash $(CODE)/srn_createDuplicates.sh
	bash $(CODE)/srn_create_pbs.sh
	bash $(CODE)/srn_qsubmission.sh

# This code combines all the runs from the advanced adenoma models and aggregates them
# together.  I also collects relevant information (e.g. AUCs).
# It also grabs the most important OTUs based on MDA and
# frequency they've occured in the 100 different runs.

exploratory/srn_rocs.RData : code/Run_srn_Combine_Testing_pull_imp_OTUs.R
	R -e "source('code/Run_srn_Combine_Testing_pull_imp_OTUs.R')"

# This code creates a 100 different 80/20 splits but with only the most
# important OTUs.  Each of the reduced advanced adenoma models are stored as .RData
# files in the exploratory directory.

# Sets up the target file names
RED_SRN_TITLE=$(addprefix exploratory/srn_Reducedfeatures_RF_model_,$(MODEL_NUMBER))
RED_SRN_MODELS=$(addsuffix .RData,$(RED_SRN_TITLE))

$(RED_SRN_MODELS) : $(TABLES)/srn_full_test_data.csv\
$(TABLES)/srn_rf_wCV_imp_vars_summary.csv code/srn_RF_reduced_vars_reference.pbs\
code/srn_reference_run_reduced_feature_RF.R code/Run_srn_reduce_feature_lesion_model.R\
$(CODE)/srn_createDuplicates_reducedVars.sh $(CODE)/srn_create_reducedVars_pbs.sh\
$(CODE)/srn_qsubmission_reducedVars.sh
	mkdir $(CODE)/reduced_srn
	R -e "source('code/Run_srn_reduce_feature_lesion_model.R')"
	bash $(CODE)/srn_createDuplicates_reducedVars.sh
	bash $(CODE)/srn_create_reducedVars_pbs.sh
	bash $(CODE)/srn_qsubmission_reducedVars.sh

# This code gathers all the data together from the 100 different reduced advanced adenoma
# model runs.  It also stores the MDA infomration for the OTUs used in this reduced model.

$(TABLES)/srn_reduced_test_data_splits.csv\
$(TABLES)/srn_reduced_test_tune_data.csv\
$(TABLES)/srn_Reduced_ROC_model_summary.csv\
$(TABLES)/srn_reduced_test_data_roc.csv\
$(TABLES)/srn_reduced_auc_summary.csv\
$(TABLES)/srn_reduced_model_top_vars_MDA_Summary.csv\
$(TABLES)/srn_reduced_lesion_model_top_vars_MDA.csv : code/Run_combine_aggregate_reduced_srn_model.R
	#Collects the needed data to generate figure 3
	R -e "source('code/Run_combine_aggregate_reduced_srn_model.R')"

# This code uses the entire normal vs advanced adenoma cohort to generate the best model
# for the reduced advanced adenoma model.

$(TABLES)/srn_reduced_test_data_roc.csv\
$(TABLES)/srn_reduced_follow_up_probability_summary.csv : $(TABLES)/srn_reduced_test_tune_data.csv\
$(TABLES)/srn_Reduced_ROC_model_summary.csv $(TABLES)/srn_reduced_test_data_roc.csv\
$(TABLES)/srn_reduced_auc_summary.csv $(PROC)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/Run_srn_reduced_best_model.R
	R -e "source('code/Run_srn_reduced_best_model.R')"


####################################################################################
#																				   #
# Model building CRC Treatment   												   #
#																				   #
#																				   #
####################################################################################

# Set up and run treatment models
exploratory/crc_treatment_model.RData : $(PROC)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/CRC_treatment_model.R
	R -e "source('code/CRC_treatment_model.R')"

$(TABLES)/crc_treatment_test_tune_data.csv\
$(TABLES)/crc_treatment_ROC_model_summary.csv\
$(TABLES)/crc_treatment_test_data_roc.csv\
$(TABLES)/crc_treatment_raw_mda_values.csv\
$(TABLES)/crc_treatment_MDA_Summary.csv : exploratory/crc_treatment_model.RData\
code/CRC_combine_agg_treat_model.R
	R -e "source('code/CRC_combine_agg_treat_model.R')"



####################################################################################
#																				   #
# Model building CRC Treatment Random Labels	 								   #
#																				   #
#																				   #
####################################################################################

# Set up and run treatment models
exploratory/crc_randomized_treatment_model.RData : $(PROC)/mod_metadata/metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/CRC_randomization_treatment_model.R
	R -e "source('code/CRC_randomization_treatment_model.R')"

$(TABLES)/crc_randomization_treatment_test_tune_data.csv\
$(TABLES)/crc_randomization_treatment_ROC_model_summary.csv\
$(TABLES)/crc_randomization_treatment_test_data_roc.csv\
$(TABLES)/crc_randomization_treatment_imp_vars_summary.csv\
$(TABLES)/crc_randomization_treatment_top_vars_MDA_Summary.csv\
$(TABLES)/crc_randomization_treatment_top_vars_MDA_full_data.csv : exploratory/crc_randomized_treatment_model.RData\
source('code/CRC_randomization_combine_agg_treat_model.R
	R -e "source('code/CRC_randomization_combine_agg_treat_model.R')"


###################################################################################
#																			 	  #
# Model building  CRC 														 	  #
#																			 	  #
#																			 	  #
###################################################################################

# Sets up the target file names
CRC_TITLE=$(addprefix exploratory/crc_RF_model_,$(MODEL_NUMBER))
CRC_MODELS=$(addsuffix .RData,$(CRC_TITLE))

$(CRC_MODELS) : $(PROC)/final.0.03.subsample.shared\
$(PROC)/mod_metadata/metaI_final.csv $(PROC)/mod_metadata/good_metaf_final.csv\
code/crc_reference_run_RF.R code/crc_RF_reference.pbs code/setup_crc_RF_test.R\
$(CODE)/crc_createDuplicates.sh $(CODE)/crc_create_pbs.sh $(CODE)/crc_qsubmission.sh
	mkdir $(CODE)/crc
	R -e "source('code/setup_crc_RF_test.R')"
	bash $(CODE)/crc_createDuplicates.sh
	bash $(CODE)/crc_create_pbs.sh
	bash $(CODE)/crc_qsubmission.sh

# This code combines all the runs from the carcinoma models and aggregates them
# together.  I also collects relevant information (e.g. AUCs).
# It also grabs the most important OTUs based on MDA and
# frequency they've occured in the 100 different runs.

exploratory/crc_rocs.RData : code/Run_Combine_Testing_pull_imp_OTUs.R
	R -e "source('code/Run_crc_Combine_Testing_pull_imp_OTUs.R')"

# This code creates a 100 different 80/20 splits but with only the most
# important OTUs.  Each of the reduced lcarcinomamodels are stored as .RData
# files in the exploratory directory.

# Sets up the target file names
RED_CRC_TITLE=$(addprefix exploratory/crc_Reducedfeatures_RF_model_,$(MODEL_NUMBER))
RED_CRC_MODELS=$(addsuffix .RData,$(RED_CRC_TITLE))

$(RED_CRC_MODELS) : $(TABLES)/crc_full_test_data.csv\
$(TABLES)/crc_rf_wCV_imp_vars_summary.csv code/crc_RF_reduced_vars_reference.pbs\
code/crc_reference_run_reduced_feature_RF.R code/Run_crc_reduce_feature_lesion_model.R\
$(CODE)/crc_createDuplicates_reducedVars.sh $(CODE)/crc_create_reducedVars_pbs.sh\
$(CODE)/crc_qsubmission_reducedVars.sh
	mkdir $(CODE)/reduced_crc
	R -e "source('code/Run_crc_reduce_feature_lesion_model.R')"
	bash $(CODE)/crc_createDuplicates_reducedVars.sh
	bash $(CODE)/crc_create_reducedVars_pbs.sh
	bash $(CODE)/crc_qsubmission_reducedVars.sh


# This code gathers all the data together from the 100 different reduced carcinoma
# model runs. It also stores the MDA infomration for the OTUs used in this reduced model.

$(TABLES)/crc_reduced_test_data_splits.csv\
$(TABLES)/crc_reduced_test_tune_data.csv\
$(TABLES)/crc_Reduced_ROC_model_summary.csv\
$(TABLES)/crc_reduced_test_data_roc.csv\
$(TABLES)/crc_reduced_auc_summary.csv\
$(TABLES)/crc_reduced_model_top_vars_MDA_Summary.csv\
$(TABLES)/crc_reduced_lesion_model_top_vars_MDA.csv : code/Run_combine_aggregate_reduced_crc_model.R
	R -e "source('code/Run_combine_aggregate_reduced_crc_model.R')"

# This code uses the entire normal vs carcinoma cohort to generate the best model for the
# reduced carcinoma model.

$(TABLES)/crc_reduced_test_data_roc.csv\
$(TABLES)/crc_reduced_follow_up_probability_summary.csv : $(TABLES)/crc_reduced_test_tune_data.csv\
$(TABLES)/crc_Reduced_ROC_model_summary.csv $(TABLES)/crc_reduced_test_data_roc.csv\
$(TABLES)/crc_reduced_auc_summary.csv $(PROC)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/Run_crc_reduced_best_model.R
	R -e "source('code/Run_crc_reduced_best_model.R')"



########################################################################################
#																					   #
# Generalized Analysis & Common Comparisons 										   #
#																					   #
#																					   #
########################################################################################


# This code runs comparisons for initial versus follow up differences for
# adenoma, advanced adenoma, and carcinoma.

$(TABLES)/OTU_paired_wilcoxson_test.csv : $(TABLES)/full_test_data.csv\
$(PROC)/final.0.03.subsample.shared $(PROC)/mod_metadata/good_metaf_final.csv\
code/Run_wilcoxson_all.R $(PROC)/mod_metadata/metaI_final.csv\
	R -e "source('code/Run_wilcoxson_all.R')"


# This code runs comparisons on the positive probability for initial versus follow up
# samples for all reduced final models.

$(TABLES)/all_crc_srn_adn_models_wilcox_paired_pvalue_summary.csv\
$(TABLES)/all_crc_srn_adn_models_confusion_summary.csv\
$(TABLES)/all_crc_srn_adn_models_summary_info.csv : $(TABLES)/adn_reduced_follow_up_probability_summary.csv\
$(TABLES)/crc_reduced_follow_up_probability_summary.csv\
$(TABLES)/srn_reduced_follow_up_probability_summary.csv\
$(PROC)/mod_metadata/good_metaf_final.csv code/Run_adn_srn_crc_probs_comparison.R
	R -e "source('code/Run_adn_srn_crc_probs_comparison.R')"


# The generation and storage of the taxonomies for the OTUs used in the
# reduced adenoma, advanced adenoma, and carcinoma models.

$(TABLES)/adn_rf_otu_tax.csv\
$(TABLES)/srn_rf_otu_tax.csv\
$(TABLES)/crc_rf_otu_tax.csv : $(PROC)/final.taxonomy $(TABLES)/IF_rf_wCV_imp_vars_summary.csv\
$(TABLES)/rf_wCV_imp_vars_summary.csv code/Run_ID_imp_OTUs.R
	R -e "source('code/Run_ID_imp_OTUs.R')"


# This code IDs the common OTUs between the three different models and
# also runs a comparison for differences between initial and follow up
# samples.

$(TABLES)/pvalue_adn_srn_crc_common_imp_vars.csv : $(TABLES)/adn_rf_otu_tax.csv\
$(TABLES)/srn_rf_otu_tax.csv $(TABLES)/crc_rf_otu_tax.csv\
$(PROC)/final.0.03.subsample.shared $(PROC)/mod_metadata/good_metaf_final.csv\
code/Run_adn_crc_Compare_models.R
	R -e "source('code/Run_adn_crc_Compare_models.R')"


# The comparisons for differences in initial and follow up samples based on whether
# radiation or chemotherapy fro caricnoma and whether surgery was used for the adneoma
# group.

$(TABLES)/crc_probs_chemo_rad_pvalue_summary.csv\
$(TABLES)/adn_combined_probs_surgery_pvalue_summary.csv\
$(TABLES)/crc_chemo_rad_summary.csv\
$(TABLES)/adn_combined_surgery_summary.csv : $(TABLES)/adn_reduced_follow_up_probability_summary.csv\
$(TABLES)/srn_reduced_follow_up_probability_summary.csv\
$(TABLES)/crc_reduced_follow_up_probability_summary.csv $(TABLES)/difference_table.csv\
$(PROC)/mod_metadata/good_metaf_final.csv $(PROC)/final.groups.ave-std.summary\
code/Run_adn_crc_Test_Chemo_Rad.R
	R -e "source('code/Run_adn_crc_Test_Chemo_Rad.R')"

# The comparisons for differences in initial and follow up samples based on whether
# radiation or chemotherapy fro caricnoma and whether surgery was used for the adneoma
# group. Main diffference from above is that it looks at common OTUs only.

$(TABLES)/chemo_rads_treatment_pvalue_summary.csv\
$(TABLES)/all_adn_surg_pvalue_summary.csv : $(TABLES)/crc_probs_chemo_rad_pvalue_summary.csv\
$(PROC)/mod_metadata/good_metaf_final.csv $(TABLES)/pvalue_adn_srn_crc_common_imp_vars.csv\
$(TABLES)/crc_chemo_rad_summary.csv $(PROC)/final.shared code/Run_treatment_common_otus.R
	R -e "source('code/Run_treatment_common_otus.R')"



###################################################################################
#																				  #
# Part 4: Run Figures                                                             #
#																				  #
#	Run scripts to generate figures 											  #
#																				  #
###################################################################################

# This figure looks at difference between initial and follow based on
# whether the individual had an adenoma, advanced adenoma, or carcinoma.

$(FIGS)/Figure1.pdf : $(TABLES)/difference_table.csv\
$(TABLES)/change_theta_fit_summary.csv $(TABLES)/thetayc_adn_IF.csv\
$(TABLES)/thetayc_crc_IF.csv $(TABLES)/thetayc_srn_IF.csv\
$(TABLES)/beta_diver_summary.csv code/Run_Figure1.R
	R -e "source('code/Run_Figure1.R')"


# This figure explores the MDA of the 10 OTUs in the  adenoma,
# advanced adenoma, and carincoma treatment models.
$(FIGS)/Figure2.pdf : $(TABLES)/reduced_adn_treatment_top_vars_MDA_Summary.csv\
$(TABLES)/reduced_srn_treatment_top_vars_MDA_Summary.csv\
$(TABLES)/reduced_crc_treatment_top_vars_MDA_Summary.csv\
$(TABLES)/reduced_adn_treatment_top_vars_MDA_full_data.csv\
$(TABLES)/reduced_srn_treatment_top_vars_MDA_full_data.csv\
$(TABLES)/reduced_crc_treatment_top_vars_MDA_full_data.csv\
code/treatment_otu_graph.R
	R -e "source('code/treatment_otu_graph.R')"


# This figure explores how the rank importance of common OTUs between the adenoma,
# advanced adenoma, and carincoma models change depending on the model.

$(FIGS)/Figure3.pdf : $(TABLES)/adn_reduced_model_top_vars_MDA_Summary.csv\
$(TABLES)/srn_reduced_model_top_vars_MDA_Summary.csv\
$(TABLES)/crc_reduced_model_top_vars_MDA_Summary.csv\
code/common_all_models.R code/common_otu_graph.R
	R -e "source('code/common_otu_graph.R')"


# This figure looks at classification of intial and follow up samples
# based on adenoma, advanced adenoma, or carcinoma.

$(FIGS)/Figure4.pdf : $(TABLES)/adn_reduced_follow_up_probability_summary.csv\
$(TABLES)/srn_reduced_follow_up_probability_summary.csv\
$(TABLES)/crc_reduced_follow_up_probability_summary.csv\
code/Run_Figure2.R
	R -e "source('code/Run_Figure4.R')"


# This graph highlights the ROC curve for each of the models used.  This
# visualizes the models normal vs adenoma, normal vs advanced adenoma, and
# normal vs carcinoma.

$(FIGS)/FigureS1.pdf : $(TABLES)/adn_reduced_test_data_roc.csv\
$(TABLES)/srn_reduced_test_data_roc.csv\
$(TABLES)/crc_reduced_test_data_roc.csv code/Run_FigureS1.R
	R -e "source('code/Run_FigureS1.R')"
	tiff2pdf -z -o results/figures/FigureS1.pdf results/figures/FigureS1.tiff
	rm results/figures/FigureS1.tiff


# This figure explores the most important OTUs for each model and displays
# them by median MDA.

$(FIGS)/FigureS2.pdf : $(TABLES)/adn_reduced_model_top_vars_MDA_Summary.csv\
$(TABLES)/adn_reduced_lesion_model_top_vars_MDA.csv\
$(TABLES)/srn_reduced_model_top_vars_MDA_Summary.csv\
$(TABLES)/srn_reduced_model_top_vars_MDA.csv\
$(TABLES)/reduced_crc_model_top_vars_MDA_Summary.csv\
$(TABLES)/crc_reduced_lesion_model_top_vars_MDA.csv code/Run_FigureS2.R
	R -e "source('code/Run_FigureS2.R')"


# This figure plots the oral specific crc-releated OTUs identified in the carcinoma model.

$(FIGS)/FigureS3.pdf : $(PROC)/final.shared\
$(TABLES)/reduced_crc_model_top_vars_MDA_Summary.csv\
$(PROC)/mod_metadata/good_metaf_final.csv\
$(TABLES)/crc_rf_otu_tax.csv code/Run_FigureS3.R
	R -e "source('code/Run_FigureS3.R')"


#####################################################################################
#																					#
# Part 5: Pull it all together 														#
#																					#
# Render the manuscript 															#
#																					#
#####################################################################################


write.paper : $(FINAL)/manuscript.Rmd\
		$(FINAL)/supplement.Rmd\
		results/tables/Table1.Rmd results/tables/Table2.Rmd\
		$(TABLES)/mod_metadata/good_metaf_final.csv\
		$(FIGS)/Figure1.pdf $(FIGS)/Figure2.pdf\
		$(FIGS)/Figure3.pdf $(FIGS).Figure4.pdf\
		$(FIGS)/FigureS2.pdf\ $(FIGS)/FigureS3.pdf 
		$(FIGS)/FigureS4.pdf code/Run_render_paper.R
	R -e "source('code/Run_render_paper.R')"

write.revision1.paper : $(FINAL)/manuscript_R1.Rmd\
		$(FIGS)/FigureS1.pdf code/Run_render_revision1_paper.R
	R -e "source('code/Run_render_revision1_paper.R')"

write.r1.marked.up : $(FINAL)/manuscript.tex\
		$(FINAL)/manuscript_R1.tex
	latexdiff $(FINAL)/manuscript.tex $(FINAL)/manuscript_R1.tex > $(FINAL)/manuscript_R1_markedup.tex
	pdflatex -output-directory=$(FINAL) $(FINAL)/manuscript_R1_markedup.tex
	rm $(FINAL)/manuscript_R1_markedup.aux 
	rm $(FINAL)/manuscript_R1_markedup.log
	rm $(FINAL)/manuscript_R1_markedup.out}
