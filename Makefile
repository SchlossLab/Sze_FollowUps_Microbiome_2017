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

################################################################################
#
# Part 1: Get the references
#
# We will need several reference files to complete the analyses including the
# SILVA reference alignment and RDP reference taxonomy.
#
################################################################################

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

################################################################################
#
# Part 2: Run data through mothur
#
#	Process fastq data through the generation of files that will be used in the
# overall analysis.
#
################################################################################

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

################################################################################
#
# Part 3: Run Analysis
#
#	Run scripts to generate tables needed for figures, tables, and manuscript
#
################################################################################

# This modifies the meta data files by adding necessary categories (e.g. lesion)
# for files that will be used for all downstream analysis.

$(TABLES)/mod_metadata/good_metaf_final.csv : $(METADATA)/followUps_metadata.txt\
$(METADATA)/initials_metadata.tsv $(METADATA)/followUp_outcome_data.csv\
code/make_metadata_tables.R 
	R -e "source('code/make_metadata_tables.R')"

# This analyzes and compares all alpha diversity metrics for lesion, adenoma, 
# and carcinoma for initial and follow up samples.

$(TABLES)/alpha_table_summary.csv : $(PROC)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.groups.ave-std.summary code/Run_Alpha_Diversity_tests.R
	R -e "source('code/Run_Alpha_Diversity_tests.R')"

# This code runs the comparison of initial and follow up for the 
# adenoma and carcinoma with respect to FIT and thetayc distances.

$(TABLES)/difference_table.csv : $(PROC)/final.thetayc.0.03.lt.ave.dist\
$(PROC)/mod_metadata/metaF_final.csv code/Run_change_theta_Fit.R
	R -e "source('code/Run_change_theta_Fit.R')"

# This code creates the NMDS needed and also runs the PERMANOVA
# analysis for initial and follow up for either adenoma or 
# carcinoma.

$(TABLES)/thetayc_% : $(PROC)/final.thetayc.0.03.lt.ave.dist\
$(PROC)/mod_metadata/metaF_final.csv code/Run_Beta_Diversity_tests.R
	R -e "source('code/Run_Beta_Diversity_tests.R')"

# This code grabs the previously CRC associated bacteria and runs
# a comparison for initial and follow up for either adenoma or
# carcinoma.  

$(TABLES)/adn_crc_maybe_pvalue_summary.csv : $(PROC)/final.taxonomy\
$(PROC)/mod_metadata/good_metaf_final.csv $(PROC)/final.shared\
code/Run_potential_cancer_specific_OTUs.R
	R -e "source('code/Run_potential_cancer_specific_OTUs.R')"

# This code runs comparisons checking for differences in time between
# initial and follow up samples for adenoma or carcinoma.

$(TABLES)/time_pvalues.csv : $(PROC)/final.thetayc.0.03.lt.ave.dist\
$(PROC)/mod_metadata/metaI_final.csv $(PROC)/mod_metadata/metaF_final.csv\
$(PROC)/mod_metadata/good_metaf_final.csv code/Run_Supplemental_time_table.R
	R -e "source('code/Run_Supplemental_time_table.R')"



########## Model building for normal versus carcinoma ##########################
exploratory/crc_RF_model_100.RData : $(PROC)/final.0.03.subsample.shared\
$(PROC)/mod_metadata/metaI_final.csv $(PROC)/mod_metadata/good_metaf_final.csv\
code/crc_reference_run_RF.R code/crc_RF_reference.pbs code/setup_crc_RF_test.R\
$(CODE)/crc_createDuplicates.sh $(CODE)/crc_create_pbs.sh $(CODE)/crc_qsubmission.sh
	mkdir $(CODE)/crc
	R -e "source('code/setup_crc_RF_test.R')"
	bash $(CODE)/crc_createDuplicates.sh
	bash $(CODE)/crc_create_pbs.sh
	bash $(CODE)/crc_qsubmission.sh

# This code combines all the runs from the lesion models and aggregates them
# together.  I also collects relevant information (e.g. AUCs).  
# It also grabs the most important OTUs based on MDA and 
# frequency they've occured in the 100 different runs.

exploratory/crc_rocs.RData : code/Run_Combine_Testing_pull_imp_OTUs.R
	R -e "source('code/Run_crc_Combine_Testing_pull_imp_OTUs.R')"

# This code creates a 100 different 80/20 splits but with only the most
# important OTUs.  Each of the reduced lesion models are stored as .RData
# files in the exploratory directory.

exploratory/crc_Reducedfeatures_RF_model_100.RData : $(TABLES)/crc_full_test_data.csv\
$(TABLES)/crc_rf_wCV_imp_vars_summary.csv code/crc_RF_reduced_vars_reference.pbs\
code/crc_reference_run_reduced_feature_RF.R code/Run_crc_reduce_feature_lesion_model.R\
$(CODE)/crc_createDuplicates_reducedVars.sh $(CODE)/crc_create_reducedVars_pbs.sh\
$(CODE)/crc_qsubmission_reducedVars.sh
	mkdir $(CODE)/reduced_crc
	R -e "source('code/Run_crc_reduce_feature_lesion_model.R')"
	bash $(CODE)/crc_createDuplicates_reducedVars.sh
	bash $(CODE)/crc_create_reducedVars_pbs.sh
	bash $(CODE)/crc_qsubmission_reducedVars.sh


# This code gathers all the data together from the 100 different reduced lesion model runs.
# It also stores the MDA infomration for the OTUs used in this reduced model.

$(TABLES)/reduced_crc_model_top_vars_MDA_Summary.csv : code/Run_combine_aggregate_reduced_crc_model.R
	#Collects the needed data to generate figure 3
	R -e "source('code/Run_combine_aggregate_reduced_crc_model.R')"

# This code uses the entire 423-person cohort to generate the best model for the 
# reduced lesion model.

$(TABLES)/reduced_crc_follow_up_probability_summary.csv : $(TABLES)/crc_reduced_test_tune_data.csv\
$(TABLES)/crc_Reduced_ROC_model_summary.csv $(TABLES)/crc_reduced_test_data_roc.csv\
$(TABLES)/crc_reduced_auc_summary.csv $(PROC)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/Run_crc_reduced_best_model.R
	R -e "source('code/Run_crc_reduced_best_model.R')"


########## Model building for normal versus adenoma
exploratory/adn_RF_model_100.RData : $(PROC)/final.0.03.subsample.shared\
$(PROC)/mod_metadata/metaI_final.csv $(PROC)/mod_metadata/good_metaf_final.csv\
code/adn_reference_run_RF.R code/adn_RF_reference.pbs code/setup_adn_RF_test.R\
$(CODE)/adn_createDuplicates.sh $(CODE)/adn_create_pbs.sh $(CODE)/adn_qsubmission.sh
	mkdir $(CODE)/adn
	R -e "source('code/setup_adn_RF_test.R')"
	bash $(CODE)/adn_createDuplicates.sh
	bash $(CODE)/adn_create_pbs.sh
	bash $(CODE)/adn_qsubmission.sh

# This code combines all the runs from the lesion models and aggregates them
# together.  I also collects relevant information (e.g. AUCs).  
# It also grabs the most important OTUs based on MDA and 
# frequency they've occured in the 100 different runs.

exploratory/adn_rocs.RData : code/Run_adn_Combine_Testing_pull_imp_OTUs.R
	R -e "source('code/Run_adn_Combine_Testing_pull_imp_OTUs.R')"

# This code creates a 100 different 80/20 splits but with only the most
# important OTUs.  Each of the reduced lesion models are stored as .RData
# files in the exploratory directory.

exploratory/adn_Reducedfeatures_RF_model_100.RData : $(TABLES)/adn_full_test_data.csv\
$(TABLES)/adn_rf_wCV_imp_vars_summary.csv code/adn_RF_reduced_vars_reference.pbs\
code/adn_reference_run_reduced_feature_RF.R code/Run_adn_reduce_feature_lesion_model.R\
$(CODE)/adn_createDuplicates_reducedVars.sh $(CODE)/adn_create_reducedVars_pbs.sh\
$(CODE)/adn_qsubmission_reducedVars.sh
	mkdir $(CODE)/reduced_adn
	R -e "source('code/Run_adn_reduce_feature_lesion_model.R')"
	bash $(CODE)/adn_createDuplicates_reducedVars.sh
	bash $(CODE)/adn_create_reducedVars_pbs.sh
	bash $(CODE)/adn_qsubmission_reducedVars.sh

# This code gathers all the data together from the 100 different reduced lesion model runs.
# It also stores the MDA infomration for the OTUs used in this reduced model.

$(TABLES)/reduced_adn_model_top_vars_MDA_Summary.csv : code/Run_combine_aggregate_reduced_adn_model.R
	#Collects the needed data to generate figure 3
	R -e "source('code/Run_combine_aggregate_reduced_adn_model.R')"

# This code uses the entire 423-person cohort to generate the best model for the 
# reduced lesion model.

$(TABLES)/reduced_adn_follow_up_probability_summary.csv : $(TABLES)/adn_reduced_test_tune_data.csv\
$(TABLES)/adn_Reduced_ROC_model_summary.csv $(TABLES)/adn_reduced_test_data_roc.csv\
$(TABLES)/adn_reduced_auc_summary.csv $(PROC)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/Run_adn_reduced_best_model.R
	R -e "source('code/Run_adn_reduced_best_model.R')"

########## Model building for normal versus SRN
exploratory/srn_RF_model_100.RData : $(PROC)/final.0.03.subsample.shared\
$(PROC)/mod_metadata/metaI_final.csv $(PROC)/mod_metadata/good_metaf_final.csv\
code/srn_reference_run_RF.R code/srn_RF_reference.pbs code/setup_srn_RF_test.R\
$(CODE)/srn_createDuplicates.sh $(CODE)/srn_create_pbs.sh $(CODE)/srn_qsubmission.sh
	mkdir $(CODE)/srn
	R -e "source('code/setup_srn_RF_test.R')"
	bash $(CODE)/srn_createDuplicates.sh
	bash $(CODE)/srn_create_pbs.sh
	bash $(CODE)/srn_qsubmission.sh

# This code combines all the runs from the lesion models and aggregates them
# together.  I also collects relevant information (e.g. AUCs).  
# It also grabs the most important OTUs based on MDA and 
# frequency they've occured in the 100 different runs.

exploratory/srn_rocs.RData : code/Run_srn_Combine_Testing_pull_imp_OTUs.R
	R -e "source('code/Run_srn_Combine_Testing_pull_imp_OTUs.R')"

# This code creates a 100 different 80/20 splits but with only the most
# important OTUs.  Each of the reduced lesion models are stored as .RData
# files in the exploratory directory.

exploratory/srn_Reducedfeatures_RF_model_100.RData : $(TABLES)/srn_full_test_data.csv\
$(TABLES)/srn_rf_wCV_imp_vars_summary.csv code/srn_RF_reduced_vars_reference.pbs\
code/srn_reference_run_reduced_feature_RF.R code/Run_srn_reduce_feature_lesion_model.R\
$(CODE)/srn_createDuplicates_reducedVars.sh $(CODE)/srn_create_reducedVars_pbs.sh\
$(CODE)/srn_qsubmission_reducedVars.sh
	mkdir $(CODE)/reduced_srn
	R -e "source('code/Run_srn_reduce_feature_lesion_model.R')"
	bash $(CODE)/srn_createDuplicates_reducedVars.sh
	bash $(CODE)/srn_create_reducedVars_pbs.sh
	bash $(CODE)/srn_qsubmission_reducedVars.sh

# This code gathers all the data together from the 100 different reduced lesion model runs.
# It also stores the MDA infomration for the OTUs used in this reduced model.

$(TABLES)/reduced_srn_model_top_vars_MDA_Summary.csv : code/Run_combine_aggregate_reduced_srn_model.R
	#Collects the needed data to generate figure 3
	R -e "source('code/Run_combine_aggregate_reduced_srn_model.R')"

# This code uses the entire 423-person cohort to generate the best model for the 
# reduced lesion model.

$(TABLES)/reduced_srn_follow_up_probability_summary.csv : $(TABLES)/srn_reduced_test_tune_data.csv\
$(TABLES)/srn_Reduced_ROC_model_summary.csv $(TABLES)/srn_reduced_test_data_roc.csv\
$(TABLES)/srn_reduced_auc_summary.csv $(PROC)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/Run_srn_reduced_best_model.R
	R -e "source('code/Run_srn_reduced_best_model.R')"


# This code runs comparisons for initial versus follow up differences for
# lesion, adenoma only, and carcinoma only.

$(TABLES)/OTU_paired_wilcoxson_test.csv : $(TABLES)/full_test_data.csv\
$(PROC)/final.0.03.subsample.shared $(PROC)/mod_metadata/good_metaf_final.csv\
code/Run_wilcoxson_all.R $(PROC)/mod_metadata/metaI_final.csv\
	R -e "source('code/Run_wilcoxson_all.R')"



# This code runs comparisons on the positive probability for initial versus follow up samples
# for both the reduced lesion and initial sample models.

adn.srn.crc.probs.comparison : $(TABLES)/adn_reduced_follow_up_probability_summary.csv\
$(TABLES)/crc_reduced_follow_up_probability_summary.csv\
$(TABLES)/srn_reduced_follow_up_probability_summary.csv\
$(PROC)/mod_metadata/good_metaf_final.csv code/Run_adn_srn_crc_probs_comparison.R
	R -e "source('code/Run_adn_srn_crc_probs_comparison.R')"


# The generation and storage of the taxonomies for the OTUs used in either the 
# reduced lesion or reduced initial sample model.

get.model.imp.taxa : $(PROC)/final.taxonomy $(TABLES)/IF_rf_wCV_imp_vars_summary.csv\
$(TABLES)/rf_wCV_imp_vars_summary.csv code/Run_ID_imp_OTUs.R
	R -e "source('code/Run_ID_imp_OTUs.R')"


# This code IDs the common OTUs between the two different models and
# also runs a comparison for differences between initial and follow up 
# samples.

get.common.otus : $(TABLES)/adn_rf_otu_tax.csv\
$(TABLES)/crc_rf_otu_tax.csv $(PROC)/final.0.03.subsample.shared\
$(PROC)/mod_metadata/good_metaf_final.csv code/Run_adn_crc_Compare_models.R
	R -e "source('code/Run_adn_crc_Compare_models.R')"


# The comparisons for differences in initial and follow up samples based on whether
# radiation or chemotherapy was used was completed with the below R script.

adn.crc.chemo.rads.comparison : $(TABLES)/adn_reduced_follow_up_probability_summary.csv\
$(TABLES)/crc_reduced_follow_up_probability_summary.csv $(TABLES)/difference_table.csv\
$(PROC)/mod_metadata/good_metaf_final.csv $(PROC)/final.groups.ave-std.summary\
code/Run_adn_crc_Test_Chemo_Rad.R
	R -e "source('code/Run_adn_crc_Test_Chemo_Rad.R')"


# This code runs comparisons on the positive probability for initial versus follow up samples
# for both the reduced lesion and initial sample models.

$(TABLES)/all_models_wilcox_paired_pvalue_summary.csv : $(TABLES)/follow_up_probability_summary.csv\
$(TABLES)/reduced_follow_up_probability_summary.csv $(TABLES)/IF_follow_up_probability_summary.csv\
$(TABLES)/reduced_IF_follow_up_probability_summary.csv $(PROC)/mod_metadata/good_metaf_final.csv\
code/Run_probs_comparison.R
	R -e "source('code/Run_probs_comparison.R')"



################################################################################
#
# Part 4: Run Figures
#
#	Run scripts to generate figures
#
################################################################################

# This figure looks at difference between initial and follow based on 
# whether the individual had an adenoma or carcinoma.

$(FIGS)/Figure1.pdf : $(TABLES)/difference_table.csv\
$(TABLES)/change_theta_fit_summary.csv $(TABLES)/thetayc_adn_IF.csv\
$(TABLES)/thetayc_crc_IF.csv $(TABLES)/beta_diver_summary.csv\
code/Run_Figure1.R
	R -e "source('code/Run_Figure1.R')"

# This figure looks specifically at 4 commonly CRC associated bacteria
# and whether there are diferences in initial versus follow ups based on
# adenoma or carcinoma.

$(FIGS)/Figure2.pdf : $(TABLES)/adn_crc_maybe_diff.csv code/Run_Figure2.R
	R -e "source('code/Run_Figure2.R')"

# This figure explores the reduced lesion model.  Overall model
# performance, MDA of OTUs in it, and how it performs in adenoma 
# and carcinoma for initial and follow up samples.

$(FIGS)/Figure3.pdf : $(TABLES)/reduced_test_data_roc.csv\
$(TABLES)/reduced_lesion_model_top_vars_MDA_Summary.csv\
$(TABLES)/reduced_lesion_model_top_vars_MDA.csv\
$(TABLES)/reduced_follow_up_probability_summary.csv $(PROC)/final.taxonomy code/Run_Figure3.R
	#Creates the actual Figure 3
	R -e "source('code/Run_Figure3.R')"
	tiff2pdf -z -o results/figures/Figure3.pdf results/figures/Figure3.tiff
	rm results/figures/Figure3.tiff

# This figure explores the reduced initial sample model.  Overall model
# performance, MDA of OTUs in it, and how it performs in adenoma 
# and carcinoma for initial and follow up samples.

$(FIGS)/Figure4.pdf : $(TABLES)/reduced_IF_test_data_roc.csv\
$(TABLES)/reduced_IF_model_top_vars_MDA_Summary.csv\
$(TABLES)/reduced_IF_model_top_vars_MDA.csv\
$(TABLES)/reduced_IF_follow_up_probability_summary.csv $(PROC)/final.taxonomy code/Run_Figure4.R
	R -e "source('code/Run_Figure4.R')"
	tiff2pdf -z -o results/figures/Figure4.pdf results/figures/Figure4.tiff
	rm results/figures/Figure4.tiff

# This figure summarizes the p-values found for all OTUs in comparing
# initial versus follow up for lesion, adenoma only, and carcinoma only.

$(FIGS)/FigureS1.pdf : $(TABLES)/OTU_paired_wilcoxson_test.csv\
code/Run_FigureS1.R
	R -e "source('code/Run_FigureS1.R')"

# This graph highlights how the thetayc's don't vary that much even
# though there is a marked difference in time of collection of follow up 
# sample based on adenoma or carcinoma.

$(FIGS)/FigureS2.pdf : $(TABLES)/time_datatable.csv\
code/Run_FigureS2.R
	R -e "source('code/Run_FigureS2.R')"



#exploratory/CommonFeatures_RF_model_100.RData : 
#	mkdir $(CODE)/common
#	R -e "source('code/Run_common_feature_model.R')"
#	bash $(CODE)/createDuplicates_commonVars.sh
#	bash $(CODE)/create_commonVars_pbs.sh
#	bash $(CODE)/qsubmission_commonVars.sh




################################################################################
#
# Part 5: Pull it all together
#
# Render the manuscript
#
################################################################################


write.paper : $(FINAL)/manuscript_outline_20161024.Rmd\
		$(TABLES)/mod_metadata/good_metaf_final.csv\
		$(TABLES)/alpha_table_summary.csv\
		$(TABLES)/time_pvalues.csv\ 
		$(FIGS)/Figure1.pdf $(FIGS)/Figure2.pdf\
		$(FIGS)/Figure3.pdf $(FIGS)/Figure4.pdf\
		$(FIGS)/FigureS1.pdf $(FIGS)/FigureS2.pdf\
		$(FIGS)/FigureS3.pdf $(FIGS)/FigureS4.pdf\
		$(FIGS)/FigureS5.pdf
	R -e "source('code/Run_render_paper.R')"

