# Set local variables
REFS = data/references
FIGS = results/figures
TABLES = results/tables
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

$(PROC)/unmatched.% : 
	bash $(CODE)/mothur.batch

$(PROC)/final.% : 
	bash $(CODE)/mothurCluster.batch

################################################################################
#
# Part 3: Figure and table generation
#
#	Run scripts to generate figures and tables
#
################################################################################

$(TABLES)/mod_metadata/good_metaf_final.csv : $(METADATA)/followUps_metadata.txt\
$(METADATA)/initials_metadata.tsv $(METADATA)/followUp_outcome_data.csv\
code/make_metadata_tables.R 
	R -e "source('code/make_metadata_tables.R')"

$(TABLES)/alpha_table_summary.csv : $(TABLES)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.groups.ave-std.summary code/Run_Alpha_Diversity_tests.R
	R -e "source('code/Run_Alpha_Diversity_tests.R')"

$(TABLES)/difference_table.csv : $(PROC)/final.thetayc.0.03.lt.ave.dist\
$(TABLES)/mod_metadata/metaF_final.csv code/Run_change_theta_Fit.R
	R -e "source('code/Run_change_theta_Fit.R')"

$(TABLES)/thetayc_% : $(PROC)/final.thetayc.0.03.lt.ave.dist\
$(TABLES)/mod_metadata/metaF_final.csv code/Run_Beta_Diversity_tests.R
	R -e "source('code/Run_Beta_Diversity_tests.R')"

$(TABLES)/adn_crc_maybe_pvalue_summary.csv : $(PROC)/final.taxonomy\
$(TABLES)/mod_metadata/good_metaf_final.csv $(PROC)/final.shared\
code/Run_potential_cancer_specific_OTUs.R
	R -e "source('code/Run_potential_cancer_specific_OTUs.R')"

$(TABLES)/time_pvalues.csv : $(PROC)/final.thetayc.0.03.lt.ave.dist\
$(TABLES)/mod_metadata/metaI_final.csv $(TABLES)/mod_metadata/metaF_final.csv\
$(TABLES)/mod_metadata/good_metaf_final.csv code/Run_Supplemental_time_table.R
	R -e "source('code/Run_Supplemental_time_table.R')"

exploratory/RF_model_100.RData : $(PROC)/final.0.03.subsample.shared\
$(TABLES)/mod_metadata/metaI_final.csv $(TABLES)/mod_metadata/good_metaf_final.csv\
code/reference_run_RF.R code/RF_reference.pbs code/setup_RF_test.R\
$(CODE)/createDuplicates.sh $(CODE)/create_pbs.sh $(CODE)/qsubmission.sh
	mkdir $(CODE)/full
	R -e "source('code/setup_RF_test.R')"
	bash $(CODE)/createDuplicates.sh
	bash $(CODE)/create_pbs.sh
	bash $(CODE)/qsubmission.sh

exploratory/rocs.RData : code/Run_Combine_Testing_pull_imp_OTUs.R
	R -e "source('code/Run_Combine_Testing_pull_imp_OTUs.R')"

exploratory/Reducedfeatures_RF_model_100.RData : $(TABLES)/full_test_data.csv\
$(TABLES)/rf_wCV_imp_vars_summary.csv code/RF_reduced_vars_reference.pbs\
code/reference_run_reduced_feature_RF.R code/Run_reduce_feature_lesion_model.R\
$(CODE)/createDuplicates_reducedVars.sh $(CODE)/create_reducedVars_pbs.sh\
$(CODE)/qsubmission_reducedVars.sh
	mkdir $(CODE)/reduced
	R -e "source('code/Run_reduce_feature_lesion_model.R')"
	bash $(CODE)/createDuplicates_reducedVars.sh
	bash $(CODE)/create_reducedVars_pbs.sh
	bash $(CODE)/qsubmission_reducedVars.sh

exploratory/RF_model_Imp_OTU.RData : $(TABLES)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/Run_Get_Imp_OTUs.R
	R -e "source('code/Run_Get_Imp_OTUs.R')"

$(TABLES)/IF_model_top_vars_MDA_Summary.csv : exploratory/RF_model_Imp_OTU.RData\
code/Run_combine_IF_aggregate_model.R
	R -e "source('code/Run_combine_IF_aggregate_model.R')"

$(TABLES)/IF_follow_up_probability_summary.csv : $(TABLES)/IF_test_tune_data.csv\
$(TABLES)/IF_ROC_model_summary.csv $(TABLES)/IF_test_data_roc.csv\
$(TABLES)/mod_metadata/good_metaf_final.csv $(PROC)/final.0.03.subsample.shared\
code/Run_IF_best_model.R
	R -e "source('code/Run_IF_best_model.R')"

exploratory/IF_reduced_RF_model_Imp_OTU.RData : $(TABLES)/IF_test_tune_data.csv\
$(TABLES)/IF_rf_wCV_imp_vars_summary.csv code/Run_reduce_feature_IF_model.R
	R -e "source('code/Run_reduce_feature_IF_model.R')"

$(TABLES)/reduced_IF_model_top_vars_MDA_Summary.csv : exploratory/IF_reduced_RF_model_Imp_OTU.RData\
Run_combine_reduced_IF_aggregate_model.R
	R -e "source('code/Run_combine_reduced_IF_aggregate_model.R')"

$(TABLES)/reduced_IF_follow_up_probability_summary.csv : $(TABLES)/reduced_IF_test_tune_data.csv\
$(TABLES)/reduced_IF_ROC_model_summary.csv $(TABLES)/reduced_IF_test_data_roc.csv\
$(TABLES)/mod_metadata/good_metaf_final.csv $(PROC)/final.0.03.subsample.shared\
code/Run_IF_reduced_best_model.R
	R -e "source('code/Run_IF_reduced_best_model.R')"

$(TABLES)/reduced_lesion_model_top_vars_MDA_Summary : $(TABLES)/reduced_test_tune_data.csv\
$(TABLES)/reduced_test_data_splits.csv code/Run_combine_aggregate_reduced_model.R
	#Collects the needed data to generate figure 3
	R -e "source('code/Run_combine_aggregate_reduced_model.R')"

$(TABLES)/roc_pvalue_summary.csv : exploratory/rocs.RData\
$(TABLES)/full_test_data.csv $(TABLES)/ROC_model_summary.csv $(TABLES)/test_data_roc.csv\
$(TABLES)/auc_summary.csv $(TABLES)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.0.03.subsample.shared $(TABLES)/follow_up_prediction_table.csv\
code/Run_Create_Use_Best_Model.R
	#Generates complete model built on all data and updates tables
	R -e "source('code/Run_Create_Use_Best_Model.R')"

$(TABLES)/reduced_follow_up_probability_summary.csv : $(TABLES)/reduced_test_tune_data.csv\
$(TABLES)/Reduced_ROC_model_summary.csv $(TABLES)/reduced_test_data_roc.csv\
$(TABLES)/reduced_auc_summary.csv $(TABLES)/mod_metadata/good_metaf_final.csv\
$(PROC)/final.0.03.subsample.shared code/Run_reduced_best_model.R
	R -e "source('code/Run_reduced_best_model.R')"

$(TABLES)/all_models_wilcox_paired_pvalue_summary.csv : $(TABLES)/follow_up_probability_summary.csv\
$(TABLES)/reduced_follow_up_probability_summary.csv $(TABLES)/IF_follow_up_probability_summary.csv\
$(TABLES)/reduced_IF_follow_up_probability_summary.csv $(TABLES)/mod_metadata/good_metaf_final.csv\
code/Run_probs_comparison.R
	R -e "source('code/Run_probs_comparison.R')"

$(TABLES)/OTU_paired_wilcoxson_test.csv : $(TABLES)/full_test_data.csv\
$(PROC)/final.0.03.subsample.shared $(TABLES)/mod_metadata/good_metaf_final.csv\
code/Run_wilcoxson_all.R
	R -e "source('code/Run_wilcoxson_all.R')"

$(TABLES)/probs_chemo_rad_pvalue_summary.csv : $(TABLES)/follow_up_probability_summary.csv\
$(TABLES)/reduced_follow_up_probability_summary.csv $(TABLES)/IF_follow_up_probability_summary.csv\
$(TABLES)/reduced_IF_follow_up_probability_summary.csv $(TABLES)/difference_table.csv\
$(TABLES)/mod_metadata/good_metaf_final.csv $(PROC)/final.groups.ave-std.summary\
code/Run_Test_Chemo_Rad.R
	R -e "source('code/Run_Test_Chemo_Rad.R')"

$(TABLES)/%_otu_tax.csv : $(PROC)/final.taxonomy $(TABLES)/IF_rf_wCV_imp_vars_summary.csv\
$(TABLES)/results/tables/rf_wCV_imp_vars_summary.csv code/Run_ID_imp_OTUs.R
	R -e "source('code/Run_ID_imp_OTUs.R')"

$(TABLES)/pvalue_IF_lesion_common_imp_vars.csv : $(TABLES)/rf_otu_tax.csv\
$(TABLES)/if_rf_otu_tax.csv $(PROC)/final.0.03.subsample.shared\
$(TABLES)/mod_metadata/good_metaf_final.csv code/Run_Compare_models.R
	R -e "source('code/Run_Compare_models.R')"

$(FIGS)/Figure1.pdf : $(TABLES)/difference_table.csv\
$(TABLES)/change_theta_fit_summary.csv $(TABLES)/thetayc_adn_IF.csv\
$(TABLES)/thetayc_crc_IF.csv $(TABLES)/beta_diver_summary.csv\
code/Run_Figure1.R
	R -e "source('code/Run_Figure1.R')"

$(FIGS)/Figure2.pdf : $(TABLES)/adn_crc_maybe_diff.csv code/Run_Figure2.R
	R -e "source('code/Run_Figure2.R')"

$(FIGS)/Figure3.pdf : $(TABLES)/reduced_test_data_roc.csv\
$(TABLES)/reduced_lesion_model_top_vars_MDA_Summary.csv\
$(TABLES)/reduced_lesion_model_top_vars_MDA.csv\
$(TABLES)/reduced_follow_up_probability_summary.csv $(PROC)/final.taxonomy code/Run_Figure3.R
	#Creates the actual Figure 3
	R -e "source('code/Run_Figure3.R')"
	tiff2pdf -z -o results/figures/Figure3.pdf results/figures/Figure3.tiff
	rm results/figures/Figure3.tiff

$(FIGS)/Figure4.pdf : $(TABLES)/reduced_IF_test_data_roc.csv\
$(TABLES)/reduced_IF_model_top_vars_MDA_Summary.csv\
$(TABLES)/reduced_IF_model_top_vars_MDA.csv\
$(TABLES)/reduced_IF_follow_up_probability_summary.csv $(PROC)/final.taxonomy code/Run_Figure4.R
	R -e "source('code/Run_Figure4.R')"
	tiff2pdf -z -o results/figures/Figure4.pdf results/figures/Figure4.tiff
	rm results/figures/Figure4.tiff

$(FIGS)/FigureS1.pdf : $(TABLES)/OTU_paired_wilcoxson_test.csv\
code/Run_FigureS1.R
	R -e "source('code/Run_FigureS1.R')"

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
# Part 4: Pull it all together
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

