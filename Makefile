# Set local variables
REFS = data/references
FIGS = results/figures
TABLES = results/tables
PROC = data/process
FINAL = submission
CODE = code


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
	wget -N http://www.mothur.org/w/images/8/88/Trainset14_032015.pds.tgz
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

$(TABLES)/mod_metadata/good_metaf_final.csv :
	R -e "source('code/make_metadata_tables.R')"


$(TABLES)/alpha_table_summary.csv : 
	R -e "source('code/Run_Alpha_Diversity_tests.R')"

$(FIGS)/Figure1.pdf : 
	R -e "source('code/Run_change_theta_Fit.R')"
	R -e "source('code/Run_Beta_Diversity_tests.R')"
	R -e "source('code/Run_Figure1.R')"

$(FIGS)/Figure2.pdf :
	R -e "source('code/Run_potential_cancer_specific_OTUs.R')"
	R -e "source('code/Run_Figure2.R')"

$(TABLES)/time_pvalues.csv : 
	R -e "source('code/Run_Supplemental_time_table.R')"
	
exploratory/RF_model_100.RData : 
	mkdir $(CODE)/wfit
	R -e "source('code/setup_RF_test.R')"
	bash $(CODE)/createDuplicates.sh
	bash $(CODE)/create_pbs.sh
	bash $(CODE)/qsubmission.sh

exploratory/Reducedfeatures_RF_model_100.RData : 
	R -e "source('code/Run_reduce_feature_lesion_model.R')"
	bash $(CODE)/createDuplicates_reducedVars.sh
	bash $(CODE)/create_reducedVars_pbs.sh
	bash $(CODE)/qsubmission_reducedVars.sh

$(TABLES)/reduced_IF_follow_up_probability_summary.csv : 
	R -e "source('code/Run_Get_Imp_OTUs.R')"
	R -e "source('code/Run_combine_IF_aggregate_model.R')"
	R -e "source('code/Run_IF_best_model.R')"
	R -e "source('code/Run_reduce_feature_IF_model.R')"
	R -e "source('code/Run_combine_reduced_IF_aggregate_model.R')"
	R -e "source('code/Run_IF_reduced_best_model.R')"

$(FIGS)/Figure3.pdf : 
	#Collects the needed data to generate figure 3
	R -e "source('code/Run_Combine_Testing_pull_imp_OTUs.R')"
	R -e "source('code/Run_combine_aggregate_reduced_model.R')"
	#Generates complete model built on all data and updates tables
	R -e "source('code/Run_Create_Use_Best_Model.R')"
	R -e "source('code/Run_reduced_best_model.R')"
	#Creates the actual Figure 3
	R -e "source('code/Run_Figure3.R')"
	tiff2pdf -z -o results/figures/Figure3.pdf results/figures/Figure3.tiff
	rm results/figures/Figure3.tiff

$(FIGS)/Figure4.pdf : 
	R -e "source('code/Run_Figure4.R')"
	rm Rplots.pdf

$(FIGS)/FigureS1.pdf :
	R -e "source('code/Run_wilcoxson_all.R')"
	R -e "source('code/Run_FigureS1.R')"

$(FIGS)/FigureS2.pdf : 
	R -e "source('code/Run_FigureS2.R')"

$(TABLES)/pvalue_IF_lesion_common_imp_vars.csv : 
	R -e "source('code/Run_ID_imp_OTUs.R')"
	R -e "source('code/Run_Compare_models.R')"

$(FIGS)/FigureS3.pdf : 
	R -e "source('code/Run_FigureS3.R')"

$(FIGS)/FigureS4.pdf : 
	R -e "source('code/Run_FigureS4.R')"

$(FIGS)/FigureS5.pdf : 
	R -e "source('code/Run_FigureS5.R')"


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
		$(FIGS)/Figure5.pdf $(FIGS)/Figure6.pdf\
		$(FIGS)/FigureS1.pdf $(FIGS)/FigureS2.pdf\
		$(FIGS)/FigureS3.pdf $(FIGS)/FigureS4.pdf\
		$(FIGS)/FigureS5.pdf
	R -e "source('code/Run_render_paper.R')"

