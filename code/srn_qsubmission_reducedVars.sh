#!/bin/sh

submission="srn_reduced_run_100_RF.pbs srn_reduced_run_10_RF.pbs srn_reduced_run_11_RF.pbs srn_reduced_run_12_RF.pbs 
srn_reduced_run_13_RF.pbs srn_reduced_run_14_RF.pbs srn_reduced_run_15_RF.pbs srn_reduced_run_16_RF.pbs srn_reduced_run_17_RF.pbs 
srn_reduced_run_18_RF.pbs srn_reduced_run_19_RF.pbs srn_reduced_run_1_RF.pbs srn_reduced_run_20_RF.pbs srn_reduced_run_21_RF.pbs 
srn_reduced_run_22_RF.pbs srn_reduced_run_23_RF.pbs srn_reduced_run_24_RF.pbs srn_reduced_run_25_RF.pbs srn_reduced_run_26_RF.pbs 
srn_reduced_run_27_RF.pbs srn_reduced_run_28_RF.pbs srn_reduced_run_29_RF.pbs srn_reduced_run_2_RF.pbs srn_reduced_run_30_RF.pbs 
srn_reduced_run_31_RF.pbs srn_reduced_run_32_RF.pbs srn_reduced_run_33_RF.pbs srn_reduced_run_34_RF.pbs srn_reduced_run_35_RF.pbs 
srn_reduced_run_36_RF.pbs srn_reduced_run_37_RF.pbs srn_reduced_run_38_RF.pbs srn_reduced_run_39_RF.pbs srn_reduced_run_3_RF.pbs 
srn_reduced_run_40_RF.pbs srn_reduced_run_41_RF.pbs srn_reduced_run_42_RF.pbs srn_reduced_run_43_RF.pbs srn_reduced_run_44_RF.pbs 
srn_reduced_run_45_RF.pbs srn_reduced_run_46_RF.pbs srn_reduced_run_47_RF.pbs srn_reduced_run_48_RF.pbs srn_reduced_run_49_RF.pbs 
srn_reduced_run_4_RF.pbs srn_reduced_run_50_RF.pbs srn_reduced_run_51_RF.pbs srn_reduced_run_52_RF.pbs srn_reduced_run_53_RF.pbs 
srn_reduced_run_54_RF.pbs srn_reduced_run_55_RF.pbs srn_reduced_run_56_RF.pbs srn_reduced_run_57_RF.pbs srn_reduced_run_58_RF.pbs 
srn_reduced_run_59_RF.pbs srn_reduced_run_61_RF.pbs srn_reduced_run_62_RF.pbs srn_reduced_run_63_RF.pbs srn_reduced_run_64_RF.pbs 
srn_reduced_run_65_RF.pbs srn_reduced_run_66_RF.pbs srn_reduced_run_67_RF.pbs srn_reduced_run_68_RF.pbs srn_reduced_run_69_RF.pbs 
srn_reduced_run_6_RF.pbs srn_reduced_run_70_RF.pbs srn_reduced_run_71_RF.pbs srn_reduced_run_72_RF.pbs srn_reduced_run_73_RF.pbs 
srn_reduced_run_74_RF.pbs srn_reduced_run_75_RF.pbs srn_reduced_run_5_RF.pbs srn_reduced_run_76_RF.pbs srn_reduced_run_77_RF.pbs 
srn_reduced_run_78_RF.pbs srn_reduced_run_79_RF.pbs srn_reduced_run_7_RF.pbs srn_reduced_run_80_RF.pbs srn_reduced_run_81_RF.pbs 
srn_reduced_run_82_RF.pbs srn_reduced_run_60_RF.pbs srn_reduced_run_83_RF.pbs srn_reduced_run_84_RF.pbs srn_reduced_run_85_RF.pbs 
srn_reduced_run_86_RF.pbs srn_reduced_run_87_RF.pbs srn_reduced_run_88_RF.pbs srn_reduced_run_89_RF.pbs srn_reduced_run_8_RF.pbs 
srn_reduced_run_90_RF.pbs srn_reduced_run_91_RF.pbs srn_reduced_run_92_RF.pbs srn_reduced_run_93_RF.pbs srn_reduced_run_94_RF.pbs 
srn_reduced_run_95_RF.pbs srn_reduced_run_96_RF.pbs srn_reduced_run_97_RF.pbs srn_reduced_run_98_RF.pbs srn_reduced_run_99_RF.pbs 
srn_reduced_run_9_RF.pbs"

x=1

for s in ${submission}
do
	if [ $(( $x % 2)) -eq 0 ]; then
		
		qsub -W depend=afterok:$reduced_run_job code/reduced_srn/$s
	else
		reduced_run_job=`qsub code/reduced_srn/$s`
	fi
	
	x=`expr $x + 1`
	sleep 1

done