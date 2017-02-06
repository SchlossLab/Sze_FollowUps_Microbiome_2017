#!/bin/sh

submission="reduced_run_100_RF.pbs reduced_run_10_RF.pbs reduced_run_11_RF.pbs reduced_run_12_RF.pbs 
reduced_run_13_RF.pbs reduced_run_14_RF.pbs reduced_run_15_RF.pbs reduced_run_16_RF.pbs reduced_run_17_RF.pbs 
reduced_run_18_RF.pbs reduced_run_19_RF.pbs reduced_run_1_RF.pbs reduced_run_20_RF.pbs reduced_run_21_RF.pbs 
reduced_run_22_RF.pbs reduced_run_23_RF.pbs reduced_run_24_RF.pbs reduced_run_25_RF.pbs reduced_run_26_RF.pbs 
reduced_run_27_RF.pbs reduced_run_28_RF.pbs reduced_run_29_RF.pbs reduced_run_2_RF.pbs reduced_run_30_RF.pbs 
reduced_run_31_RF.pbs reduced_run_32_RF.pbs reduced_run_33_RF.pbs reduced_run_34_RF.pbs reduced_run_35_RF.pbs 
reduced_run_36_RF.pbs reduced_run_37_RF.pbs reduced_run_38_RF.pbs reduced_run_39_RF.pbs reduced_run_3_RF.pbs 
reduced_run_40_RF.pbs reduced_run_41_RF.pbs reduced_run_42_RF.pbs reduced_run_43_RF.pbs reduced_run_44_RF.pbs 
reduced_run_45_RF.pbs reduced_run_46_RF.pbs reduced_run_47_RF.pbs reduced_run_48_RF.pbs reduced_run_49_RF.pbs 
reduced_run_4_RF.pbs reduced_run_50_RF.pbs reduced_run_51_RF.pbs reduced_run_52_RF.pbs reduced_run_53_RF.pbs 
reduced_run_54_RF.pbs reduced_run_55_RF.pbs reduced_run_56_RF.pbs reduced_run_57_RF.pbs reduced_run_58_RF.pbs 
reduced_run_59_RF.pbs reduced_run_61_RF.pbs reduced_run_62_RF.pbs reduced_run_63_RF.pbs reduced_run_64_RF.pbs 
reduced_run_65_RF.pbs reduced_run_66_RF.pbs reduced_run_67_RF.pbs reduced_run_68_RF.pbs reduced_run_69_RF.pbs 
reduced_run_6_RF.pbs reduced_run_70_RF.pbs reduced_run_71_RF.pbs reduced_run_72_RF.pbs reduced_run_73_RF.pbs 
reduced_run_74_RF.pbs reduced_run_75_RF.pbs reduced_run_5_RF.pbs reduced_run_76_RF.pbs reduced_run_77_RF.pbs 
reduced_run_78_RF.pbs reduced_run_79_RF.pbs reduced_run_7_RF.pbs reduced_run_80_RF.pbs reduced_run_81_RF.pbs 
reduced_run_82_RF.pbs reduced_run_60_RF.pbs reduced_run_83_RF.pbs reduced_run_84_RF.pbs reduced_run_85_RF.pbs 
reduced_run_86_RF.pbs reduced_run_87_RF.pbs reduced_run_88_RF.pbs reduced_run_89_RF.pbs reduced_run_8_RF.pbs 
reduced_run_90_RF.pbs reduced_run_91_RF.pbs reduced_run_92_RF.pbs reduced_run_93_RF.pbs reduced_run_94_RF.pbs 
reduced_run_95_RF.pbs reduced_run_96_RF.pbs reduced_run_97_RF.pbs reduced_run_98_RF.pbs reduced_run_99_RF.pbs 
reduced_run_9_RF.pbs"

x=1

for s in ${submission}
do
	if [ $(( $x % 2)) -eq 0 ]; then
		
		qsub -W depend=afterok:$reduced_run_job code/reduced/$s
	else
		reduced_run_job=`qsub code/reduced/$s`
	fi
	
	x=`expr $x + 1`
	sleep 1

done