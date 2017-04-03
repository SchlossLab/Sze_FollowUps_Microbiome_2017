#!/bin/sh

submission="common_run_100_RF.pbs common_run_10_RF.pbs common_run_11_RF.pbs common_run_12_RF.pbs 
common_run_13_RF.pbs common_run_14_RF.pbs common_run_15_RF.pbs common_run_16_RF.pbs common_run_17_RF.pbs 
common_run_18_RF.pbs common_run_19_RF.pbs common_run_1_RF.pbs common_run_20_RF.pbs common_run_21_RF.pbs 
common_run_22_RF.pbs common_run_23_RF.pbs common_run_24_RF.pbs common_run_25_RF.pbs common_run_26_RF.pbs 
common_run_27_RF.pbs common_run_28_RF.pbs common_run_29_RF.pbs common_run_2_RF.pbs common_run_30_RF.pbs 
common_run_31_RF.pbs common_run_32_RF.pbs common_run_33_RF.pbs common_run_34_RF.pbs common_run_35_RF.pbs 
common_run_36_RF.pbs common_run_37_RF.pbs common_run_38_RF.pbs common_run_39_RF.pbs common_run_3_RF.pbs 
common_run_40_RF.pbs common_run_41_RF.pbs common_run_42_RF.pbs common_run_43_RF.pbs common_run_44_RF.pbs 
common_run_45_RF.pbs common_run_46_RF.pbs common_run_47_RF.pbs common_run_48_RF.pbs common_run_49_RF.pbs 
common_run_4_RF.pbs common_run_50_RF.pbs common_run_51_RF.pbs common_run_52_RF.pbs common_run_53_RF.pbs 
common_run_54_RF.pbs common_run_55_RF.pbs common_run_56_RF.pbs common_run_57_RF.pbs common_run_58_RF.pbs 
common_run_59_RF.pbs common_run_61_RF.pbs common_run_62_RF.pbs common_run_63_RF.pbs common_run_64_RF.pbs 
common_run_65_RF.pbs common_run_66_RF.pbs common_run_67_RF.pbs common_run_68_RF.pbs common_run_69_RF.pbs 
common_run_6_RF.pbs common_run_70_RF.pbs common_run_71_RF.pbs common_run_72_RF.pbs common_run_73_RF.pbs 
common_run_74_RF.pbs common_run_75_RF.pbs common_run_5_RF.pbs common_run_76_RF.pbs common_run_77_RF.pbs 
common_run_78_RF.pbs common_run_79_RF.pbs common_run_7_RF.pbs common_run_80_RF.pbs common_run_81_RF.pbs 
common_run_82_RF.pbs common_run_60_RF.pbs common_run_83_RF.pbs common_run_84_RF.pbs common_run_85_RF.pbs 
common_run_86_RF.pbs common_run_87_RF.pbs common_run_88_RF.pbs common_run_89_RF.pbs common_run_8_RF.pbs 
common_run_90_RF.pbs common_run_91_RF.pbs common_run_92_RF.pbs common_run_93_RF.pbs common_run_94_RF.pbs 
common_run_95_RF.pbs common_run_96_RF.pbs common_run_97_RF.pbs common_run_98_RF.pbs common_run_99_RF.pbs 
common_run_9_RF.pbs"

x=1

for s in ${submission}
do
	if [ $(( $x % 2)) -eq 0 ]; then
		
		qsub -W depend=afterok:$common_run_job code/common/$s
	else
		common_run_job=`qsub code/common/$s`
	fi
	
	x=`expr $x + 1`
	sleep 1

done