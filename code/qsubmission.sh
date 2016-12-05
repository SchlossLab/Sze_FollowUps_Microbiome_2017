#!/bin/sh

submission="run_100_RF.pbs run_10_RF.pbs run_11_RF.pbs run_12_RF.pbs run_13_RF.pbs 
run_14_RF.pbs run_15_RF.pbs run_16_RF.pbs" 

submission_two="run_17_RF.pbs run_18_RF.pbs run_19_RF.pbs 
run_1_RF.pbs run_20_RF.pbs run_21_RF.pbs run_22_RF.pbs run_23_RF.pbs" 

submission_three="run_24_RF.pbs run_25_RF.pbs run_26_RF.pbs run_27_RF.pbs run_28_RF.pbs 
run_29_RF.pbs run_2_RF.pbs run_30_RF.pbs" 

submission_four="run_31_RF.pbs run_32_RF.pbs run_33_RF.pbs run_34_RF.pbs run_35_RF.pbs 
run_36_RF.pbs run_37_RF.pbs run_38_RF.pbs" 

submission_five="run_39_RF.pbs run_3_RF.pbs run_40_RF.pbs run_41_RF.pbs 
run_42_RF.pbs run_43_RF.pbs run_44_RF.pbs run_45_RF.pbs" 

submission_six="run_46_RF.pbs run_47_RF.pbs run_48_RF.pbs run_49_RF.pbs 
run_4_RF.pbs run_50_RF.pbs run_51_RF.pbs run_52_RF.pbs run_53_RF.pbs" 

submission_seven="run_54_RF.pbs run_55_RF.pbs run_56_RF.pbs run_57_RF.pbs 
run_58_RF.pbs run_59_RF.pbs" 

submission_eight="run_61_RF.pbs run_62_RF.pbs run_63_RF.pbs 
run_64_RF.pbs run_65_RF.pbs run_66_RF.pbs run_67_RF.pbs run_68_RF.pbs"

submission_nine="run_69_RF.pbs run_6_RF.pbs run_70_RF.pbs run_71_RF.pbs 
run_72_RF.pbs run_73_RF.pbs run_74_RF.pbs run_75_RF.pbs run_5_RF.pbs" 

submission_ten="run_76_RF.pbs run_77_RF.pbs run_78_RF.pbs run_79_RF.pbs 
run_7_RF.pbs run_80_RF.pbs run_81_RF.pbs run_82_RF.pbs run_60_RF.pbs" 

submission_eleven="run_83_RF.pbs run_84_RF.pbs run_85_RF.pbs run_86_RF.pbs 
run_87_RF.pbs run_88_RF.pbs run_89_RF.pbs run_8_RF.pbs"

submission_twelve="run_90_RF.pbs run_91_RF.pbs run_92_RF.pbs run_93_RF.pbs run_94_RF.pbs run_95_RF.pbs 
run_96_RF.pbs run_97_RF.pbs run_98_RF.pbs run_99_RF.pbs run_9_RF.pbs"

x=1

for s in ${submission}
do
	if [ $(( $x % 2)) -eq 0 ]; then
		echo $x
		qsub -W depend=afterok:$run_job wfit/$s
	else
		run_job=`qsub wfit/$s`
	fi
	
	x=`expr $x + 1`
	sleep 1
done