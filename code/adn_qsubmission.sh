#!/bin/sh

submission="adn_run_100_RF.pbs adn_run_10_RF.pbs adn_run_11_RF.pbs adn_run_12_RF.pbs 
adn_run_13_RF.pbs adn_run_14_RF.pbs adn_run_15_RF.pbs adn_run_16_RF.pbs adn_run_17_RF.pbs 
adn_run_18_RF.pbs adn_run_19_RF.pbs adn_run_1_RF.pbs adn_run_20_RF.pbs adn_run_21_RF.pbs 
adn_run_22_RF.pbs adn_run_23_RF.pbs adn_run_24_RF.pbs adn_run_25_RF.pbs adn_run_26_RF.pbs 
adn_run_27_RF.pbs adn_run_28_RF.pbs adn_run_29_RF.pbs adn_run_2_RF.pbs adn_run_30_RF.pbs 
adn_run_31_RF.pbs adn_run_32_RF.pbs adn_run_33_RF.pbs adn_run_34_RF.pbs adn_run_35_RF.pbs 
adn_run_36_RF.pbs adn_run_37_RF.pbs adn_run_38_RF.pbs adn_run_39_RF.pbs adn_run_3_RF.pbs 
adn_run_40_RF.pbs adn_run_41_RF.pbs adn_run_42_RF.pbs adn_run_43_RF.pbs adn_run_44_RF.pbs 
adn_run_45_RF.pbs adn_run_46_RF.pbs adn_run_47_RF.pbs adn_run_48_RF.pbs adn_run_49_RF.pbs 
adn_run_4_RF.pbs adn_run_50_RF.pbs adn_run_51_RF.pbs adn_run_52_RF.pbs adn_run_53_RF.pbs 
adn_run_54_RF.pbs adn_run_55_RF.pbs adn_run_56_RF.pbs adn_run_57_RF.pbs adn_run_58_RF.pbs 
adn_run_59_RF.pbs adn_run_61_RF.pbs adn_run_62_RF.pbs adn_run_63_RF.pbs adn_run_64_RF.pbs 
adn_run_65_RF.pbs adn_run_66_RF.pbs adn_run_67_RF.pbs adn_run_68_RF.pbs adn_run_69_RF.pbs 
adn_run_6_RF.pbs adn_run_70_RF.pbs adn_run_71_RF.pbs adn_run_72_RF.pbs adn_run_73_RF.pbs 
adn_run_74_RF.pbs adn_run_75_RF.pbs adn_run_5_RF.pbs adn_run_76_RF.pbs adn_run_77_RF.pbs 
adn_run_78_RF.pbs adn_run_79_RF.pbs adn_run_7_RF.pbs adn_run_80_RF.pbs adn_run_81_RF.pbs 
adn_run_82_RF.pbs adn_run_60_RF.pbs adn_run_83_RF.pbs adn_run_84_RF.pbs adn_run_85_RF.pbs 
adn_run_86_RF.pbs adn_run_87_RF.pbs adn_run_88_RF.pbs adn_run_89_RF.pbs adn_run_8_RF.pbs 
adn_run_90_RF.pbs adn_run_91_RF.pbs adn_run_92_RF.pbs adn_run_93_RF.pbs adn_run_94_RF.pbs 
adn_run_95_RF.pbs adn_run_96_RF.pbs adn_run_97_RF.pbs adn_run_98_RF.pbs adn_run_99_RF.pbs 
adn_run_9_RF.pbs"

x=1

for s in ${submission}
do
	if [ $(( $x % 2)) -eq 0 ]; then
		
		qsub -W depend=afterok:$run_job code/adn/$s
	else
		run_job=`qsub code/adn/$s`
	fi
	
	x=`expr $x + 1`
	sleep 1

done