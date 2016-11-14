#!/bin/sh

for j in {1..100}
do
	cp RF_test_1to10.pbs run_${j}_RF.pbs
	sed -i -e "s/run_1_RF\.R/run_${j}_RF\.R/g" run_${j}_RF.pbs
done