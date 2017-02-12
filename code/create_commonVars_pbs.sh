#!/bin/sh

for j in {1..100}
do
	cp code/RF_common_vars_reference.pbs code/common/common_run_${j}_RF.pbs
	sed -i -e "s/run_1_RF\.R/run_${j}_RF\.R/g" code/common/common_run_${j}_RF.pbs
done