#!/bin/sh

for j in {1..100}
do
	cp code/srn_RF_reference.pbs code/srn/srn_run_${j}_RF.pbs
	sed -i -e "s/srn_run_1_RF\.R/srn_run_${j}_RF\.R/g" code/srn/srn_run_${j}_RF.pbs
done