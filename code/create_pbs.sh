#!/bin/sh

for j in {1..100}
do
	cp code/RF_reference.pbs code/full/run_${j}_RF.pbs
	sed -i -e "s/run_1_RF\.R/run_${j}_RF\.R/g" code/full/run_${j}_RF.pbs
done