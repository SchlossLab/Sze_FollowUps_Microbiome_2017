#!/bin/sh

for j in {1..100}
do
	cp code/crc_RF_reference.pbs code/crc/crc_run_${j}_RF.pbs
	sed -i -e "s/crc_run_1_RF\.R/crc_run_${j}_RF\.R/g" code/crc/crc_run_${j}_RF.pbs
done