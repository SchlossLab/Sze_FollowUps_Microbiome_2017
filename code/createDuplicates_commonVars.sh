#!/bin/sh

for j in {1..100}
do
	cp code/reference_run_common_feature_RF.R code/common/common_run_${j}_RF.R
	sed -i -e "s/i\s=\s1/i = ${j}/g" code/common/common_run_${j}_RF.R
done