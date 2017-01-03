#!/bin/sh

for j in {1..100}
do
	cp code/reference_run_reduced_feature_RF.R code/wfit/reduced_run_${j}_RF.R
	sed -i -e "s/i\s=\s1/i = ${j}/g" code/wfit/reduced_run_${j}_RF.R
done