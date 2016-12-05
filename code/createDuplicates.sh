#!/bin/sh

for j in {1..100}
do
	cp code/reference_run_RF.R code/wfit/run_${j}_RF.R
	sed -i -e "s/i\s=\s1/i = ${j}/g" code/wfit/run_${j}_RF.R
done