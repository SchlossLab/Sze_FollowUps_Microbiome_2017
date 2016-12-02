#!/bin/sh

for j in {1..100}
do
	cp exploratory/run_1to10_RF.R exploratory/run_${j}_RF.R
	sed -i -e "s/i\s=\s1/i = ${j}/g" exploratory/run_${j}_RF.R
done