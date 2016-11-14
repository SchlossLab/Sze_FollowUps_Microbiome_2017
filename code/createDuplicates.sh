#!/bin/sh

for j in {1..100}
do
	cp run_1to10_RF.R run_${j}_RF.R
	sed -i -e "s/i\s=\s1/i = ${j}/g" run_${j}_RF.R
done