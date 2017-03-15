#!/bin/sh

for j in {1..100}
do
	cp code/adn_reference_run_RF.R code/adn/adn_run_${j}_RF.R
	sed -i -e "s/i\s=\s1/i = ${j}/g" code/adn/adn_run_${j}_RF.R
done