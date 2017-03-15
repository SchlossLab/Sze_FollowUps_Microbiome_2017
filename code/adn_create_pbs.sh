#!/bin/sh

for j in {1..100}
do
	cp code/adn_RF_reference.pbs code/adn/adn_run_${j}_RF.pbs
	sed -i -e "s/adn_run_1_RF\.R/adn_run_${j}_RF\.R/g" code/adn/adn_run_${j}_RF.pbs
done