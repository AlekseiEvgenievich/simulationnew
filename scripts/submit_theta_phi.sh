#!/bin/bash

counts=0
feng=10000
fgap=10000

#for (( phi=0;phi<=360;phi+=45 ))
for (( phi=0;phi<=0;phi+=45 ))
do

for (( theta=0;theta<=90;theta+=5 ))
do
angle="${theta}_${phi}"
macname="circle_theta_phi/gamma_circle_$angle.mac"

for (( i=0;i<9;i++))
do
		(( counts++ ))
		eng=$(( feng + i*fgap ))
#nohup ./CrystalEye ../macros/${macname} ${eng} "${eng}_${angle}" > ../../../Log/Eng_${eng}_${angle}.log 2>&1 &
		echo "Counts:$counts Submit Energy: ${eng} keV!  Theta: ${theta}  Phi: ${phi}  ......."

done

done

done
