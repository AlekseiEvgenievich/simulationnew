#!/bin/bash

couts=0
for (( t=0;t<=0;t+=5))
do

fmin=10
fmax=100

macname="circle_theta_Pow/gamma_circle_Theta${t}_E${fmin}_${fmax}.mac"
for (( i=0;i<1;i++))
do
		(( couts++ ))
		argv="Theta${t}_E${fmin}_${fmax}_${i}"
		nohup ./CrystalEye ../macros/${macname} ${argv} ${argv} > ../../../Log/Eng_${argv}.log 2>&1 &
		echo "Counts:$couts Submit Energy: ${argv} keV!  Theta: $t  ......."

done

done
