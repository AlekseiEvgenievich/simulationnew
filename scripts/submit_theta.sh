#!/bin/bash

couts=0
for (( t=0;t<=90;t+=5))
do
macname="circle_theta/gamma_circle_$t.mac"

#feng=10
#fgap=10
#for (( i=0;i<9;i++))
#do
#		(( couts++ ))
#		argv=$(( feng + i*fgap ))
#		nohup ./CrystalEye ../macros/${macname} ${argv} "${argv}_$t" > ../../../Log/Eng_${argv}_$t.log 2>&1 &
#		echo "Counts:$couts Submit Energy: ${argv} keV!  Theta: $t  ......."
#
#done

feng=100
fgap=100
for (( i=0;i<9;i++))
do
		(( couts++ ))
		argv=$(( feng + i*fgap ))
		nohup ./CrystalEye ../macros/${macname} ${argv} "${argv}_$t" > ../../../Log/Eng_${argv}_$t.log 2>&1 &
		echo "Counts:$couts Submit Energy: ${argv} keV!  Theta: $t  ......."

done

#feng=1000
#fgap=1000
#for (( i=0;i<9;i++))
#do
#		argv=$(( feng + i*fgap ))
#		nohup ./CrystalEye ../macros/${macname} ${argv} "${argv}_$t" > ../../../Log/Eng_${argv}_$t.log 2>&1 &
#		echo "Counts:$couts Submit Energy: ${argv} keV!  Theta: $t  ......."
#
#done

#feng=10000
#fgap=10000
#for (( i=0;i<=9;i++))
#do
#		argv=$(( feng + i*fgap ))
#		nohup ./CrystalEye ../macros/${macname} ${argv} "${argv}_$t" > ../../../Log/Eng_${argv}_$t.log 2>&1 &
#		echo "Counts:$couts Submit Energy: ${argv} keV!  Theta: $t  ......."
#
#done

done
