#!/bin/bash
phi=0
counts=0
#for (( phi=0;phi<=360;phi+=45 ))
#for (( phi=0;phi<=15;phi+=5 ))
#do

for (( theta=0;theta<=90;theta+=5 ))
do
angle="${theta}_${phi}"
#angle="${theta}"
#macname="circle_theta_phi_uniform/gamma_circle_Angle${angle}_E10_50000.mac"
#macname="circle_theta_phi_pow/gamma_circle_Angle${angle}_E10_50000.mac"
#macname="alexei/gamma_pow_Angle${angle}_Energy10_50000.mac"
#macname="aleksei/gamma_circle_${angle}.mac"
macname="aleksei/60kev/gamma_circle_${angle}.mac"

for (( i=0;i<1;i++))
do
		(( counts++ ))
    nohup ./CrystalEye ../macros/${macname}  "${angle}_${i}" > ../Log/Angle_${angle}_${i}.log 2>&1 &
#echo "Counts:$counts Submit Theta: ${theta}  Phi: ${phi}  ......."

done

done

#done
