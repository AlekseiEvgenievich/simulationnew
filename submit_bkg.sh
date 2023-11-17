#!/bin/bash

counts=0
macname="gamma_pow"
#macname="neutron_pow"
#macname="protonSAA_pow"
#macname="protonPrimary_pow"

for (( i=0;i<10;i++ ))
do
 (( counts ))
 nohup ./CrystalEye ../macros/background/${macname}.mac  "${macname}_${i}" > ../../../Log/Bkg_${macname}_${i}.log 2>&1 &
 echo "Submit Job: ${macname}_$i ......"

done

