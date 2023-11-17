#!/bin/bash

for ((i=0;i<112;i++))
do
echo $i
mac_name="Pixel_$i.mac"
energy=10000000 #10GeV
nohup ./CrystalEye ../macros/Pixel/${mac_name} ${energy} "Pixel_${i}" > ../../../Log/Pixel_$i.log 2>&1 &
echo "Submit: ${mac_name}"

done
