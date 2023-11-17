#!/bin/bash


#eng=10 # 10keV
#macname="Pixel_gamma/Pixel_109.mac"
#end_name="Pixel_109_${eng}"

#eng=10
#macname="albedo_10s_10keV.mac"
#end_name="albedo_10s_10keV"
#nohup ./CrystalEye ../macros/${macname} ${eng} "${end_name}" > ../../../Log/${end_name}.log 2>&1 &
#echo " Submit Energy: ${eng} keV!   ......."

#eng=10
#macname="cosmic_10s_10keV_2.mac"
#end_name="cosmic_10s_10keV_2"
#nohup ./CrystalEye ../macros/background/${macname} ${eng} "${end_name}" > ../../../Log/${end_name}.log 2>&1 &
#echo " Submit Energy: ${eng} keV!   ......."

#eng=10
#macname="albedo_10s_10keV_2.mac"
#end_name="albedo_10s_10keV_2"
#nohup ./CrystalEye ../macros/background/${macname} ${eng} "${end_name}" > ../../../Log/${end_name}.log 2>&1 &
#echo " Submit Energy: ${eng} keV!   ......."

eng=10
macname="gamma_uniformity.mac"
end_name="gamma_uniformity"
nohup ./CrystalEye ../macros/background/${macname} ${eng} "${end_name}" > ../../../Log/${end_name}.log 2>&1 &
echo " Submit Energy: ${eng} keV!   ......."
