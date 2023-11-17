#!/bin/bash

MacName="../macros/CNAF/Gamma_Map_Pow/Gamma_Pow_"
JobN=1

if [ ! -d "./run" ]; then
  mkdir run
fi
if [ ! -d "./out" ]; then
  mkdir out
fi
if [ ! -d "./err" ]; then
  mkdir err
fi
if [ ! -d "./cmd" ]; then
  mkdir cmd
fi
if [ ! -d "./log" ]; then
  mkdir log
fi

condor_out="./out"
condor_err="./err"
condor_log="./log"
condor_cmd="./cmd"
condor_sh="./run"

n=0
fTime=`date +%s`
echo "===>  Time Code is ${fTime}  <==="

for (( phi=0;phi<=360;phi+=5 ))
do
for (( theta=0;theta<=90;theta+=5 ))
do
angle="T${theta}_P${phi}"
Macros="${MacName}${angle}_Energy10_50000.mac"

for (( jobi=0;jobi<${JobN};jobi++))
do
fSeed=`expr $fTime + $n`
outName=${angle}_${jobi}

#--------Creat cmd file------------------------------
if test -s $condor_cmd/${outName}.cmd; then
#echo "CMD File is not Empty"
rm  $condor_cmd/${outName}.cmd
fi
cat >> $condor_cmd/${outName}.cmd <<EOF
Universe             = vanilla
Notification         = Never
GetEnv               = True
Executable           = $condor_sh/${outName}.sh
Output               = $condor_out/${outName}.out
Error                = $condor_err/${outName}.err
Log                  = $condor_log/${outName}.log

# File transfer behavior
ShouldTransferFiles  = Yes
WhenToTransferOutput = ON_EXIT

# Run job once
queue 1

EOF
#--------Creat cmd file end------------------------------

#--------Creat .sh file ------------------------------
if test -s $condor_sh/${outName}.sh; then
#echo "Execution File is not Empty"
rm  $condor_sh/${outName}.sh
fi
cat >> $condor_sh/${outName}.sh <<EOF
#!/bin/bash
source /cvmfs/dampe.cern.ch/centos7/etc/setup.sh
source /cvmfs/dampe.cern.ch/centos7/opt/DMPSW/trunk/bin/thisdmpsw.sh
source /cvmfs/dampe.cern.ch/centos7/opt/root-6-22-06_python3/bin/thisroot.sh
source /storage/gpfs_data/dampe/users/libowu/software/CADMesh/bashrc_cadmesh

cd /storage/gpfs_data/dampe/users/libowu/CrystalEye/CrystalEyeCode/BkgSimulationCode/build/
./CrystalEye ${Macros} ${outName} ${fSeed}
EOF
#--------Creat .sh file end------------------------------

  ((n++))
	chmod 777 $condor_sh/${outName}.sh
	chmod 777 $condor_cmd/${outName}.cmd

	condor_submit -spool -name sn-02.cr.cnaf.infn.it ${condor_cmd}/${outName}.cmd 

	echo "Submit job_${n}: ${outName} with seed of ${fSeed}"

	done #Loop: jobs
	done #Loop: theta
	done #Loop: Phi

