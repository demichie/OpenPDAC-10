module load openFOAM-10
source $FOAM_BASHRC

foamCleanCase

cd preprocessing
python createSTL.py
cd ..


blockMesh 
checkMesh -allTopology -allGeometry

snappyHexMesh -overwrite
extrudeMesh
changeDictionary

# sbatch MPIJob_snappy.script
# extrudeMesh
# changeDictionary

rm -rf 0
cp -r org.0 0
cp 0/alpha.air.init 0/alpha.air
cp 0/alpha.particles.init 0/alpha.particles
cp 0/T.air.init 0/T.air
cp 0/T.particles.init 0/T.particles
cp 0/U.air.init 0/U.air
cp 0/U.particles.init 0/U.particles


cp ./system/controlDict.init ./system/controlDict
cp ./system/fvSolution.init ./system/fvSolution

#FOR PARALLEL RUN:
#sbatch MPIJob_init.script
#squeue

#FOR SCALAR RUN:
OpenPDAC


cp 0/alpha.air.run 0/alpha.air
cp 0/alpha.particles.run 0/alpha.particles
cp 0/T.air.run 0/T.air
cp 0/T.particles.run 0/T.particles
cp 0/U.air.run 0/U.air
cp 0/U.particles.run 0/U.particles

cp ./system/controlDict.run system/controlDict
cp ./system/fvSolution.run system/fvSolution

#FOR PARALLEL RUN:
#sbatch MPIJob_run.script
#squeue

#FOR SCALAR RUN:
OpenPDAC

