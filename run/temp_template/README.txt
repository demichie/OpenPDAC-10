foamCleanCase
blockMesh 
snappyHexMesh -overwrite
extrudeMesh
checkMesh -allTopology -allGeometry
changeDictionary 
cp 0/T.air.orig 0/T.air
cp 0/T.particles.orig 0/T.particles
cp 0/ph_rgh.orig 0/ph_rgh
cp 0/p.orig 0/p
cp 0/rho.orig 0/rho

FOR PARALLEL RUN:
decomposePar
mpirun -np 8 OpenPDAC -parallel
reconstructPar -newTimes

FOR SCALAR RUN:
OpenPDAC



