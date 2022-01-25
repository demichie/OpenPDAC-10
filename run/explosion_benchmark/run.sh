foamCleanCase
blockMesh 
checkMesh -allTopology -allGeometry
changeDictionary
cp 0/alpha.air.orig 0/alpha.air
cp 0/alpha.particles.orig 0/alpha.particles
cp 0/T.air.orig 0/T.air
cp 0/T.particles.orig 0/T.particles
cp 0/ph_rgh.orig 0/ph_rgh
cp 0/p.orig 0/p
cp 0/rho.orig 0/rho

cp system/controlDict_init system/controlDict
cp system/fvSolution_init system/fvSolution

#FOR PARALLEL RUN:
#decomposePar
#mpirun -np 8 myMultiphaseEulerFoamNew -parallel
#reconstructPar -newTimes

#FOR SCALAR RUN:
myMultiphaseEulerFoamNew

setFields
cp system/controlDict_run system/controlDict
cp system/fvSolution_run system/fvSolution

#FOR PARALLEL RUN:
#decomposePar
#mpirun -np 8 myMultiphaseEulerFoamNew -parallel
#reconstructPar -newTimes

#FOR SCALAR RUN:
myMultiphaseEulerFoamNew





