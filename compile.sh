# 0 : local (Mac)
# 1 : local (unix)
machine="0"
#echo ""
#echo "=============================================================================="
#echo "=                                 Running on                                 ="
if [ "$machine" == "0" ]; then
echo "=                                                                            ="
echo "=                            Local machine (OS X)                            ="
echo "=                                                                            ="
cp ./os_x/makefile_osx ./source_code/Makefile

elif [ "$machine" == "1" ]; then
echo "=                                                                            ="
echo "=                                   Local UNIX                               ="
echo "=                                                                            ="
cp ./local/makefile_local ./source_code/Makefile

fi
######################################################################################
#### PARAMETERS OF THE RUN
#number of MPI tasks

num_task="2"
#with particles, at least two tasks in parts (parts tasks = numprocs/4)

split="0"
#fluid, stats and particles split of the communicator.
#If there are no particles:
#split=0: All MPI tasks solve the fields and do the stats.
#split=1: 2/3 MPI tasks solve the fields and 1/3 do the stats.
#If there are particles
#split=0: Does not work.
#split=1: 1/2 MPI tasks fields, 1/4 stats and 1/4 LPT
#see m_openmpi.f90 for details

nx="64"
#number of nodes (nx=ny=nz)

first_iteration="0"
#fresh strart?

#maximum time steps number
last_iteration="20000"

#frequency for stats 
dump_statistics="10"

#frequency for restart
dump_restart="10"

#frequency for output
dump_output="100"

#time step (double)
dt="0.0005"

#final time (double)
maxtime="10.0"

#Reynolds number (double)
re="40."

##### END OF PARAMETERS DEFINITION
######################################################################################

# define machine
# 0 : local (Mac)
cd set_run
#create the folder (if missing)
mkdir -p sc_compiled
mkdir -p results
mkdir -p paraview_fields
mkdir -p paraview_particles
#clean the folders
rm -r sc_compiled/*
rm -r results/*
rm -r paraview_fields/*
rm -r paraview_particles/*
cp -r ../paraview_fields/*    ./paraview_fields
cp -r ../paraview_particles/* ./paraview_particles
cd sc_compiled
cp -r ../../source_code/* ./
rm hit36
rm *.o
rm *.mod


#before compiling, setting the split flag
# only for Mac
if [ "$machine" == "0" ]; then
  if [ "$split" == "0" ]; then
    echo "i'm here"
    sed -i "" "s/splitflag/never/g" ./m_openmpi.f90
  fi
  if [ "$split" == "1" ]; then
    sed -i "" "s/splitflag/split/g" ./m_openmpi.f90
  fi
#for the rest of the machines
else
  if [ "$split" == "0" ]; then
    sed -i "s/splitflag/never/g" ./m_openmpi.f90
  fi
  if [ "$split" == "1" ]; then
    sed -i "s/splitflag/split/g" ./m_openmpi.f90
  fi
fi

#edit input paramters in set_run
# MAC OS 
if [ "$machine" == "0" ]; then
  sed -i "" "s/nxnynz/$nx/g" ./input.f90
  sed -i "" "s/firstiter/$first_iteration/g" ./input.f90
  sed -i "" "s/laststep/$last_iteration/g" ./input.f90
  sed -i "" "s/dumpstats/$dump_statistics/g" ./input.f90
  sed -i "" "s/dumprestart/$dump_restart/g" ./input.f90
  sed -i "" "s/dumpoutput/$dump_output/g" ./input.f90
  sed -i "" "s/itmax/$maxtime/g" ./input.f90
  sed -i "" "s/reeeeee/$re/g" ./input.f90
  sed -i "" "s/dtttttt/$dt/g" ./input.f90
else
# Linux 
  sed -i "s/nxnynz/$nx/g" ./input.f90
  sed -i "s/nxnynz/$nx/g" ./input.f90
  sed -i "s/firstiter/$first_iteration/g" ./input.f90
  sed -i "s/laststep/$last_iteration/g" ./input.f90
  sed -i "s/dumpstats/$dump_statistics/g" ./input.f90
  sed -i "s/dumprestart/$dump_restart/g" ./input.f90
  sed -i "s/dumpoutput/$dump_output/g" ./input.f90
  sed -i "s/itmax/$maxtime/g" ./input.f90
  sed -i "s/reeeeee/$re/g" ./input.f90
  sed -i "s/dtttttt/$dt/g" ./input.f90
fi

#compiling
make clean
make
rm *.o

#clear
#running the code
echo "Compilation completed, ready to run (if compilation was ok)"

mpirun -np $num_task ./hit36

#!  make sure DT is appropriate for scalars
#!  DT < 0.09 * dx^2*Pe
