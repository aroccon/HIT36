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
echo "=                                   Local UNIX                               ="
fi
######################################################################################
#sim name, at least 10 characters
SIM_NAME="nameofthesim"
SPLIT="0" #fluid and particles split, to be removed from command line
# define machine
# 0 : local (Mac)
cd set_run
#create the folder (if missing)
mkdir -p sc_compiled
mkdir -p results
mkdir -p paraview_vtk
#clean the folders
rm -r sc_compiled/*
rm -r results/*
rm -r paraview_vtk/*
cp -r ../paraview_vtk/* ./paraview_vtk


cd sc_compiled
cp -r ../../source_code/* ./
rm hit36
rm *.o
rm *.mod
make clean
make
rm *.o

#clear
#running the code
mpirun -np 4 ./hit36 $SIM_NAME $SPLIT

#!  make sure DT is appropriate for scalars
#!  DT < 0.09 * dx^2*Pe
