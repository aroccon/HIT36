SIM_NAME="BENCHMARK1"
#mkdir set_run
# define machine
# 0 : local (Mac)
machine="0"
#echo ""
#echo "=============================================================================="
#echo "=                                 Running on                                 ="
if [ "$machine" == "0" ]; then
echo "=                            Local machine (OS X)                            ="
cp ./os_x/makefile_osx ./source_code/Makefile

elif [ "$machine" == "1" ]; then
echo "=                                   Local UNIX                               ="
fi
######################################################################################
cd set_run
#create the folder (if missing)
mkdir -p sc_compiled
mkdir -p results
#clean the folders
rm -r sc_compiled/*
rm -r results/*

cd sc_compiled
cp -r ../../source_code/* ./
rm *.o
rm *.mod
make clean
make
rm *.o
