#!/usr/bin/tcsh
#PBS -l walltime=15:00:00
#PBS -l mem=2gb

# copy files on the temporary directory from your directory
#cp ../auxiliary_files
# ...
# cp $DIR/infilen .

# execution of hte program

./playground

# copy of the results in a proper directory
cp output.txt Output
# ...
# cp outfilen $DIR/.


