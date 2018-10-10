#Runs the selected shock tube example in MHDShockTube.lua

#Clean the directory of old data
rm -r flow/*
rm -r grid/*
rm -r hist/*
rm -r plot/*
rm -r solid/*
rm -r loads/*

#Create the gas file
prep-gas ideal-air.inp ideal-air-gas-model.lua

#Run the analytical version- this is only useful for Sod's test
python ShockTubeCalc.py

#Execute the eilmer file with custom post-processing

e4shared --prep --job=MHDShockTube
if [ "$?" -ne "0" ] ; then
    echo "e4prep-clean ended abnormally."
fi

e4shared --run --job=MHDShockTube --verbosity=1 --max-cpus=2
if [ "$?" -ne "0" ] ; then
    echo "e4main-clean ended abnormally."
fi

e4shared --post --job=MHDShockTube --tindx-plot=all --vtk-xml

e4shared --custom-post --script-file=extract_tfinal.lua

#This is currently only applied to Sod's test, but with some adjustments could be applied to the MHD examples
python error_check.py 

gnuplot ShockTubePlot.gnuplot
