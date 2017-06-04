[![Build Status](https://travis-ci.org/DlangScience/OpenMPI.svg?branch=master)](https://travis-ci.org/DlangScience/OpenMPI)
# OpenMPI
D bindings to OpenMPI.

How to use:
==========
This is a slightly unusual dub package, in that it requires some extra steps to use. The first time you try and use this package as a dependency (unless you specify the "noLibs" configuration) you will get an error message saying that it has not been configured. Follow the instructions and hopefully everything will just work. Essentially all that is necessary is to run ```gen/setup.sh```.

Here's what to do if it doesn't work out:

Option 1:
Manually tweak the dub.json file to get it to work with your OpenMPI system. This probably involves changing the ```"dflags"``` and ```"lflags"``` fields. Instead of editing the file in the ```dub``` cache (e.g. ```~/.dub/packages/```), you might want to make your own copy to make it more permanant. Using ```dub add-local``` or ```dub add-path``` will allow you to locally register this copy with ```dub``` so it works automatically for all your projects.

Option 2:
Use the "noLibs" configuration and work out how to link to openmpi yourself (```mpicc --showme:link``` will show you what is necessary)

Status:
=======
At the moment everything should work fine on OpenMPI versions 1.4.3 upwards, not including pre-release versions. Other versions may work, but haven't been tested. However, please do report any problems you have, on any version of OpenMPI and any system, however old or new.
