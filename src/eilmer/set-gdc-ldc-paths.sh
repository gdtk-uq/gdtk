# set-gdc-ldc-paths.sh
#
# When experimenting with the alternative D-language compilers,
# it is convenient to place their complete distribution of files
# into a directory within your $HOME directory.
# One of the following commands can be executed manually (cut and paste)
# to set the PATH to point to the appropriate directory and then 
# you'll be ready to compile.
#
# Updated to the packages that are current at 2016-03-29.
#
export PATH=$PATH:$HOME/i686-pc-linux-gnu/bin
export PATH=$PATH:$HOME/x86_64-pc-linux-gnu/bin
export PATH=$PATH:$HOME/ldc2-0.17.1-linux-x86/bin
export PATH=$PATH:$HOME/ldc2-0.17.1-linux-x86_64/bin
