#!/usr/bin/env bash
set -e

function quote()
{
    A=("${@/#/\"}")
    A=("${A[@]/%/\",}")
    echo "${A[@]}"
}

linkList=$(quote $(mpicc --showme:link))
linkListClean=$(echo $linkList | sed 's/-Wl,//g' | sed 's/-pthread/-lpthread/g')
IFS='%'
configStr="\
            \"name\": \"with-libs\",\n\
            \"dflags-gdc\": [$linkList],\n\
            \"lflags-dmd\": [$linkListClean],\n\
            \"lflags-ldc\": [$linkListClean]"

echo $configStr

if [ ! -f dub.json.back ]
then
    cp dub.json dub.json.back
fi

FLAG=0
while read line; do
    if [ $FLAG == 1 ]
    then
        if [[ $line == *error-noconfig* ]]
        then
            FLAG=0
            continue
        else
            continue
        fi
    fi
    if [[ $line == *error-noconfig* ]]
    then
        FLAG=1
        printf $configStr'\n' >> dub.json.part
        continue
    fi
    echo $line >> dub.json.part
done < dub.json.back

mv dub.json.part dub.json

unset IFS

dub build mpi:configure
dub build mpi:splice
./mpi_splice ./source/mpi/package.d.in <(./mpi_configure $(./gen/get_mpi.h.sh)) > ./source/mpi/package.d
