#! /usr/bin/env python3
"""
Calcualtes forces and moments acting on a group.

Author: IJ 17/4/2020
"""

from getopt import getopt
import numpy as np
import os
import sys


def getforces(fileName_abs, axis_location=[0., 0.]):
    """
    Extract inviscid and viscous forces and moments from a single file in /loads.

    Inputs:
        fileName_abs  - filenName, including absolute path
        axis_location - [x, y] Moment will be evaluated about Z-axis going through
                                this point. Defaults to [0., 0.]

    Outputs:
        Inviscid and Viscous force and moment components.
    """
    # data structures
    x = []  # x-cooridnate of face centre
    y = []  # y-cooridnate of face centre
    p = []  # pressure (Pa)
    A = []  # surface element area (m^2)
    nx = []  # x-component of normal vector
    ny = []  # y-component of normal vector
    tau_w = []  # wall shear stress (Pa)
    lx = []  # x-direction cosine for shear stress
    ly = []  # y-direction cosine for shear stress
    outsign = []  # interface outsign
    n = 0  # number of surface elements

    # read surface data file
    with open(fileName_abs, 'r') as f:
        for line in f:
            dat = line.split()
            if (dat[0] != "#"):  # skip commented lines
                dat = line.split()
                x.append(float(dat[0]))
                y.append(float(dat[1]))
                p.append(float(dat[11]))
                A.append(float(dat[3]))
                nx.append(float(dat[12]))
                ny.append(float(dat[13]))
                tau_w.append(float(dat[7]))
                lx.append(float(dat[8]))
                ly.append(float(dat[8]))  # TODO: Need to check that this is correct.
                outsign.append(float(dat[19]))
                n += 1

    # compute force along flat plate
    Fxi = 0.0  # N
    Fyi = 0.0  # N
    Mzi = 0.0  # Nm
    Fxv = 0.0  # N
    Fyv = 0.0  # N
    Mzv = 0.0  # Nm
    for i in range(0, n):
        Fxi += p[i]*A[i]*(nx[i]*outsign[i])      # inviscid contribution
        Fyi += p[i]*A[i]*(ny[i]*outsign[i])      # inviscid contribution
        Mzi += - p[i]*A[i]*(nx[i]*outsign[i]) * (y[i]-axis_location[1]) \
            + p[i]*A[i]*(ny[i]*outsign[i]) * (x[i]-axis_location[0])
        Fxv += tau_w[i]*A[i]*(lx[i]*outsign[i])  # viscous contribution
        Fyv += tau_w[i]*A[i]*(ly[i]*outsign[i])  # viscous contribution
        Mzv = - tau_w[i]*A[i]*(lx[i]*outsign[i]) * (y[i]-axis_location[1]) \
            + tau_w[i]*A[i]*(ly[i]*outsign[i]) * (x[i]-axis_location[0])
    return Fxi, Fyi, Mzi, Fxv, Fyv, Mzv


def generate_timeList(working_dir_abs, jobfile, uodict, verbosity):
    """Generate the timeList for extracting data."""
    # read time tindx-plot
    if "--tindx-plot" not in uodict:
        time_flag = 'last'
    else:
        if uodict["--tindx-plot"] == "all":
            time_flag = 'all'
        elif uodict["--tindx-plot"] == "last":
            time_flag = 'last'
        else:
            try:
                time_index = int(uodict["--tindx-plot"])
            except:
                raise CustomError("Value entered for --tindx-plot={} not usable.".format(uodict["--tindx-plot"]))
            time_flag = 'index'

    # read file and generate list.
    timesFile = os.path.join(working_dir_abs, "loads", "{}-loads.times".format(jobfile))
    timeList = []
    with open(timesFile, 'r') as f:
        for line in f:
            if line[0] == '#':
                pass
            else:
                temp = line.split(" ")
                i = int(temp[0])
                t = float(temp[1])
                if time_flag == 'all':
                    timeList.append([i, t])
                elif time_flag == 'index':
                    if i == time_index:
                        timeList.append([i, t])
                        break
                    else:
                        pass
                elif time_flag == 'last':
                    if len(timeList) == 0:
                        timeList.append([i, t])
                    else:
                        timeList[-1] = [i, t]
    timeList = np.array(timeList)

    return timeList


def main(uodict):
    """Main Code."""
    if "--verbosity" not in uodict:
        verbosity = 0
    else:
        verbosity = int(uodict["--verbosity"])

    jobfile = uodict["--job"]
    if verbosity > 0:
        print("jobfileName = {}".format(jobfile))

    # get working directory
    if "--absolute-path" in uodict:
        path = uodict["--absolute-path"]
        if path[0] == "~" or not path[0] == os.sep:
            raise CustomError("Path set for '--absolutepath' should be of format '/cwd' .")
        working_dir_abs = os.path.abspath(path)
    elif "--relative-path" in uodict:
        relative_dir = uodict["--relative-path"].split(os.sep)
        working_dir_abs = os.path.abspath(os.path.join(*relative_dir))
    else:
        working_dir_abs = os.path.dirname(os.path.realpath(__file__))

    # generate timeList
    timeList = generate_timeList(working_dir_abs, jobfile, uodict, verbosity)
    if verbosity > 0:
        print("")
        print("timeList has been generated for processing:", timeList)

    # Define set axis loaction from uodict
    if "--axis_location" in uodict:
        raise CustomError("Option '--axis-location' not yet implemented.")
    # Write data to file
    axis_location = [0., 0.]

    data = []
    for i, (tindx, time) in enumerate(timeList):
        if verbosity > 0:
            print("")
            print("Processing forces/moments at tindx={}, time={}[s]".format(tindx, time))
        path = os.path.join(working_dir_abs, "loads", "t{:04d}".format(int(tindx)))
        temp = os.listdir(path)

        # only keep correct group tags
        fileList = []
        if "--compute-loads-on-group" in uodict:
            tag = uodict["--compute-loads-on-group"]
            for file in temp:
                if tag in file:
                    fileList.append(file)
        else:  # if no group tag specified sum across all tags
            fileList = temp
        fileList = sorted(fileList)
        if verbosity > 0:
            print("fileList for processing:", fileList)
        fileList_abs = [os.path.join(path, file) for file in fileList]

        Forces = np.zeros((len(fileList_abs), 6))
        for b, fileName_abs in enumerate(fileList_abs):
            Fxi, Fyi, Mzi, Fxv, Fyv, Mzv = getforces(fileName_abs, axis_location=axis_location)
            Forces[b, :] = [Fxi, Fyi, Mzi, Fxv, Fyv, Mzv]
        Force = np.sum(Forces, axis=0)
        if verbosity > 0:
            # print("Forces:", Forces)
            print("Force:", Force)
        data.append([i, tindx, time, Force])

    if verbosity > 0 or "--output-file" not in uodict:  # print data to screen.
        print("")
        print("Extracted Forces and Moments:")
        for k, d in enumerate(data):
            if k % 20 == 0:
                print(' {:^6}  {:^5}  {:^12}  {:^12}  {:^12}  {:^12}  {:^12}  {:^12}  {:^12} '
                      .format("Index", "tindx", "Time[s]",  "Fxi[N]", "Fyi[N]", "Mzi[Nm]",
                              "Fxv[N]", "Fyv[N]", "Mzv[Nm]"))
            print(' {:>6d}  t{:04d}  {:>12.3e} {:>12.3e} {:>12.3e} {:>12.3e} {:>12.3e} {:>12.3e} {:>12.3e} '
                  .format(d[0], int(d[1]), d[2], d[3][0], d[3][1], d[3][2], d[3][3], d[3][4], d[3][5]))

    if "--output-file" in uodict:
        outFileName = uodict["--output-file"]
        if verbosity > 0:
            print("")
            print("Writing data to '{}' .".format(outFileName))
        outFileName_abs = os.path.join(working_dir_abs, outFileName)
        with open(outFileName_abs, 'w+') as f:
            # Write Comments
            if "--compute-loads-on-group" in uodict:
                tag = uodict["--compute-loads-on-group"]
            else:
                tag = "all"
            f.write("# Group-tag= {}\n".format(str(tag)))
            f.write("# jobfile= {}\n".format(jobfile))
            f.write("# Case-directory= {}\n".format(working_dir_abs))
            # Write Header
            f.write("pos0:Time[s] pos1:Fxi[N] pos2:Fyi[N] pos3:Mzi[Nm] pos4:Fxv[N] pos5:Fyv[N] pos6:Mzv[Nm]\n")
            for d in data:
                f.write("{0:.8e} {1:.8e} {2:.8e} {3:.8e} {4:.8e} {5:.8e} {6:.8e}\n"
                        .format(d[2], d[3][0], d[3][1], d[3][2], d[3][3], d[3][4], d[3][5]))
    return 0


class CustomError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def check_inputs(uo_dict):
    """Check all mandatory options have been provided."""
    reqd_options = ["--job"]
    for op in reqd_options:
        if op not in uo_dict:
            raise CustomError("".join(("compute_forces_and_moments.py requires argument '",
                              op, "', but this argument was not provided.\n")))

    # Check that jobfile exists
    # if not os.path.isfile(uo_dict["--job"]):
    #    raise CustomError("".join(("Jobfile '", uo_dict["--job"], "' not found,",
    #                      " check that file path is correct.\n")))


def print_usage():
    """Show usage instructions."""
    print("")
    print("Function to extract total forces and moments  ")
    print("Argument:                            Comment:")
    print("---------------------------------------------------------------------")
    print("  --job=<string>                     file names built from this string")
    print("  --tindx-plot=<int>|all|last        defaults to last")
    print("  --output-file=<string>             defaults to stdout")
    print('  --compute-loads-on-group=""        group tag')
    print('                                     (if not specified, tags will be summed)')
    print('  --relative-path="."                relative path to Case directory')
    print('  --absolute-path=" "                absolute path to Case director, takes precedence over relative path')
    print('  --axis_location="x,y"              location for Z-axis (default [0., 0.])' )
    print("  --verbosity                        0, 1, 2 - sets level of output")
    print("")


short_options = ""
long_options = ["help", "job=", "tindx-plot=", "output-file=", "compute-loads-on-group=", "relative-path=", "absolute-path=", "verbosity="]


if __name__ == '__main__':
    user_options = getopt(sys.argv[1:], short_options, long_options)
    uo_dict = dict(user_options[0])

    if len(user_options[0]) == 0 or "--help" in uo_dict:
        print_usage()
        sys.exit(1)

    else:
        check_inputs(uo_dict)
        main(uo_dict)
        print("\n")
        print("SUCCESS.")
        print("\n")
