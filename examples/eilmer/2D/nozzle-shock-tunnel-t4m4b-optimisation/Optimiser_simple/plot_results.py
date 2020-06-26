#! /usr/bin/env python3
# optimize.py
# Automate the running of the flow simulations while searching
# for the optimum parameters (angles) that define the nozzle shape.
#
# PJ, 2018-03-04, take bits from nenzfr

from optimize import *

if __name__ == "__main__":
    print("Now plotting the results")
    #first want to plot the contours
    # bezCtrlPts_orig = np.loadtxt("Bezier-control-pts-t4m4b-initial.data",skiprows=1)
    # bezCtrlPts_opt = np.loadtxt("Bezier-control-pts-t4m4b.opt.data",skiprows=1)
    # tlist = np.linspace(0,1,100)
    # xlist_orig = []
    # ylist_orig = []
    # xlist_opt = []
    # ylist_opt = []     
    # for t in tlist:
    #     points = eval_Bezier(bezCtrlPts_orig, t)
    #     xlist_orig.append(points[0])
    #     ylist_orig.append(points[1])

    #     points = eval_Bezier(bezCtrlPts_opt, t)
    #     xlist_opt.append(points[0])
    #     ylist_opt.append(points[1])     

    # plt.figure()
    # plt.plot(xlist_opt, ylist_opt,label="opt")
    # plt.plot(xlist_orig, ylist_orig,label="orig")  
    # plt.legend()

    #Now want to plot the outflows
    original_data = np.loadtxt("{0}-exit-initial.data".format(jobname),skiprows=1)
    optimized_data = np.loadtxt("{0}-exit.data".format(jobname),skiprows=1)

    M_opt = optimized_data[:,18]
    y_opt = optimized_data[:,1]

    M_orig = original_data[:,18]
    y_orig = original_data[:,1]

    plt.figure()
    plt.plot(y_opt, M_opt,label="opt")
    plt.plot(y_orig, M_orig,label="orig")
    plt.legend()
    plt.show()     
