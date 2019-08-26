import numpy as np
import pandas as pd
import Nucleaseq_data_processing as processing

def combined_average_data(rep_on,file_clv,path_on,path_clv):
    
    xdata, ydata, yerr = processing.prepare_multiprocessing_combined(rep_on,file_clv,path_on,path_clv,True,False,True,True)
    
    xdata_average = []
    ydata_average = []
    yerr_average = []
    
    
    kclv = []
    errclv = []
    kon = []
    erron = []
    i = 1
    while i < 20+1:
        for a in range(len(xdata)):
            if len(xdata[a])==1 and xdata[a][0]==i:
                kclv.append(ydata[a][0][0])
                errclv.append(yerr[a][0][0])
                kon.append((ydata[a][1][0]))
                erron.append((yerr[a][1][0]))
                i = i + 1

    kclv.append(ydata[0][0][0])
    kclv = np.array(kclv)
    errclv.append(yerr[0][0][0])
    errclv = np.array(errclv)
    kon.append((ydata[0][1][0]))
    kon = np.array(kon)
    erron.append(yerr[0][1][0])
    erron = np.array(erron)


    k_double = np.zeros([20,20])
    k_err = np.zeros([20,20])
    on_double = np.zeros([20,20])
    on_err = np.zeros([20,20])

    for a in range(len(xdata)):
        if len(xdata[a])==2:
            i = xdata[a][0]
            j = xdata[a][1]
            k_double[j-1,i-1] = (ydata[a][0][0])
            k_err[j-1,i-1] = yerr[a][0][0]
            if len(ydata[a][1])>0:
                on_double[j-1,i-1] = (ydata[a][1][0])
                on_err[j-1,i-1] = yerr[a][1][0]
            else:
                on_double[j-1,i-1] = np.nan
                on_err[j-1,i-1] = np.nan
    
    #OT
    xdata_average.append(xdata[0])
    ydata_average.append(ydata[0])
    yerr_average.append(yerr[0])
    
    #smm
    #1-7
    xdata_average.append([1])
    weightsclv, errorclv = processing.weighting(errclv[0:7])
    weightson, erroron = processing.weighting(np.append(erron[0],erron[2:5])) #eliminating 2,6,7
    ydata_average.append([[np.average(kclv[0:7],weights=weightsclv)],[np.average(np.append(kon[0],kon[2:5]),weights=weightson)]])
    yerr_average.append([[errorclv],[erroron]])

    #8-11
    xdata_average.append([8])
    weightsclv, errorclv = processing.weighting(errclv[7:11])
    weightson, erroron = processing.weighting(erron[7:11])
    ydata_average.append([[np.average(kclv[7:11],weights=weightsclv)],[np.average(kon[7:11],weights=weightson)]])
    yerr_average.append([[errorclv],[erroron]])

    #12-17
    xdata_average.append([12])
    weightsclv, errorclv = processing.weighting(errclv[11:17])
    weightson, erroron = processing.weighting(erron[11:17])
    ydata_average.append([[np.average(kclv[11:17],weights=weightsclv)],[np.average(kon[11:17],weights=weightson)]])
    yerr_average.append([[errorclv],[erroron]])

    #18-20
    xdata_average.append([18])
    weightsclv, errorclv = processing.weighting(errclv[17:20])
    weightson, erroron = processing.weighting(erron[17:20])
    ydata_average.append([[np.average(kclv[17:20],weights=weightsclv)],[np.average(kon[17:20],weights=weightson)]])
    yerr_average.append([[errorclv],[erroron]])

    #dmm, numbered 1-10
    #1
    xdata_average.append([1,2])
    
    tempclv = []
    temperrclv = []
    for j in range(0,7):
        for i in range(0,j):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(2,5):
        for i in range(0,1)+range(2,j):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #2
    xdata_average.append([1,8])
    
    tempclv = []
    temperrclv = []
    for j in range(7,11):
        for i in range(0,7):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(7,11):
        for i in range(0,1)+range(2,5):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #3
    xdata_average.append([1,12])
    
    tempclv = []
    temperrclv = []
    for j in range(11,17):
        for i in range(0,7):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(11,17):
        for i in range(0,1)+range(2,5):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #4
    xdata_average.append([1,18])
    
    tempclv = []
    temperrclv = []
    for j in range(17,20):
        for i in range(0,7):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(17,20):
        for i in range(0,1)+range(2,5):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #5
    xdata_average.append([8,9])
    
    tempclv = []
    temperrclv = []
    for j in range(7,11):
        for i in range(7,j):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(7,11):
        for i in range(7,j):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #6
    xdata_average.append([8,12])
    
    tempclv = []
    temperrclv = []
    for j in range(11,17):
        for i in range(7,11):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(11,17):
        for i in range(7,11):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #7
    xdata_average.append([8,18])
    
    tempclv = []
    temperrclv = []
    for j in range(17,20):
        for i in range(7,11):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(17,20):
        for i in range(7,11):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #8
    xdata_average.append([12,13])
    
    tempclv = []
    temperrclv = []
    for j in range(11,17):
        for i in range(11,j):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(11,17):
        for i in range(11,j):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #9
    xdata_average.append([12,18])
    
    tempclv = []
    temperrclv = []
    for j in range(17,20):
        for i in range(11,17):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(17,20):
        for i in range(11,17):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #10
    xdata_average.append([18,19])
    
    tempclv = []
    temperrclv = []
    for j in range(17,20):
        for i in range(17,j):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(17,20):
        for i in range(17,j):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    
    
    return xdata_average, ydata_average, yerr_average

def combined_average_data_aba(file_aba,file_clv,path_aba,path_clv):
    
    xdata, ydata, yerr = processing.prepare_multiprocessing_combined_aba(file_aba,file_clv,path_aba,path_clv,True)
    
    xdata_average = []
    ydata_average = []
    yerr_average = []
    
    
    kclv = []
    errclv = []
    kon = []
    erron = []
    i = 1
    while i < 20+1:
        for a in range(len(xdata)):
            if len(xdata[a])==1 and xdata[a][0]==i:
                kclv.append(ydata[a][0][0])
                errclv.append(yerr[a][0][0])
                kon.append((ydata[a][1][0]))
                erron.append((yerr[a][1][0]))
                i = i + 1

    kclv.append(ydata[0][0][0])
    kclv = np.array(kclv)
    errclv.append(yerr[0][0][0])
    errclv = np.array(errclv)
    kon.append((ydata[0][1][0]))
    kon = np.array(kon)
    erron.append(yerr[0][1][0])
    erron = np.array(erron)


    k_double = np.zeros([20,20])
    k_err = np.zeros([20,20])
    on_double = np.zeros([20,20])
    on_err = np.zeros([20,20])

    for a in range(len(xdata)):
        if len(xdata[a])==2:
            i = xdata[a][0]
            j = xdata[a][1]
            k_double[j-1,i-1] = (ydata[a][0][0])
            k_err[j-1,i-1] = yerr[a][0][0]
            if len(ydata[a][1])>0:
                on_double[j-1,i-1] = (ydata[a][1][0])
                on_err[j-1,i-1] = yerr[a][1][0]
            else:
                on_double[j-1,i-1] = np.nan
                on_err[j-1,i-1] = np.nan
    
    #OT
    xdata_average.append(xdata[0])
    ydata_average.append(ydata[0])
    yerr_average.append(yerr[0])
    
    #smm
    #1-7
    xdata_average.append([1])
    weightsclv, errorclv = processing.weighting(errclv[0:7])
    weightson, erroron = processing.weighting(np.append(erron[0],erron[2:5])) #eliminating 2,6,7
    ydata_average.append([[np.average(kclv[0:7],weights=weightsclv)],[np.average(np.append(kon[0],kon[2:5]),weights=weightson)]])
    yerr_average.append([[errorclv],[erroron]])

    #8-11
    xdata_average.append([8])
    weightsclv, errorclv = processing.weighting(errclv[7:11])
    weightson, erroron = processing.weighting(erron[7:11])
    ydata_average.append([[np.average(kclv[7:11],weights=weightsclv)],[np.average(kon[7:11],weights=weightson)]])
    yerr_average.append([[errorclv],[erroron]])

    #12-17
    xdata_average.append([12])
    weightsclv, errorclv = processing.weighting(errclv[11:17])
    weightson, erroron = processing.weighting(erron[11:17])
    ydata_average.append([[np.average(kclv[11:17],weights=weightsclv)],[np.average(kon[11:17],weights=weightson)]])
    yerr_average.append([[errorclv],[erroron]])

    #18-20
    xdata_average.append([18])
    weightsclv, errorclv = processing.weighting(errclv[17:20])
    weightson, erroron = processing.weighting(erron[17:20])
    ydata_average.append([[np.average(kclv[17:20],weights=weightsclv)],[np.average(kon[17:20],weights=weightson)]])
    yerr_average.append([[errorclv],[erroron]])

    #dmm, numbered 1-10
    #1
    xdata_average.append([1,2])
    
    tempclv = []
    temperrclv = []
    for j in range(0,7):
        for i in range(0,j):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(2,5):
        for i in range(0,1)+range(2,j):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #2
    xdata_average.append([1,8])
    
    tempclv = []
    temperrclv = []
    for j in range(7,11):
        for i in range(0,7):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(7,11):
        for i in range(0,1)+range(2,5):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #3
    xdata_average.append([1,12])
    
    tempclv = []
    temperrclv = []
    for j in range(11,17):
        for i in range(0,7):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(11,17):
        for i in range(0,1)+range(2,5):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #4
    xdata_average.append([1,18])
    
    tempclv = []
    temperrclv = []
    for j in range(17,20):
        for i in range(0,7):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(17,20):
        for i in range(0,1)+range(2,5):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #5
    xdata_average.append([8,9])
    
    tempclv = []
    temperrclv = []
    for j in range(7,11):
        for i in range(7,j):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(7,11):
        for i in range(7,j):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #6
    xdata_average.append([8,12])
    
    tempclv = []
    temperrclv = []
    for j in range(11,17):
        for i in range(7,11):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(11,17):
        for i in range(7,11):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #7
    xdata_average.append([8,18])
    
    tempclv = []
    temperrclv = []
    for j in range(17,20):
        for i in range(7,11):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(17,20):
        for i in range(7,11):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #8
    xdata_average.append([12,13])
    
    tempclv = []
    temperrclv = []
    for j in range(11,17):
        for i in range(11,j):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(11,17):
        for i in range(11,j):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #9
    xdata_average.append([12,18])
    
    tempclv = []
    temperrclv = []
    for j in range(17,20):
        for i in range(11,17):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(17,20):
        for i in range(11,17):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #10
    xdata_average.append([18,19])
    
    tempclv = []
    temperrclv = []
    for j in range(17,20):
        for i in range(17,j):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
            
    tempon = []
    temperron = []
    for j in range(17,20):
        for i in range(17,j):
            tempon.append(on_double[j][i])
            temperron.append(on_err[j][i])
    weights,erron = processing.weighting(temperron)
    avgon = np.average(tempon,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    
    
    return xdata_average, ydata_average, yerr_average


def cleavage_average_data(file_clv,path_clv):
    
    xdata, ydata, yerr = processing.prepare_multiprocessing_nucleaseq_log(file_clv,path_clv,True)
    
    xdata_average = []
    ydata_average = []
    yerr_average = []
    
    
    kclv = []
    errclv = []
    i = 1
    while i < 20+1:
        for a in range(len(xdata)):
            if len(xdata[a])==1 and xdata[a][0]==i:
                kclv.append(ydata[a][0])
                errclv.append(yerr[a][0])
                i = i + 1

    kclv.append(ydata[0][0])
    kclv = np.array(kclv)
    errclv.append(yerr[0][0])
    errclv = np.array(errclv)


    k_double = np.zeros([20,20])
    k_err = np.zeros([20,20])

    for a in range(len(xdata)):
        if len(xdata[a])==2:
            i = xdata[a][0]
            j = xdata[a][1]
            k_double[j-1,i-1] = (ydata[a][0])
            k_err[j-1,i-1] = yerr[a][0]
            
    #OT
    xdata_average.append(xdata[0])
    ydata_average.append(ydata[0])
    yerr_average.append(yerr[0])
    
    #smm
    #1-7
    xdata_average.append([1])
    weightsclv, errorclv = processing.weighting(errclv[0:7])
    ydata_average.append([np.average(kclv[0:7],weights=weightsclv)])
    yerr_average.append([errorclv])

    #8-11
    xdata_average.append([8])
    weightsclv, errorclv = processing.weighting(errclv[7:11])
    ydata_average.append([np.average(kclv[7:11],weights=weightsclv)])
    yerr_average.append([errorclv])

    #12-17
    xdata_average.append([12])
    weightsclv, errorclv = processing.weighting(errclv[11:17])
    ydata_average.append([np.average(kclv[11:17],weights=weightsclv)])
    yerr_average.append([errorclv])

    #18-20
    xdata_average.append([18])
    weightsclv, errorclv = processing.weighting(errclv[17:20])
    ydata_average.append([np.average(kclv[17:20],weights=weightsclv)])
    yerr_average.append([errorclv])

    #dmm, numbered 1-10
    #1
    xdata_average.append([1,2])
    
    tempclv = []
    temperrclv = []
    for j in range(0,7):
        for i in range(0,j):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #2
    xdata_average.append([1,8])
    
    tempclv = []
    temperrclv = []
    for j in range(7,11):
        for i in range(0,7):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #3
    xdata_average.append([1,12])
    
    tempclv = []
    temperrclv = []
    for j in range(11,17):
        for i in range(0,7):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #4
    xdata_average.append([1,18])
    
    tempclv = []
    temperrclv = []
    for j in range(17,20):
        for i in range(0,7):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #5
    xdata_average.append([8,9])
    
    tempclv = []
    temperrclv = []
    for j in range(7,11):
        for i in range(7,j):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #6
    xdata_average.append([8,12])
    
    tempclv = []
    temperrclv = []
    for j in range(11,17):
        for i in range(7,11):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #7
    xdata_average.append([8,18])
    
    tempclv = []
    temperrclv = []
    for j in range(17,20):
        for i in range(7,11):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #8
    xdata_average.append([12,13])
    
    tempclv = []
    temperrclv = []
    for j in range(11,17):
        for i in range(11,j):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #9
    xdata_average.append([12,18])
    
    tempclv = []
    temperrclv = []
    for j in range(17,20):
        for i in range(11,17):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    #10
    xdata_average.append([18,19])
    
    tempclv = []
    temperrclv = []
    for j in range(17,20):
        for i in range(17,j):
            tempclv.append(k_double[j][i])
            temperrclv.append(k_err[j][i])
    weights,errclv = processing.weighting(temperrclv)
    avgclv = np.average(tempclv,weights=weights)
    
    ydata_average.append([avgclv])
    yerr_average.append([errclv])
    
    
    
    return xdata_average, ydata_average, yerr_average