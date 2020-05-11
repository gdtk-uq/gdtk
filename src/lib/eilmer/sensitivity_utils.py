"""
sensitivity_utils.py -- Small untility functions needed by the
               sensitivity  programs.

.. Author: Luke Doherty (luke.doherty@eng.ox.ac.uk)
           Oxford Thermofluids Institute
           The University of Oxford

  Version: 04-May-2020
"""
#--------------------------------------------------------------------------
def get_values(dict, propertyList):
    """
    Adapted from Peter Jacob's function "get_fractions" which
    may be found within "cea2_gas.py".
    """
    valueList = []
    #print propertyList
    for s in propertyList:
        #print s
        if s in dict.keys():
            #print s
            #print dict[s]
            valueList.append(dict[s])
        else:
            valueList.append(0.0)
            print "WARNING: "+s+"was not found in the current case dictionary."
    return valueList


def write_sensitivity_summary(sensitivity, perturbedVariables, \
                              exitVar, type):
    """
    Write out a file summarising the sensitivity of each exit flow 
    parameter to each of the input parameters.
    """
     
    titleFormatDict = {k : '{0:{fill}>'+"{0}".format(max(11,len(k))+2)+"}" \
                           for k in perturbedVariables}
    #print 'titleFormatDict: ',titleFormatDict
    
    formatDict = {k : '{0:>'+"{0}".format(max(11,len(k))+2)+".5g}" \
                      for k in perturbedVariables} 
    #print 'formatDict: ',formatDict
    
    # Old nenzfr-specific header formatting dictionaries
    #titleFormatDict = {'p1':'{0:{fill}>13}', 'T1':'{0:{fill}>13}',
    #                   'Vs':'{0:{fill}>13}', 'pe':'{0:{fill}>13}',
    #                   'Tw':'{0:{fill}>13}', 'BLTrans':'{0:{fill}>13}',
    #                   'TurbVisRatio':'{0:{fill}>14}',
    #                   'TurbInten':'{0:{fill}>13}',
    #                   'CoreRadiusFraction':'{0:{fill}>20}'}
    #formatDict = {'p1':'{0:13.5g}', 'T1':'{0:>13.5g}',
    #              'Vs':'{0:>13.5g}', 'pe':'{0:>13.5g}',
    #              'Tw':'{0:>13.5g}', 'BLTrans':'{0:>13.5g}',
    #              'TurbVisRatio':'{0:>14.5g}',
    #              'TurbInten':'{0:>13.5g}',
    #              'CoreRadiusFraction':'{0:>20.5g}'}
    
    # Open the file
    if type in ['relative']:
        #fout = open('sensitivities_rel.dat','w')
        fout = open('sensitivities_rel.dat','w')
    elif type in ['absolute']:
        fout = open('sensitivities_abs.dat','w')
        

    # Write header information
    #fout.write('{0:}\n'.format(type+' sensitivities'))
    fout.write('{0:}\n'.format(type+' sensitivities'))
    #fout.write('{0:>13}'.format('variable'))
    fout.write('{0:>13}'.format('variable'))
    #
    for k in perturbedVariables:
        #fout.write(titleFormatDict[k].format(k,fill=''))
        fout.write(titleFormatDict[k].format(k,fill=''))
    fout.write('\n')
    #
    for k in perturbedVariables:
        #fout.write(titleFormatDict[k].format('-',fill='-'))
        fout.write(titleFormatDict[k].format('-',fill='-'))
    fout.write('{0:->13}'.format('-'))
    fout.write('\n')
    
    # Now write out all the data
    for k in exitVar:
        #fout.write('{0:>13}'.format(k))
        fout.write('{0:>13}'.format(k))
        for kk in perturbedVariables:
            #X_kk = sensitivity[kk][exitVar.index(k)]
            X_kk = sensitivity[kk][k]

            fout.write(formatDict[kk].format(X_kk))
        #
        fout.write('\n')
    #
    fout.close()
    
    if type in ['relative']:
       print 'File "sensitivities_rel.dat" written'
    elif type in ['absolute']:
       print 'File "sensitivities_abs.dat" written'

    return    
    
def write_uncertainty_summary(uncertainty, perturbedVariables, \
                              exitVar, inputUncertainties):
    """
    Write out a file summarising the sensitivity of each exit flow
    parameter to each of the input parameters.
    """

    titleFormatDict = {k : '{0:{fill}>'+"{0}".format(max(8,len(k))+2)+"}" \
                           for k in perturbedVariables}
    #print 'titleFormatDict: ',titleFormatDict
    
    formatDict = {k : '{0:>'+"{0}".format(max(8,len(k))+2)+".2%}" \
                      for k in perturbedVariables}
    #print 'formatDict: ',formatDict

    # Old nenzfr-specific header formatting dictionaries
    #titleFormatDict = {'p1':'{0:{fill}>10}', 'T1':'{0:{fill}>10}',
    #                   'Vs':'{0:{fill}>10}', 'pe':'{0:{fill}>10}',
    #                   'Tw':'{0:{fill}>10}', 'BLTrans':'{0:{fill}>10}',
    #                   'TurbVisRatio':'{0:{fill}>14}',
    #                   'TurbInten':'{0:{fill}>11}',
    #                   'CoreRadiusFraction':'{0:{fill}>20}'}
    #formatDict = {'p1':'{0:10.2%}', 'T1':'{0:>10.2%}',
    #              'Vs':'{0:>10.2%}', 'pe':'{0:>10.2%}',
    #              'Tw':'{0:>10.2%}', 'BLTrans':'{0:>10.2%}',
    #              'TurbVisRatio':'{0:>14.2%}',
    #              'TurbInten':'{0:>11.2%}',
    #              'CoreRadiusFraction':'{0:>20.2%}'}

    fout = open('uncertainties.dat','w')
    
    # Write header information
    fout.write('{0:>19}'.format('variable'))
    for k in perturbedVariables:
        fout.write(titleFormatDict[k].format(k,fill=''))
    fout.write('{0:>10}'.format('total'))
    fout.write('\n')
    
    # Write out the uncertainty for each perturbed variable
    fout.write('{0:>19}'.format('(input uncertainty)'))
    for k in perturbedVariables:
        fout.write(formatDict[k].format(inputUncertainties[k]))
    fout.write('{0:>10}'.format(''))
    fout.write('\n')
    for k in perturbedVariables:
        fout.write(titleFormatDict[k].format('-',fill='-'))
    fout.write('{0:->29}'.format('-'))
    fout.write('\n')    
    
    # Now write out all the data. We also calculate and write out the 
    # total uncertainty in each exit flow variable due to the contributions
    # from all input values.
    for k in exitVar:
        fout.write('{0:>19}'.format(k))
        X_total = 0.0
        for kk in perturbedVariables:
            #X_kk = uncertainty[kk][exitVar.index(k)]
            X_kk = uncertainty[kk][k]
            fout.write(formatDict[kk].format(X_kk))
            X_total += X_kk**2 # Calculate total
        fout.write('{0:10.2%}'.format(X_total**0.5))
        fout.write('\n')

    fout.close()

    print 'File "uncertainties.dat" written'
    return


