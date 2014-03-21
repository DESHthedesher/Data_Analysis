# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 15:08:12 2013

This program takes a list of force curves,
finds the fractional extension,
corrects for the number of monomers,
and produces a figure of the energy per monomer for each 
data set delivered.

Forces and errors are from the data sheet "all data summary", page 2

@author: Duncan
"""
import matplotlib.pyplot as plt
import numpy as np

Lmax = 0.0

def fractional_extend_calc(force):
    ###This function assumes a WLC extension with the given parameters.
    ###Takes the Equation for force assuming WLC behaviour:
    ### F = KT/L( 1/4*(1-phi)**-2 - 1/4 + phi) (f is force, L is persitence length, phi is extension ratio)
    ### Transformed into the cubic polynomial paramaterized with A = FL/KT + 1/4:
    ### 0 = Phi**3 + Phi**2(-2-A) +Phi(1+2A) - A + 1/4
    ### Finds the roots, and returns the one real root

    force = force*10**-12
    persistence_length = 0.92*10**-9
    temp = 300
    boltzmann = 1.3806488*10**-23
    
    A = (force*persistence_length)/(temp*boltzmann) + 0.25
    
    poly = [1, -2-A, 1 + 2*A, 0.25 - A]
    
    roots = np.roots(poly) 
    
    
    return roots[2]
   
def energy_per_monomer(force, phi):
    ###Makes use of the formula E/N = F*D/N
    ### Works out to be E/N = d/Lc * Lm *F = Phi*Lm*F
    force = force*10**-12        
    monomer_length = 0.24*10**-9    
    lm = monomer_length
    k = 1.3806488*10**-23
    t = 300


    EpM = phi*lm*force/(k*t)
    
    return EpM
    
    
def linreg_waterzero(x,y,ID):
     
    
    sumxy = 0.0
    sumxx = 0.0
    sumx = 0.0

    b = 3.993850

    for i in range(0, len(x)):
        sumxy += x[i]*y[i]
        sumxx += x[i]**2
        sumx += x[i]

    slope = (sumxy - b*sumx)/float(sumxx)
    xfit = np.arange(0, max(x), 0.01)
    yfit = []
    for i in xfit:
        newy = b + slope*i
        yfit.append(newy)
        
    SSres = 0.0
    SStot = 0.0
    ybar = 0.0    
    
    for i in y:
        ybar += i
    ybar = ybar/len(y)
    

    for i in range(0, len(x)):
        SSres += (slope*x[i] + 3.993850 - y[i])**2
        SStot += (ybar - y[i])**2
        
    #NOTE this rsquare does not take error into account    
    rsquare = 1.0 - SSres/SStot        
 
    #print '\n r^2 of',rsquare
    
    print '\n',ID,'Has a slope of', slope
        

    return(xfit, yfit)
   
def linreg_statfit(x,y,e,ID):
    #uses equations 6-12 from bevington for slope paramters
    #equations 6-25 for error parameters
    
    
    sumxys = 0.0
    sumxxs = 0.0
    sumxs = 0.0
    sumys = 0.0
    sums = 0.0    
    
    for i in range(0,len(x)):
        sums += 1.0/float(e[i]**2)
        sumxs += (x[i]/float(e[i]**2))
        sumys += (y[i]/(float(e[i]**2)))
        sumxxs += (x[i]**2)/(float(e[i]**2))
        sumxys += ((y[i])*(x[i]))/(float(e[i]**2))

    delta = sums*sumxxs - sumxs**2
    a = (1/delta)*(sumxxs*sumys - sumxs*sumxys)
    b = (1/delta)*(sums*sumxys - sumxs*sumys)
    
    SSres = 0.0
    SStot = 0.0
    ybar = 0.0    
    
    for i in y:
        ybar += i
    ybar = ybar/len(y)
    

    for i in range(0, len(x)):
        SSres += (b*x[i] + a - y[i])**2
        SStot += (ybar - y[i])**2
        
    #NOTE this rsquare does not take error into account    
    rsquare = 1.0 - SSres/SStot        
 
 
    intercept_error = ((1/delta)*(sumxxs))**0.5
    slope_error = ((sums)/float(delta))**0.5
    
 
 
    print '\n',ID,"Error corrected slope of",b,"and an error corrected intercept of",a
    print 'r^2 of',rsquare
    print "intercept +/-",intercept_error
    print "slope +/-", slope_error
    
    
    
    
    xfit = np.arange(0, max(x), 0.01)
    yfit = []
    for i in xfit:
        yfit.append(b*i+a)
    
    return (xfit, yfit)
    
def Qdefine(x,y,e,delta,k):
    #Gives the error-weighted sum of squares using the langmuir equation    
    
    Q = 0.000
    ybar = 0.00
    pbar = 0.00
    
    for i in range(1, len(x)):
        Q += (1/(2*(e[i]**2)))*(3.993850+(delta*k*x[i])/(1+k*x[i]) - y[i])**2
        pbar += 3.993850+(delta*k*x[i])/(1+k*x[i])
        ybar += y[i]
        
        
    if pbar < ybar:
        sign = 1
    elif pbar >= ybar:
        sign = -1
    
    
    
    return Q,sign


def langguess_statfit(x,y,e,ID):
    #gradient based minimization to the langmuir isotherm
    global Lmax
        
    kold = 0
    knew = 0.01
    step_size = 0.001
    precision = 0.00000001    
        
    avgY = 0.0
    for i in y:
        avgY += i - y[0]
    
    '''if avgY < 0:
        knew = knew*-1
        step_size = step_size*-1'''
    
    
    while abs(knew - kold) > precision:
        kold = knew
        knew = kold + step_size*change_in_error(kold,x,y,e)
        
    

    print "k for %s ="%(ID),kold
    print "Lmax for %s ="%(ID),Lmax   
    
    #Lmax = 8.0
    #kold = 0.5    
    
    xfit,yfit = langguess_visual(x,y,kold)
    
    
    #duncan style minimization
    '''oldk = 0
    newk = 0.05
    stepsize = 0.001


    while lang_error(x,y,e,oldk) > lang_error(x,y,e,newk):
        oldk = newk
        newk = oldk + stepsize
    
    print "k for %s ="%(ID),newk   
    
    xfit,yfit = langguess_visual(x,y,oldk)'''
    
    return xfit,yfit
        

    
def langguess_visual(x,y,k):
    global Lmax
    
    xfit = np.arange(0, max(x), 0.01)
    yfit = []
    for i in xfit:
        yfit.append(y[0] + Lmax*(k*i)/(1+k*i))
        
    return xfit,yfit
    
def change_in_error(k,x,y,e):
    #calculates the change in error caused by change in k for the langmuir isotherm (dE/dK)    
    global Lmax    
    dE = 0.0
    Ltop = 0.0
    Lbot = 0.0    
    
    
    for i in range(0, len(x)):
        Ltop += (1/e[i]**2)*((k*x[i])/(k*x[i]+1))*(y[i]-y[0])
        Lbot += (1/e[i]**2)*((k*x[i])/(k*x[i]+1))**2
    
    Lmax = float(Ltop)/float(Lbot)    
    
    for i in range(0,len(x)):
        dE += (2/(e[i]**2))*((Lmax*k*x[i])/(1+k*x[i]) + y[0] -y[i])*((Lmax*x[i])/(1+k*x[i]))*(1-k/(1+k*x[i]**2))        
        #dE += (1/e[i])*((k*x[i])/(1+k*x[i])+y[0]-y[i])*(x[i])/(1+k*x[i])-(k*(x[i]**2))/(1+k*x[i])**2
    #print dE
    
    return dE
     
def figmake():
    ##Makes the figure and handles all of the data.
    ##Array names for specie nnnn use the following nomenclature: 
    ##nnnnX = conc
    ##nnnnY = force
    ##nnnnE = error in force
    ##nnnnExt = phi
    ##nnnn_EPM = energy per monomer
    ##nnnn_EErr = error in energy

    plt.figure(figsize=(12.0,7.50)) 
    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}

    plt.rc('font', **font)
    

    ##makes the tmaofig
    TMAOX = [0,1,2,3,4,5,6]
    TMAOY = [78.5, 80.81, 82.53, 83.29, 82.97, 84.91, 83.97]
    TMAOE = [0.2,2.77,1.272,3.42,2.22,2.80,3.44]
    TMAOExt = []
    TMAO_EPM = []
    TMAO_ERR = []

    for i in TMAOY:
        TMAOExt.append(fractional_extend_calc(i))

    for i in range(0, len(TMAOY)):
        TMAO_EPM.append(energy_per_monomer(TMAOY[i],TMAOExt[i]))
    
    for i in range(0,len(TMAOY)):
        TMAO_ERR.append(TMAO_EPM[i]/TMAOY[i]*TMAOE[i])
        

    Txfit, Tyfit = linreg_waterzero(TMAOX,TMAO_EPM,'TMAO')
    Tx2, Ty2 = linreg_statfit(TMAOX, TMAO_EPM, TMAO_ERR, 'TMAO')
    TxLang, TyLang = langguess_statfit(TMAOX, TMAO_EPM, TMAO_ERR, "TMAO")
        
    plt.errorbar(TMAOX,TMAO_EPM,yerr=TMAO_ERR,fmt='o',color = '#00FF00',markersize = 10)
    #plt.plot(Txfit, Tyfit, color = '#00FF00',label = 'TMAO')
    #plt.plot(Tx2, Ty2, color = "#00AA33", label = 'TMAO2')
    plt.plot(TxLang, TyLang, color = "#00FF00", label = "TMAO")
    
    
    
    


    #Ureafig
    UreaX = [0,1, 2,4,6,10]
    UreaY = [78.5, 75.27, 74.36 , 74.720,72.271,73.96]
    UreaE = [0.2, 3.124, 2.312, 2.647, 2.802,2.92]
    UreaExt = []
    Urea_EPM = []
    Urea_ERR = []

    for i in TMAOY:
        UreaExt.append(fractional_extend_calc(i))

    for i in range(0, len(UreaY)):
        Urea_EPM.append(energy_per_monomer(UreaY[i],UreaExt[i]))
    
    for i in range(0,len(UreaY)):
        Urea_ERR.append(Urea_EPM[i]/UreaY[i]*UreaE[i])

    Uxfit,Uyfit = linreg_waterzero(UreaX,Urea_EPM,'Urea')
    Ux2, Uy2 = linreg_statfit(UreaX, Urea_EPM, Urea_ERR, 'Urea')
    UxLang,UyLang = langguess_statfit(UreaX, Urea_EPM, Urea_ERR, "UREA")

    plt.errorbar(UreaX,Urea_EPM,yerr=Urea_ERR,fmt='^',color = '#FF0000',markersize = 10)
    #plt.plot(Uxfit, Uyfit, color = '#FF0000', linestyle = ':', label = "Urea")
    plt.plot(UxLang,UyLang, color = '#FF0000', label = "Urea")
    #plt.plot(Ux2, Uy2, color = "#AA0033", label = "Urea")

    ##makes the NaClfig
    NACLX = [0,1,2,5]
    NACLY = [78.5,81.08223395,83.6644707,89]
    NACLExt = []
    NACL_EPM = []  

    for i in NACLY:
        NACLExt.append(fractional_extend_calc(i))

    for i in range(0, len(NACLY)):
        NACL_EPM.append(energy_per_monomer(NACLY[i],NACLExt[i]))
    
    Nxfit, Nyfit = linreg_waterzero(NACLX,NACL_EPM, 'NaCl')
    
    plt.plot(NACLX, NACL_EPM,'rs', color = '#0000FF',markersize = 10)    
    plt.plot(Nxfit, Nyfit, color = '#0000FF', linestyle = '--',label = "NaCl")    
    
    
    #final commands
    
    plt.legend(bbox_to_anchor=(0, 0, 1, 1), loc = 1)
    plt.xlabel('Concentration of Solute (mol/L)')
    plt.ylabel('Extension Energy (KT)')
    plt.axis([0,11,3.4,4.6])
    
    plt.show()
    #plt.savefig('all_energies')
    
    print 'TERR', TMAO_ERR
    print 'TEng:%s\n' %(TMAO_EPM)
    print 'UERR', Urea_ERR
    print 'Ueng:%s\n' %(Urea_EPM)
    print 'NEng:%s\n' %(NACL_EPM)
       
      

figmake()
