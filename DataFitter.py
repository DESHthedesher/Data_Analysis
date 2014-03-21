
class DataFitter():

    def __init__(self, *args, **kwargs):
        global Lmax
        Lmax = 0
        print "DataFitter initiated"
        
        
      
    def linefit_simple(self,x,y):
        ### Returns Chi-Square fit parameters to errorless data  
        slope = 0
        intercept = 0
        chi2 = 0
        
        return slope, intercept, chi2




    
    def linefit_errorbars(self,x,y,e):
        # Returns Chi-Square fir parameters to data including error
        #uses equations 6-12 from bevington for slope paramters
        #equations 6-25 for error parameters
    
        ##compute the parameters
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
        intercept = (1/delta)*(sumxxs*sumys - sumxs*sumxys)
        slope = (1/delta)*(sums*sumxys - sumxs*sumys)
        
        
        
        #compute error and chi2
        intercepterror = ((1/delta)*(sumxxs))**0.5
        slopeerror = ((sums)/float(delta))**0.5
        
        chi2 = 0.0
        
        for i in range(0,len(x)):
            fi = intercept + slope*x[i]
            chi2 += (((y[i]-fi)**2)/e[i]**2)
        
        
        return slope, slopeerror, intercept, intercepterror, chi2

    
    def langmuir_errorfit(self, x,y,e):
        #compute the least-squares fit to Langmuir Isotherm by analytical and gradient based methods
        global Lmax
        
        kold = 0
        knew = 0.01
        step_size = 0.001
        precision = 0.00000001    
        
        avgY = 0.0
        for i in y:
            avgY += i - y[0]
    
        while abs(knew - kold) > precision:
            kold = knew
            de = self.change_in_error(kold,x,y,e)
            knew = kold + step_size*de
            
        chi2 = 0.0    
        for i in range(0,len(x)):
            fi = y[0] + Lmax*(kold*x[i])/(1+kold*x[i])
            chi2 += (((y[i]-fi)**2)/e[i]**2)
            
        return kold, Lmax, chi2 
    
        
    def change_in_error(self, k,x,y,e):
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
    
        return dE


    def pearson_r2(self,x,y,slope,intercept):
        #computes the pearson coefficient for data
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
       


if __name__ == "__main__":
    DataFitter(None)
        

    
