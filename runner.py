import Constants
import DataReader
import DataFitter

DR = DataReader.DataReader(None)
F = DataFitter.DataFitter(None)

def main():
    
    runfor("tmao.dat")
    runfor("urea.dat")
    #runfor("nacl.dat")
    

def runfor(ID):
    print "\n"
    print "***************************************************"
    print "Fits for %s" %(ID)
    
    ##get data
    x,y,e = DR.xye_in(ID)

    ##get linefit
    slope,sloperror,intercept,intercepterror,chi2 = F.linefit_errorbars(x,y,e)
    
    print "Linear fit data for %s:"%(ID)
    print "best slope is %s +/- %s" %(slope, sloperror)
    print "best intercept is %s +/- %s" %(intercept,intercepterror)
    print "chi2 is %s" %(chi2)
    print "\n"
   
   
    ##get langfit
    k, Lmax ,chi2 = F.langmuir_errorfit(x,y,e)
    
    print "Saturation fit data for %s:"%(ID)
    print "k =",k
    print "Lmax =",Lmax   
    print "chi2 is %s" %(chi2) 
    print "\n"
    
main()
