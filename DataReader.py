import os

class DataReader():
    
    def __init__(self, *args, **kwargs):
        print "DataReader initiated"
    
    def xye_in(self, fin):
        dataFile = open(fin)
        
        x = []
        y = []
        e = []
        
        for i in dataFile:
            try:
                line = i.split()
                x.append(float(line[0]))
                y.append(float(line[1]))
                e.append(float(line[2]))
                
            except(IndexError):
                print "not all columns exist in file %s" %(fin)

        
        return x, y, e
    
    
    
    
    
if __name__ == "__main__":
    tmp = DataReader(None)