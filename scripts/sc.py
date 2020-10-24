import numpy as np
import sys
import math


for i in range(len(sys.argv)-1):

    a = float(sys.argv[i+1]) / 2.0832417
    b = 0.7071/ math.sqrt(float(sys.argv[i+1]))

    n = 1.0/math.sqrt(2-2*math.exp(-2*a*b*b))

    print(2)
    print(a,"\t",n,"\t",b,"\t",0,"\t",0)
    print(a,"\t",-n,"\t",-b,"\t",0,"\t",0)
    print(2)
    print(a,"\t",n,"\t",0,"\t",b,"\t",0)
    print(a,"\t",-n,"\t",0,"\t",-b,"\t",0)
    print(2)
    print(a,"\t",n,"\t",0,"\t",0,"\t",b)
    print(a,"\t",-n,"\t",0,"\t",0,"\t",-b)



