import numpy as np
import sys
import math


for i in range(len(sys.argv)-1):

    a = float(sys.argv[i+1]) *1.08/0.8
    b = math.sqrt(0.8)*0.58/ math.sqrt(float(sys.argv[i+1]))

    n = 1.0/math.sqrt(4-8*math.exp(-2*a*b*b) + 4*math.exp(-4*a*b*b))

    b2 = b * math.sqrt(2)

    print(4)
    print(a,"\t",n,"\t",b,"\t",b,"\t",0)
    print(a,"\t",-n,"\t",-b,"\t",b,"\t",0)
    print(a,"\t",-n,"\t",b,"\t",-b,"\t",0)
    print(a,"\t",n,"\t",-b,"\t",-b,"\t",0)
    print(4)
    print(a,"\t",n,"\t",0,"\t",b,"\t",b)
    print(a,"\t",-n,"\t",0,"\t",b,"\t",-b)
    print(a,"\t",-n,"\t",0,"\t",-b,"\t",b)
    print(a,"\t",n,"\t",0,"\t",-b,"\t",-b)
    print(4)
    print(a,"\t",n,"\t",b,"\t",0,"\t",b)
    print(a,"\t",-n,"\t",-b,"\t",0,"\t",b)
    print(a,"\t",-n,"\t",b,"\t",0,"\t",-b)
    print(a,"\t",n,"\t",-b,"\t",0,"\t",-b)
    print(4)
    print(a,"\t",n,"\t",b2,"\t",0,"\t",0)
    print(a,"\t",n,"\t",-b2,"\t",0,"\t",0)
    print(a,"\t",-n,"\t",0,"\t",b2,"\t",0)
    print(a,"\t",-n,"\t",0,"\t",-b2,"\t",0)
    print(4)
    print(a,"\t",n,"\t",b2,"\t",0,"\t",0)
    print(a,"\t",n,"\t",-b2,"\t",0,"\t",0)
    print(a,"\t",-n,"\t",0,"\t",0,"\t",b2)
    print(a,"\t",-n,"\t",0,"\t",0,"\t",-b2)