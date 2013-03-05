#/usr/bin/env python3
import sys
from math import sqrt
dim = 100
if len(sys.argv) > 1:
    dim = int(sys.argv[1])
def rangle():
    for m in range(dim):
        for n in range(dim):
            print("v",float(m)/dim,float(n)/dim,0)
    
    def ind(m,n): return m + n*dim+1
    for m in range(dim-1):
        for n in range(dim-1):
            print("f",ind(m,n),ind(m,n+1),ind(m+1,n))
            print("f",ind(m,n+1),ind(m+1,n+1),ind(m+1,n))

def equal():
    for m in range(dim):
        parity = (m%2 == dim%2)
        if parity:
            print("v",float(m)/dim,0,0)
            
        for n in range(dim):
            if parity:
                print("v",float(m)/dim,(float(n)+.5)/dim/(sqrt(3)/2),0)
            else:
                print("v",float(m)/dim,float(n)/dim/(sqrt(3)/2.0),0)
        if not parity:
            print("v",float(m)/dim,(1.0-.5/dim) *2/sqrt(3),0)

    def ind(m,n): 
        prevlevels = int(m/2)
        buf = 1+prevlevels * (2*dim+2)
        if dim%2==0:
            if m%2==1:
                return buf +dim+n+1
            else:
                return buf + n
        else:
            if m%2==1:
                return buf +dim+n+1
            else:
                return buf + n
    print("Faces")
    for m in range(dim-1):
        parity = (m%2 == dim%2)
        
        if dim %2!=0:
            if m%2!=0:
                print("f",ind(m,1),ind(m,0),ind(m+1,0))
            else:
                print("f",ind(m,0),ind(m+1,0),ind(m+1,1))
        else:
            if m%2==0:
                print("f",ind(m,1),ind(m,0),ind(m+1,0))
            else:
                print("f",ind(m,0),ind(m+1,0),ind(m+1,1))
        for n in range(dim-1):
            if not parity:
                print("f",ind(m,n+1),ind(m,n),ind(m+1,n+1))
                print("f",ind(m,n+1),ind(m+1,n+1),ind(m+1,n+2))
            else:
                print("f",ind(m,n+1),ind(m+1,n),ind(m+1,n+1))
                print("f",ind(m,n+2),ind(m,n+1),ind(m+1,n+1))

        if dim %2!=0:
            if m%2!=0:
                print("f",ind(m,dim),ind(m+1,dim-1),ind(m+1,dim))
            else:
                print("f",ind(m,dim),ind(m,dim-1),ind(m+1,dim))
        else:
            if m%2==0:
                print("f",ind(m,dim),ind(m+1,dim-1),ind(m+1,dim))
            else:
                print("f",ind(m,dim),ind(m,dim-1),ind(m+1,dim))




#rangle()
equal();
