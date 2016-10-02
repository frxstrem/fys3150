from matplotlib.pyplot import *
from numpy import *
S1 = loadtxt("eigenvectorsminmin.dat")
E1 = loadtxt("eigenvaluesminmin.dat")
S2 = loadtxt("eigenvectorsmin.dat")
E2 = loadtxt("eigenvaluesmin.dat")
S3 = loadtxt("eigenvectors1.dat")
E3 = loadtxt("eigenvalues1.dat")
S4 = loadtxt("eigenvectors5.dat")
E4 = loadtxt("eigenvalues5.dat")
index1 = argmin(E1)
index2 = argmin(E2)
index3 = argmin(E3)
index4 = argmin(E4)



vector1 = S1[:,index1]/sqrt(60./800.)
vector2 = S2[:,index2]/sqrt(60./800.)
vector3 = S3[:,index3]/sqrt(60./800.)
vector4 = S4[:,index4]/sqrt(60./800.)
x =linspace(0,60,len(vector1))

##############Writing to file for Latex representation
with open("vectors.dat","w") as fp:
    fp.write("rho  rhoN  omega  vec \n")
    for i in range(len(vector1)):
        fp.write("%f %d %f %f \n" %(x[i], 60, 0.01, vector1[i]))
    for i in range(len(vector2)):
        fp.write("%f %d %f %f \n" %(x[i], 60, 0.5, vector2[i]))
    
    for i in range(len(vector3)):
        fp.write("%f %d %f %f \n" %(x[i], 60, 1, vector3[i]))
    

    for i in range(len(vector4)):
        fp.write("%f %d %f %f \n" %(x[i], 60, 5, vector4[i]))
#################




#################PLOTTING
plot(x,abs(vector1/sqrt(60./800.)),label="$\\omega_r$ = 0.01")

plot(x,abs(vector2/sqrt(60./800.)),label="$\\omega_r$ = 0.5")

plot(x,abs(vector3/sqrt(60./800.)),label="$\\omega_r$ = 1")

plot(x,abs(vector4/sqrt(60./800.)),label="$\\omega_r$ = 5")
legend()
xlabel("$\\rho$")
ylabel("$\\Psi$")
title("Eigenvectors of ground state as a function of \\omega")
show()
