#include <iostream>
#include <armadillo>
#include <cmath>
#include <time.h>
#include <fstream>
using namespace std;
using namespace arma;

void RHO_A_FILL(vec &rho, mat &A, int N,double rhoN); //rho, kind of like a linespace
                                                 //A, Tridiagonal matrix
void Maxoff(mat &A, int N,int &k, int &l, double &max);        //Finds the max element of a matrix

int main(){
    int N = 100; //matrix size; N x N
    double rhoN = 6;
    double max;
    int k,l;

    mat A = mat(N,N,fill::zeros); //indexes go from (0) to (N-1)
    mat S = mat(N,N,fill::eye);
    vec rho = vec(N,fill::zeros);
    fvec eigen = fvec(N);


    RHO_A_FILL(rho,A,N,rhoN);
    Maxoff(A,N,k,l,max);
    //A.print();
    int iterations = 0;

    double eps = 1E-10;

    double tau, t, s, c, il, ik, kk, ll, s_ik, s_il;
    double start, finish;

    start = clock();                      //clock value before eigen solve
    while(max > eps){
        tau = (A(l,l)-A(k,k))/(2.*A(k,l));
        if(tau>0){
             t = 1.0/(tau + sqrt(1.0 + tau*tau));}
        else{t = -1.0/( -tau + sqrt(1.0 + tau*tau));}
     
        //cosine and sine
        c = 1./sqrt(1.+t*t);
        s = t*c;
     
        //Jacobi rotating A round theta in N-dim space
        for(int i = 0; i<N; i++){
            if ((i != k) && (i !=l)){
                ik = A(i,k)*c - A(i,l)*s;
                il = A(i,l)*c + A(i,k)*s;
                A(i,k) = ik;
                A(k,i) = ik;
                A(i,l) = il;
                A(l,i) = il;
                }
            s_ik = S(i,k);
            s_il = S(i,l);
            S(i,k) = c*s_ik - s*s_il;
            S(i,l) = c*s_il + s*s_ik; 
            }

        kk = A(k,k)*c*c - 2.*A(k,l)*c*s + A(l,l)*s*s;
        ll = A(l,l)*c*c + 2.*A(k,l)*c*s + A(k,k)*s*s;
     
        A(k,k) = kk;
        A(l,l) = ll;
        A(k,l) = 0;
        A(l,k) = 0;

        iterations++;
        Maxoff(A,N,k,l,max);


        } //end of while
    finish = clock();                    //clock value after eigen solve

    //[eigen] is now eigenvalues
    for(int i = 0; i<N;i++){
        eigen(i) = A(i,i);
        }

    eigen = sort(eigen);
    eigen.print();

    double time = (finish -start)/CLOCKS_PER_SEC/iterations;
    
    fstream outfile;
    outfile.open("eigenvectors.dat",ios::out);

    for (int i = 0; i<N; i++){
        for (int j = 0; j<N; j++){
            if(j%N == 0){outfile<<endl;}                        
            outfile << S(i,j)<<" ";
            }
        }
    outfile.close();



    //fstream Outfile;
    //Outfile.open("steps.dat",ios::out);
    //Outfile << "N     epsilon    steps  step_time"<<endl;

    //Outfile <<N <<" "<<eps<<" "<<iterations<<" "<<time;
    //Outfile<<endl;
        
    //Outfile.close();

} //end of main


void RHO_A_FILL(vec &rho, mat &A, int N,double rhoN){
    double h = rhoN/(N);
    for(int i = 0; i < N;i++){
        rho(i) = i*h;
        A(i,i) = 2./(h*h)+(rho(i)*rho(i));
        if(i<N-1){
            A(i+1,i) = -1./(h*h);
            A(i,i+1) = -1./(h*h);
            }
        }
    }

void Maxoff(mat &A, int N,int &k, int &l, double &max){
    max = 0;
    for(int i = 0; i < N ; i++){
        for(int j = i+1; j < N ; j++){
            if ( A(i,j)*A(i,j) > max){
                max =A(i,j)*A(i,j);
                k = i; 
                l = j;
                }
            }
        }
    
    }
