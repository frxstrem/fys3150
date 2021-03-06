#include <iostream>
#include <armadillo>
#include <cmath>
#include <time.h>
#include <fstream>

#include "../code-fredrik/comp_eig.hh"

using namespace std;
using namespace arma;

void RHO_A_FILL(vec &rho, mat &A, int N,double rhoN); //rho, kind of like a linespace
                                                 //A, Tridiagonal matrix
void Maxoff(mat &A, int N,int &k, int &l, double &max);        //Finds the max element of a matrix

int main(){
    int values[] = { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100 };
    fstream Outfile;
    Outfile.open("steps.dat",ios::out);
    Outfile << "N     epsilon    steps  step_time  error"<<endl;
    for(int iter = 0; iter < 20; iter++){
        int N = values[iter]; //matrix size; N x N
        double rhoN = 6;
        double eps = 1E-10;
        double max;
        int k,l;

        double total_error = 0;
        double total_time = 0;
        int total_iterations = 0;

        // run 10 times and average
        int zcount = 10;
        for(int z = 0; z < zcount; z++) {
          std::cout << N << std::endl;
          mat A(N,N,fill::zeros); //indexes go from (0) to (N-1)
          mat S(N,N,fill::eye);
          vec rho(N,fill::zeros);
          vec eigen(N);


          RHO_A_FILL(rho,A,N,rhoN);
          Maxoff(A,N,k,l,max);
          //A.print();
          int iterations = 0;

          // solve with Armadillo's eig_sym
          vec arma_eigenvalues;
          mat arma_eigenvectors;
          eig_sym(arma_eigenvalues, arma_eigenvectors, A);

          double tau, t, s, c, il, ik, kk, ll, s_ik, s_il;
          double start, finish;

          start = clock();                      //clock value before eigen solve
          while(max > eps){
              tau = (A(l,l)-A(k,k))/(2.*A(k,l));
              if(tau>0){
                  t = 1.0/(tau + sqrt(1.0 + tau*tau));
              } else{
                  t = -1.0/( -tau + sqrt(1.0 + tau*tau));
              }

              //cosine and sine
              c = 1./sqrt(1.+t*t);
              s = t*c;

              //Jacobi rotating A round theta in N-dim space
              for(int i = 0; i<N; i++){
                  if ((i != k) && (i !=l)){

                      ik = A(i,k)*c - A(i,l)*s;
                      il = A(i,l)*c + A(i,k)*s;
                      A(i,k) = ik;
                      A(i,l) = il;
                      A(k,i) = ik;
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
              A(l,l) = ll ;
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

          // eigen.print();

          // compare eigenvectors from own method with armadillo
          total_error += comp_eig(eigen, S, arma_eigenvalues, arma_eigenvectors);

          total_time += (finish -start)/CLOCKS_PER_SEC/iterations;
          total_iterations += iterations;

          /*
          fstream outfile;
          outfile.open("eigenvectors.dat",ios::out);

          for (int i = 0; i<N; i++){
              for (int j = 0; j<N; j++){
                  if(j%N == 0){outfile<<endl;}
                  outfile << S(i,j)<<" ";
                  }
              }
          outfile.close();
          */
        }


        Outfile <<N <<" "<<eps<<" "<<((double)total_iterations / zcount)<<" "<<(total_time / zcount) <<" "<<(total_error / zcount)<<endl;

    //    Outfile.close();
    }
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
