#include "matrix.h"
#include <iostream>
#include <vector>
#include <random>
#include <math.h>
#include <time.h>


QSMatrix<double> matrix_rand_ortho(int dim){ // Generate a random orthonormal matrix with Householder transformations

    QSMatrix<double> A(dim, dim, 'i'); // The random orthnormal matrix
    QSMatrix<double> d(1,dim,0);

    srand(time(NULL));
    double number = rand()/(double)RAND_MAX;
    if (number < 0.5){
        d(0,dim-1)=-1;
    }
    else{
        d(0, dim-1)=1;
    }

    double norm_x=0;
    double sign=0;
    double beta = 0;


    for (int k=dim-2; k>-1; k--){
        // Generate random Householder transformation
        QSMatrix<double> x(1, dim-k, 'r');

        // norm of x
        norm_x=0;
        for (int i=0; i<dim-k; i++){
            norm_x += x(0, i)*x(0, i);
        }
        norm_x = sqrt(norm_x);

        sign = 1;
        if (x(0,0)<0){
            sign = -1;
        }

        norm_x *= sign;
        d(0,k) = -sign;

        x(0,0) += norm_x;

        beta = norm_x * x(0,0);

        // apply the transformation
        QSMatrix<double> A_(dim-k, dim, 0);
        for (int i=0; i<dim-k; i++){
            for (int j=0; j<dim; j++){
                A_(i, j)=A(i+k, j);
            }
        }
        QSMatrix<double> y(1,dim,0);
        for (int i=0; i<dim-k; i++){
            for (int j=0; j<dim; j++){
                y(0,j)+=x(0,i)*A_(i,j);
            }
        }

        y = y/beta;

        QSMatrix<double> outer_xy(dim-k, dim, 0);
        for (int i=0; i<dim-k; i++){
            for (int j=0; j<dim; j++){
                outer_xy(i,j)=x(0,i)*y(0,j);
            }
        }
        A_-=outer_xy;

        for (int i=0; i<dim-k; i++){
            for (int j=0; j<dim; j++){
                A(i+k, j)=A_(i, j);
            }
        }
    }

    // change sign of rows
    for (int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            A(i,j) = A(i,j)*d(0,i);
        }
    }

    return A;
}


int main(int argc, char **argv) {


    // Create a var P
    QSMatrix<double> P(1,1,1);
    // Generate a random orthonormal matrix
    P = matrix_rand_ortho(10);

    // Print the matrix
    P.print_mat();

    // Pour sortir l'element à la ligne 2 et colonne 3
    std::cout<<P(1,2)<<std::endl;


    return 0;
}
