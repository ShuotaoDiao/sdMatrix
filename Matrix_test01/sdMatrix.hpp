//
//  sdMatrix.hpp
//  Matrix_test01
//
//  Created by Shuotao Diao on 2/8/17.
//  Copyright © 2017 Shuotao Diao. All rights reserved.
//  Matrix Class

#ifndef sdMatrix_hpp
#define sdMatrix_hpp

#include <stdio.h>
#include <string>
#include <limits>
class sdMatrix{
private:
    std::string name_;
    double ** entries_;//2 dimensional array
    int row_;
    int col_;
    int i_ = -1;
    int j_ = -1;
    int i_write = -1;
    int j_write = -1;
    //I need a stack
    int ii[5] = {-1,-1,-1,-1,-1};
    int jj[5] = {-1,-1,-1,-1,-1};
    int index = 0;
public:
    sdMatrix();//constructor
    sdMatrix(double *entries, int row, int col);//constructor 1
    sdMatrix(int entry, int row, int col);//constuctor 2
    sdMatrix(std::string name, double *entries, int row, int col);//constructor 3
    sdMatrix(const sdMatrix & matrix );//copy constructor 
    ~sdMatrix();//destructor
    sdMatrix operator+ (const sdMatrix& matrix);//matrix plus
    sdMatrix operator+ (double scalar);//matrix plus a scalar
    friend sdMatrix operator+ (double scalar, const sdMatrix& matrix);//a scalar plus a matrix
    sdMatrix operator- (const sdMatrix& matrix);//matrix minus
    sdMatrix operator- (double scalar);//matrix minus a scalar
    friend sdMatrix operator- (double scalar, const sdMatrix& matrix); //a scalar minus a matrix
    sdMatrix operator* (const sdMatrix& matrix);//matrix multiplication
    sdMatrix operator* (double scalar);//matrix times a scalar
    friend sdMatrix operator* (double scalar, const sdMatrix& matrix);//a scalar times a matrix
    sdMatrix operator/ (const sdMatrix& matrix);//matrix division
    sdMatrix operator/ (double scalar);//matrix divided by a constant 
    friend sdMatrix times(const sdMatrix& matrix1, const sdMatrix& matrix2);//element by element muliplication
    sdMatrix& operator= (const sdMatrix& matrix);//matrix assginment
    sdMatrix& operator= (const double element);//assign a value to one element
    sdMatrix& operator() (int row, int col);//access single element (update the set of indices)
    sdMatrix& rowSwitch(int row1, int row2);//row switching 
    double det();//compute the determinant of sqaure matrix 
    sdMatrix inv();//invert the square matrix
    sdMatrix T();
    int length();//give the length of the vector
    double** toDouble();//output the matrix in double arrays
    double element();//output the value in the desired location
    void getValue() const;//print out the matrix
    void setName(std::string name);//set the name of the matrix
};
#endif /* sdMatrix_hpp */
