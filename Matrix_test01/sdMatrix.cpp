//
//  sdMatrix.cpp
//  Matrix_test01
//
//  Created by Shuotao Diao on 2/8/17.
//  Copyright © 2017 Shuotao Diao. All rights reserved.
//
#include <iostream>
#include <string>
#include <limits>
#include "sdMatrix.hpp"

//default constructor
sdMatrix::sdMatrix(){
    name_ = "Matrix";
    entries_ = new double*[1];
    entries_[0] = new double[1];
    entries_[0][0] = 0;
    row_ = 1;
    col_ = 1;
}

//constructor 1
sdMatrix::sdMatrix(double* entries, int row, int col){
    name_ = "Matrix";
    row_ = row;//number of rows
    col_ = col;//number of columns
    //right now, I don't know how to check the whether #entries == row * col
    entries_ = new double*[row_];
    for(int i = 0; i < row_; ++i){
        entries_[i] = new double[col_];
        for(int j = 0; j < col_; ++j){//row by row 
            entries_[i][j] = entries[j + col_ * i];//assign the vlaue to each entry
        }
    }
}
//constructor 2
sdMatrix::sdMatrix(int entry, int row, int col){
    name_ = "Matrix";
    row_ = row;//number of rows
    col_ = col;//number of columns
    entries_ = new double*[row_];
    for(int i = 0; i < row_; ++i){
        entries_[i] = new double[col_];
        for(int j = 0; j < col_; ++j){
            entries_[i][j] = entry;//assignment the value to each entry
        }
    }
}
//constructor 3
sdMatrix::sdMatrix(std::string name, double *entries, int row, int col){
    name_ = name;
    row_ = row;//number of rows
    col_ = col;//number of columns
    entries_ = new double*[row_];
    for(int i = 0; i < row_; ++i){
        entries_[i] = new double[col_];
        for(int j = 0; j < col_; ++j){
            entries_[i][j] = entries[j + col_ * i];//assignment the value to each entry
        }
    }
}
//copy constructor
sdMatrix::sdMatrix(const sdMatrix & matrix){
    name_ = "Matrix";
    row_ = matrix.row_;
    col_ = matrix.col_;
    entries_ = new double*[row_];
    for(int i = 0; i < row_; ++i){
        entries_[i] = new double[col_];
        for(int j = 0; j < col_; ++j){
            entries_[i][j] = matrix.entries_[i][j];//assign the vlaue to each entry
        }
    }
}

//destrutor
sdMatrix::~sdMatrix(){
    for (int i = 0; i < row_; ++i){
        delete [] entries_[i];//clean up the memory of each row
    }
    delete [] entries_;//release the memory
}

//matrix plus
sdMatrix sdMatrix::operator+ (const sdMatrix & matrix){
    sdMatrix temp(*this);//copy the object to the temp
    if(row_ == matrix.row_ && col_ == matrix.col_){
        for(int i = 0; i < row_; ++i){
            for(int j = 0; j < col_; ++j){
                temp.entries_[i][j] += matrix.entries_[i][j];//add the values element by element
            }
        }
        return temp;//return temp
    }
    else{
        std::cout << "Warning(plus):The sizes of two matrices are not equal!" << std::endl;//print out the warning
        sdMatrix matrixNAN(0,1,1);
        matrixNAN.entries_[0][0] = std::numeric_limits<double>::quiet_NaN();
        return matrixNAN;
    }
    
}
//matrix plus a scalar
sdMatrix sdMatrix::operator+ (double scalar){
    sdMatrix temp(*this);//copy the object to the temp
    for(int i = 0; i < row_; ++i){
        for(int j = 0; j < col_; ++j){
            temp.entries_[i][j] += scalar;//add the values element by element
        }
    }
    return temp;//return temp
}
//a scalar plus a matrix
sdMatrix operator+ (double scalar, const sdMatrix & matrix){
    sdMatrix temp(matrix);
    for(int i = 0; i < temp.row_; ++i){
        for(int j = 0; j < temp.col_; ++j){
            temp.entries_[i][j] += scalar;//add the scalar element by element
        }
    }
    return temp;//return temp
}
//matrix minus a matrix
sdMatrix sdMatrix::operator- (const sdMatrix & matrix){
    sdMatrix temp(*this);//copy the object to the temp
    if(row_ == matrix.row_ && col_ == matrix.col_){
        for(int i = 0; i < row_; ++i){
            for(int j = 0; j < col_; ++j){
                temp.entries_[i][j] -= matrix.entries_[i][j];
            }
        }
        return temp;//return temp
    }
    else{
        std::cout << "Warning(minus):The sizes of two matrix are not equal!" << std::endl;//print out the warning
        sdMatrix matrixNAN(0,1,1);
        matrixNAN.entries_[0][0] = std::numeric_limits<double>::quiet_NaN();
        return matrixNAN;
    }
}
//matrix minus a scalar
sdMatrix sdMatrix::operator- (double scalar){
    sdMatrix temp(*this);//copy the object to the temp
    for(int i = 0; i < row_; ++i){
        for(int j = 0; j < col_; ++j){
            temp.entries_[i][j] -= scalar;//substract by the scalar element by element
        }
    }
    return temp;//return temp
}
//a scalar minus a matrix
sdMatrix operator- (double scalar, const sdMatrix & matrix){
    sdMatrix temp(scalar, matrix.row_, matrix.col_);//declare a matrix
    for(int i = 0; i < temp.row_; ++i){
        for(int j = 0; j < temp.col_; ++j){
            temp.entries_[i][j] -= matrix.entries_[i][j];//substract by the scalar element by element
        }
    }
    return temp;
}
//a scalar times a matrix
sdMatrix operator* (double scalar, const sdMatrix & matrix){
    sdMatrix temp(0, matrix.row_, matrix.col_);//declare a matrix
    for(int i = 0; i < matrix.row_; ++i){
        for(int j = 0; j < matrix.col_; ++j){
            temp.entries_[i][j] = matrix.entries_[i][j] * scalar;
        }
    }
    return temp;
}
//matrix multiplication
sdMatrix sdMatrix::operator* (const sdMatrix & matrix){
    sdMatrix temp(0, row_, matrix.col_);//copy the object to the temp
    if(col_ != matrix.row_){
        std::cout << "Warning(multiplication): The number of columns of the first matrix does not " << std::endl;
        std::cout << "eqaul to the number of rows of the second matrix" << std::endl;//print out the warning
        return temp;
    }
    for(int i = 0; i < row_; ++i){
        for(int j = 0; j < matrix.col_; ++j){
            for(int k = 0; k < col_; ++k){
                temp.entries_[i][j] += entries_[i][k] * matrix.entries_[k][j];
            }
        }
    }
    return temp;
}
//matrix times a scalar
sdMatrix sdMatrix::operator* (double scalar){
    sdMatrix temp(0,row_, col_ );
    for(int i = 0; i < row_; ++i){
        for(int j = 0; j < col_; ++j){
            temp.entries_[i][j] = entries_[i][j] * scalar;
        }
    }
    return temp;
}

//element by element multiplication
sdMatrix times(const sdMatrix & matrix1, const sdMatrix & matrix2){
    if(matrix1.row_ != matrix2.row_ || matrix1.row_ != matrix2.row_){
        sdMatrix temp(0, 1, 1);
        return temp;
    }
    sdMatrix temp(0,matrix1.row_, matrix1.col_);//declare a new matrix object
    for(int i = 0; i < temp.row_; ++i){
        for(int j = 0; j < temp.col_; ++j){
            temp.entries_[i][j] = matrix1.entries_[i][j] * matrix2.entries_[i][j];
        }
    }
    return temp;
}

//matrix division
sdMatrix sdMatrix::operator/ (const sdMatrix & matrix){
    if(row_ != col_){
        std::cout << "Division Warning: Input matrix is not a square matrix" << std::endl;
        sdMatrix A("null",0,1,1);
        return A;
    }
    if(matrix.row_ != matrix.col_){
        std::cout << "Division Warning: Input matrix is not a square matrix" << std::endl;
        sdMatrix A("null",0,1,1);
        return A;
    }
    if(row_ != matrix.row_){
        std::cout << "Division Warning: The sizes of two matrices don't match" << std::endl;
        sdMatrix A("null",0,1,1);
        return A;
    }
    //check if it is sqaure matrix
    if(row_ != col_){
        std::cout << "ERROR: Input should be sqaure matrix!" << std::endl;
        sdMatrix matrixNAN(0,1,1);
        matrixNAN.entries_[0][0] = std::numeric_limits<double>::quiet_NaN();
        return matrixNAN;
    }
    sdMatrix inverse(0, matrix.row_, matrix.col_);
    sdMatrix A = matrix;
    double answer = 1;
    //initialize the invserse matrix (identity matrix)
    for(int i = 0; i < matrix.row_; ++i){
        inverse.entries_[i][i] = 1;
    }
    //get the upper triangle matrix
    for (int it = 0; it < matrix.row_ - 1; ++it){
        sdMatrix L(0,matrix.row_,matrix.row_);//build a lower triangle matrix
        //check if the triangluar entry is 0 if yes, do the row switching
        if(A.entries_[it][it] == 0){
            bool test = true;
            for(int k = it + 1; k < matrix.row_ ; ++k){
                if(A.entries_[k][it] != 0){
                    A = A.rowSwitch(it,k);//
                    inverse = inverse.rowSwitch(it, k);
                    answer = -1 * answer;
                    test = false;
                }
            }
            if(test){//no entry at the col "it" that is non zero, then the det is 0
                answer = 0;
                sdMatrix matrixNAN(0,1,1);
                matrixNAN.entries_[0][0] = std::numeric_limits<double>::quiet_NaN();
                return matrixNAN;
            }
        }
        for( int i = 0; i < matrix.row_; ++i){
            for(int j = 0; j < matrix.col_; ++j){
                if(i == j) L.entries_[i][j] = 1;
                if(i > j && j == it) {
                    {
                        L.entries_[i][j] = -A.entries_[i][j] / A.entries_[it][it];
                    }
                }
                //std::cout << L.entries_[i][j] << " ";
            }
            //std::cout << std::endl;
        }
        A = L * A;//update the matrix
        inverse = L * inverse;//update the inverse matrix
        //A.getValue();
    }
    //A.getValue();
    //inverse.getValue();
    //transform the upper triangle matrix into indentity matrix
    for(int it = col_ - 1; it > 0; --it){
        sdMatrix L(0,row_,row_);//build a lower triangle matrix
        for( int i = 0; i < row_ ; ++i){
            for(int j = 0; j < col_ ; ++j){
                if(i == j) L.entries_[i][j] = 1;
                if(i < j && j == it) {
                    {
                        L.entries_[i][j] = -A.entries_[i][j] / A.entries_[it][it];
                    }
                }
                //std::cout << L.entries_[i][j] << " ";
            }
            //std::cout << std::endl;
        }
        A = L * A;//update the matrix
        inverse = L * inverse;
    }
    //A.getValue();
    //inverse.getValue();
    for(int i = 0; i < row_; ++i){
        for(int j = 0; j < col_; ++j){
            inverse.entries_[i][j] = inverse.entries_[i][j] / A.entries_[i][i];
        }
    }
    sdMatrix B = (*this) * inverse;
    return B;
}
    
sdMatrix sdMatrix::operator/ (double scalar){
    sdMatrix temp(0,row_, col_);
    for(int i = 0; i < row_; ++i){
        for(int j = 0; j < col_; ++j){
            temp.entries_[i][j] = entries_[i][j] / scalar;
        }
    }
    return temp;
}


//access the single element (update the set of indices)
sdMatrix& sdMatrix::operator() (int row, int col){
    if( row > row_ - 1 && col > col_ - 1) {
        std::cout << "Warning: Index Exceeds Matrix Dimension! " << std::endl;
        return *this;
    }
    i_ = row;
    j_ = col;
    ii[index] = row;
    jj[index] = col;
    index++;
    return *this;
}
//row switching
sdMatrix& sdMatrix::rowSwitch(int row1, int row2){
    if(row1 > row_ || row2 > row_) {
        std::cout << "Warning: Index Exceeds Matrix Dimension!" << std::endl;
        return *this;
    }
    double temp = 0;
    for(int i = 0; i < col_; ++i){
        temp = entries_[row1][i];
        entries_[row1][i] = entries_[row2][i];
        entries_[row2][i] = temp;
    }
    return *this;
}


//compute the determinant by doing row elimination
double sdMatrix::det(){
    //check if it is sqaure matrix
    if(row_ != col_){
        std::cout << "ERROR: Input should be sqaure matrix!" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    sdMatrix A = *this;//copy the matrix to Adouble
    double answer = 1;
    for (int it = 0; it < row_ - 1; ++it){
        sdMatrix L(0,row_,row_);//build a lower triangle matrix
        //check if the triangluar entry is 0 if yes, do the row switching
        if(A.entries_[it][it] == 0){
            bool test = true;
            for(int k = it + 1; k < row_ ; ++k){
                if(A.entries_[k][it] != 0){
                    A = A.rowSwitch(it,k);//
                    answer = -1 * answer;
                    test = false;
                }
            }
            if(test){//no entry at the col "it" that is non zero, then the det is 0
                answer = 0;
                return answer;
            }
        }
        for( int i = 0; i < row_; ++i){
            for(int j = 0; j < col_; ++j){
                if(i == j) L.entries_[i][j] = 1;
                if(i > j && j == it) {
                    {
                    L.entries_[i][j] = -A.entries_[i][j] / A.entries_[it][it];
                    }
                }
                //std::cout << L.entries_[i][j] << " ";
            }
            //std::cout << std::endl;
        }
        A = L * A;//update the matrix
        //A.getValue();
    }
    //compute the determinant by multiplying triangular entries of A
    for(int i = 0; i < row_; ++i){
        answer *= A.entries_[i][i];
    }
    return answer;
}

//invert the sqaure matrix
sdMatrix sdMatrix::inv(){
    //check if it is sqaure matrix
    if(row_ != col_){
        std::cout << "ERROR: Input should be sqaure matrix!" << std::endl;
        sdMatrix matrixNAN(0,1,1);
        matrixNAN.entries_[0][0] = std::numeric_limits<double>::quiet_NaN();
        return matrixNAN;
    }
    sdMatrix inverse(0, row_, col_);
    sdMatrix A = *this;
    double answer = 1;
    //initialize the invserse matrix (identity matrix)
    for(int i = 0; i < row_; ++i){
        inverse.entries_[i][i] = 1;
    }
    //get the upper triangle matrix
    for (int it = 0; it < row_ - 1; ++it){
        sdMatrix L(0,row_,row_);//build a lower triangle matrix
        //check if the triangluar entry is 0 if yes, do the row switching
        if(A.entries_[it][it] == 0){
            bool test = true;
            for(int k = it + 1; k < row_ ; ++k){
                if(A.entries_[k][it] != 0){
                    A = A.rowSwitch(it,k);//
                    inverse = inverse.rowSwitch(it, k);
                    answer = -1 * answer;
                    test = false;
                }
            }
            if(test){//no entry at the col "it" that is non zero, then the det is 0
                answer = 0;
                sdMatrix matrixNAN(0,1,1);
                matrixNAN.entries_[0][0] = std::numeric_limits<double>::quiet_NaN();
                return matrixNAN;
            }
        }
        for( int i = 0; i < row_; ++i){
            for(int j = 0; j < col_; ++j){
                if(i == j) L.entries_[i][j] = 1;
                if(i > j && j == it) {
                    {
                        L.entries_[i][j] = -A.entries_[i][j] / A.entries_[it][it];
                    }
                }
                //std::cout << L.entries_[i][j] << " ";
            }
            //std::cout << std::endl;
        }
        A = L * A;//update the matrix
        inverse = L * inverse;//update the inverse matrix
        //A.getValue();
    }
    //A.getValue();
    //inverse.getValue();
    //transform the upper triangle matrix into indentity matrix
    for(int it = col_ - 1; it > 0; --it){
        sdMatrix L(0,row_,row_);//build a lower triangle matrix
        for( int i = 0; i < row_ ; ++i){
            for(int j = 0; j < col_ ; ++j){
                if(i == j) L.entries_[i][j] = 1;
                if(i < j && j == it) {
                    {
                        L.entries_[i][j] = -A.entries_[i][j] / A.entries_[it][it];
                    }
                }
                //std::cout << L.entries_[i][j] << " ";
            }
            //std::cout << std::endl;
        }
        A = L * A;//update the matrix
        inverse = L * inverse;
    }
    //A.getValue();
    //inverse.getValue();
    for(int i = 0; i < row_; ++i){
        for(int j = 0; j < col_; ++j){
            inverse.entries_[i][j] = inverse.entries_[i][j] / A.entries_[i][i];
        }
    }
    return inverse;
}


//print out the matrix
void sdMatrix::getValue() const{
    std::cout << name_ << " = "<< std::endl;
    for(int i = 0; i < row_;++i){
        std::cout << "     ";
        for(int j = 0; j < col_; ++j){
            std::cout << entries_[i][j] << " ";//print out element by element
        }
        std::cout << std::endl;
    }
    
}
//assignemnt constructor
sdMatrix& sdMatrix::operator= (const sdMatrix & matrix){
    if(this == &matrix) return *this;//protect from self assignment
    row_ = matrix.row_;
    col_ = matrix.col_;
    entries_ = new double*[row_];
    for(int i = 0; i < row_; ++i){
        entries_[i] = new double[col_];
        for(int j = 0; j < col_; ++j){
            entries_[i][j] = matrix.entries_[i][j];//assign the vlaue to each entry
        }
    }
    return *this;
}
//assign the value to one element
sdMatrix& sdMatrix::operator= (const double element){
    //if (i_ == -1 && j_ == -1){
    if(index == 0){
        std::cout << "Incorrect Assignment!" << std::endl;
        return *this;
    }
    /*
    i_write = i_;
    j_write = j_;
     */
    i_write = ii[index - 1];
    j_write = jj[index - 1];
    this -> entries_[i_write][j_write] = element;//assign the value to the elment at the i_ row and j_ col
    //after assignment, let i_ and j_ return to the default value
    //erase all the memory
    i_ = -1;
    j_ = -1;
    i_write = -1;
    j_write = -1;
    index = 0;
    for(int i = 0; i < 5; ++i){
        ii[i] = -1;
        jj[i] = -1;
    }
    return *this;
}
sdMatrix sdMatrix::T (){
    int col = row_;
    int row = col_;
    sdMatrix transpose(0,row,col);
    for(int i = 0; i < row; ++i){
        for(int j = 0; j < col; ++j){
            transpose.entries_[i][j] = entries_[j][i];
        }
    }
    return transpose;
}
double sdMatrix::element(){
    
    //if (i_ == -1 || j_ == -1){
    if(index == 0){
        std::cout << "ELEMENT ERROR: Index is wrong" << std::endl;
        return 0;
    }
    else{
        index--;
        int i = ii[index];
        int j = jj[index];
        //double a = entries_[i_][j_];
        double a = entries_[i][j];
        //erase the memory
        ii[index + 1] = -1;
        jj[index + 1] = -1;
        i_ = -1;
        j_ = -1;
        return a;
    }
}
//set the name of matrix
void sdMatrix::setName(std::string name){
    name_ = name;//set the name
}
//output the length of the vector
int sdMatrix::length(){
    return row_;
}
//Output the matrix in double arrays
double** sdMatrix::toDouble(){
    double ** temp_ = new double*[row_];//must be dynamic allocation 
    for(int i = 0; i < row_ ; ++i){
        temp_[i] = new double[col_];
        for(int j = 0; j < col_; ++j){
            temp_[i][j] = entries_[i][j];
        }
    }
    return temp_;
}

