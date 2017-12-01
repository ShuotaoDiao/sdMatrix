//
//  main.cpp
//  Matrix_test01
//
//  Created by Shuotao Diao on 2/8/17.
//  Copyright © 2017 Shuotao Diao. All rights reserved.
//  This program is made to discover how to build matrix class

#include <iostream>
#include <string>
#include "sdMatrix.hpp"
sdMatrix f(sdMatrix, sdMatrix);

using namespace std;
int main(int argc, const char * argv[]) {
    
    double a[6] = {1, 2, 3, 4, 5, 6};
    double b[6] = {2, 4, 6, 8, 10, 12};
    sdMatrix A(a, 2, 3);
    //test on assigning the name to the matrix object
    A.setName("A");
    A.getValue();
    //test on one element assignment
    cout << "Test on one-element assignment" << endl;
    cout << "A(1,2) = 25" << endl;
    A(1,2) = 25;
    A.getValue();
    A(0,1) = 100;
    cout << "A(0,1) = 100" << endl;
    A.getValue();
    //test on robustness
    cout << "Test on robustness " << endl;
    cout << "A(100,20) = 99" << endl;
    A(100,20) = 99;//robustness test
    A.getValue();
    sdMatrix B(b, 3, 2);
    B.setName("B");
    B.getValue();
    sdMatrix C;
    C.setName("C");
    //test on matrix multiplication
    cout << "Test on matrix multiplication" << endl;
    cout << "C = A * B " << endl;
    C = A * B;
    A.getValue();
    B.getValue();
    C.getValue();
    sdMatrix D;
    D.setName("D");
    //test on addition
    cout << "Test on addition" << endl;
    cout << "D = 1 + A" << endl;
    D = 1 + A;
    A.getValue();
    D.getValue();
    D = D + A;
    cout << "D = D + A" << endl;
    D.getValue();
    cout << "D = D + B" << endl;
    B.getValue();
    D = D + B;
    D.getValue();
    sdMatrix E;
    E.setName("E");
    //test on subtraction
    cout << "Test on subtraction" << endl;
    cout << "E = 2 - B" << endl;
    E = 2 - B;
    B.getValue();
    E.getValue();
    //E.setName("EE");
    //E.getValue();
    //test on element by element multiplication
    cout << "Test on element by element multiplication" << endl;
    cout << "m3 = times(m1,m2)" << endl;
    sdMatrix m1(a,2,3);
    sdMatrix m2(b,2,3);
    m1.setName("m1");
    m1.getValue();
    m2.setName("m2");
    m2.getValue();
    sdMatrix m3 = times(m1, m2);
    m3.setName("m3");
    m3.getValue();
    //test on computing determinant
    cout << "Test on computing determinant" << endl;
    double c[9] = {-2,2,-3,-1,1,3,2,0,-1};
    sdMatrix M1("M1", c, 3, 3);
    double ans = M1.det();
    M1.getValue();
    cout << "The determinant of M1 is " << ans << endl;
    double d[9] = {0,0,1,0,1,0,1,0,0};
    sdMatrix M2("M2", d, 3, 3);
    M2.getValue();
    cout << "The determinant of M2 is " << M2.det() << endl;
    cout << "Test on robustness of computing determinant" << endl;
    cout << A.det() << endl;
    //test on inverting the matrix
    cout << "test on inverting the matrix " << endl;
    sdMatrix M1Inv = M1.inv();
    M1.getValue();
    cout << "The inverse of M1 is " << endl;
    M1Inv.getValue();
    cout << "The determinant of the inverse of M1 is " << M1Inv.det() << endl;
    sdMatrix M3 = M1 / M2;
    cout << "M3 = M1 / M2 is " << endl;
    M3.getValue();
    cout << "The determniant of the M3 is " << M3.det() << endl;
    double** MM3 = M3.toDouble();
    cout << MM3[0][0] << endl;
    MM3[0][0] = 2;
    cout << "MM3[0][0] = 2" << endl;
    cout << "MM3[0][0] = " << M3.toDouble()[0][2] << endl;
    M3.getValue();
    sdMatrix AA(a, 2, 3);
    sdMatrix BB(b, 2, 3);
    sdMatrix F;
    F = f(AA,BB);
    F.setName("F");
    F.getValue();
    cout << F(0,0).element() + F(0,1).element() << endl;
    (F(0,0)) = (F(0,0).element() + F(0,1).element());
    F.getValue();
    cout << "Hello, World!\n";
    
    return 0;
}

sdMatrix f(sdMatrix A, sdMatrix B){
    return A / 3+ 2 * B;
}
