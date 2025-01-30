#include <iostream>
#include "MatVec.h"

int main(int argc, char **argv)
{
    // Test initialization
    CVector myVec(5, 5.0);

    // Test copy
    CVector myVec2(myVec);

    // Test initiate
    std::cout << "initiate()" << std::endl;
    CVector myVec4;
    myVec4.Initialize(5, 7.0);
    myVec4.print();

    // Test getSize()
    std::cout << "getSize()" << std::endl;
    std::cout << myVec.GetSize() << std::endl;
    std::cout << myVec2.GetSize() << std::endl;

    // Test dotProduct()
    std::cout << "dot()" << std::endl;
    CVector myVec3(5, 1.0);
    double dotProduct(0.0), dotProduct2(0.0);
    dotProduct = myVec.dotProd(myVec3);
    dotProduct2 = myVec3.dotProd(myVec);
    std::cout << dotProduct << std::endl;
    std::cout << dotProduct2 << std::endl;

    // Test norm()
    std::cout << "norm()" << std::endl;
    std::cout << myVec.norm() << std::endl;

    // Test display()
    std::cout << "display()" << std::endl;
    myVec.print();

    // Test reset();
    std::cout << "reset()" << std::endl;
    myVec2.Reset();
    myVec2.print();

    // Test operator[] as member
    std::cout << "operator[]" << std::endl;
    std::cout << myVec[3] << std::endl;
    myVec[3] = 100.0;
    myVec.print();

    // Test operator = as member
    std::cout << "operator=" << std::endl;
    myVec2 = myVec;
    myVec2.print();

    // Test operator += as member
    std::cout << "operator+=" << std::endl;
    myVec.SetAllValues(2.0);
    myVec2.SetAllValues(3.0);
    myVec2 += myVec;
    myVec2.print();

    // Test operator -= as member
    std::cout << "operator-=" << std::endl;
    myVec2 -= myVec;
    myVec2.print();

    // Test operator *= as member
    std::cout << "operator*=" << std::endl;
    myVec.SetAllValues(2.0);
    myVec *= 2.0;
    myVec.print();

    // Test operator /= as member
    std::cout << "operator/=" << std::endl;
    myVec /= 2.0;
    myVec.print();

    // Test addition
    std::cout << "addition" << std::endl;
    myVec.SetAllValues(2.0);
    myVec2.SetAllValues(3.0);
    CVector vecSum(5, 0.0);
    vecSum = myVec + myVec2;
    vecSum.print();

    // Test multiplication
    std::cout << "multiplication" << std::endl;
    myVec.SetAllValues(2.0);
    CVector vecMul(5, 0.0);
    vecMul = myVec * 2.0;
    vecMul = 3.0 * myVec;
    vecMul.print();

    // Test linear combination
    std::cout << "linear combination" << std::endl;
    myVec.SetAllValues(1.0);
    myVec2.SetAllValues(1.0);
    myVec3.SetAllValues(1.0);
    myVec4 = myVec + (2.0 * myVec2 + 5.0 * myVec3) * 2.0;
    myVec4.print();

    std::cout << "Pointers..." << std::endl;
    CVector *myVecPoint = NULL;
    CVector *myVecPoint2 = NULL;
    myVecPoint = new CVector(5, 10.0);
    myVecPoint->print();
    std::cout << (*myVecPoint)[3] << std::endl;
    myVecPoint2 = new CVector(*myVecPoint);
    myVecPoint2->print();
    myVecPoint->Reset();
    *myVecPoint2 = *myVecPoint;
    myVecPoint2->print();

    if (myVecPoint != NULL)
        delete myVecPoint;
    if (myVecPoint2 != NULL)
        delete myVecPoint2;

    return 0;
}
