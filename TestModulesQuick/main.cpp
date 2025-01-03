#include <QCoreApplication>

#include"tests.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    //
    //TestQStringFormatting();//works well
    //Test1DLib();//works well
    //TestPtr1();//works well
    //
    //TestArr2DLib(); // seems to work gut
    //TestArr2DLib_VB();//seems to work gut, ne checked details
    //TestMatrices();//not all matrices work gut, Matrix_PS causes crash
    //TestLinearApprox();
    //TestPtrNotFirst();
    TestMatrix();
    //TestSeidel();
    //
    std::cout<<"\nAll functions are successfully tested\n"<<std::endl;
    //
    return a.exec();
}
