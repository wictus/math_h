#include <iostream>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
//This is an example how to use
//values with uncertainties
int main()
{
    //how to declare values with uncertainties
    const value<> A(10.,1.0),B(3.,0.2);
    //how to get value's components
    cout<<"A: value = "<<A.val()<<endl;
    cout<<"A: delta = "<<A.uncertainty()<<endl;
    cout<<"A: epsilon = "<<A.epsilon()<<endl;
    cout<<"A: min = "<<A.min()<<endl;
    cout<<"A: max = "<<A.max()<<endl;
    cout<<"B: value = "<<B.val()<<endl;
    cout<<"B: delta = "<<B.uncertainty()<<endl;
    cout<<"B: epsilon = "<<B.epsilon()<<endl;
    cout<<"B: min = "<<B.min()<<endl;
    cout<<"B: max = "<<B.max()<<endl;
    //how to perform arithmetic operations
    cout<<"A+B = "<<A+B<<endl;
    cout<<"A-B = "<<A-B<<endl;
    cout<<"A*B = "<<A*B<<endl;
    cout<<"A/B = "<<A/B<<endl;
    const auto func=[](double x)->double{return log(x);};
    cout<<"log(A) = "<< value<>(func,A) <<endl;
}
