// this file is distributed under
// LGPLv3 license
#include <iostream>
#include <gtest/gtest.h>
#include <math_h/integrate.h>
#include <math_h/functions.h>
using namespace std;
using namespace MathTemplates;
#define _EQ(a,b) EXPECT_TRUE(pow(a-b,2)<0.00001)
TEST(Sympson, BaseTest)
{
    auto F = [](double x) {
        return x;
    };
    _EQ(0.5, Sympson(F, 0.0, 1.0, 0.1));
    _EQ(0.5, Sympson(F, 0.0, 1.0, -0.1));
    _EQ(-0.5, Sympson(F, 1.0, 0.0, 0.1));
    _EQ(-0.5, Sympson(F, 1.0, 0.0, -0.1));
}
TEST(Sympson, BaseTest2)
{
    _EQ(0.0, Sympson([](double) {
        return 0.0;
    }, 0.0, 1.0, 0.00001));
    _EQ(1.0, Sympson([](double) {
        return 1.0;
    }, 0.0, 1.0, 0.00001));
    _EQ(0.5, Sympson([](double x) {
        return x;
    }, 0.0, 1.0, 0.00001));
    _EQ(0.333333, Sympson([](double x) {
        return x * x;
    }, 0.0, 1.0, 0.00001));
    _EQ(0.25, Sympson([](double x) {
        return x * x * x;
    }, 0.0, 1.0, 0.00001));
    _EQ(1.0, Sympson([](double x) {
        return Gaussian(x, 5.0, 1.0);
    }, 0.0, 10.0, 0.0001));
}
TEST(Int_Trapez_Table, BasicTest)
{
    SortedPoints<double> func;
    for (double x = 0; x <= 1; x += 0.0001)
        func << point<double>(x, x * x);
    auto res = Int_Trapez_Table(func);
    EXPECT_EQ(func.left().X(), res.left().X());
    EXPECT_EQ(func.right().X(), res.right().X());
    EXPECT_EQ(0, res.left().Y());
    _EQ(0.33333333, res.right().Y());
    cout << res.right().Y() << endl;
}
TEST(Int_Trapez_Table_PositiveStrict, BasicTest)
{
    SortedPoints<double> func;
    for (double x = 0; x <= 1; x += 0.0001)
        func << point<double>(x, x * x);
    auto res = Int_Trapez_Table_PositiveStrict(func);
    EXPECT_EQ(func.left().X(), res.left().X());
    EXPECT_EQ(func.right().X(), res.right().X());
    EXPECT_EQ(0, res.left().Y());
    _EQ(0.333333, res.right().Y());
    cout << res.right().Y() << endl;
}
TEST(Convolution, BasicTest)
{
    auto F1 = [](double x) {
        return x;
    };
    auto F2 = [](double x) {
        return x * x;
    };
    Convolution<double> C(F1, F2, 0, 1, 0.01);
    double x = 0.1;
    double test_value = Sympson([x, F1, F2](double k) {
        return F1(k) * F2(x - k);
    }, 0.0, 1.0, 0.01);
    EXPECT_EQ(test_value, C(x));
}
TEST(Convolution, BasicTest2)
{
    auto F1 = [](const double & x) {
        return x;
    };
    auto F2 = [](const double & x) {
        return x * x;
    };
    auto C = make_convolution(F1, F2, 0., 1., 0.01);
    double x = 0.1;
    double test_value = Sympson([x, F1, F2](double k) {
        return F1(k) * F2(x - k);
    }, 0.0, 1.0, 0.01);
    EXPECT_EQ(test_value, C(x));
}
