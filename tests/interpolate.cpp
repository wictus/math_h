// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math.h>
#include <math_h/interpolate.h>
using namespace std;
using namespace MathTemplates;
TEST(LinearInterpolation, Create)
{
    LinearInterpolation<double> F;
    EXPECT_EQ(0, F.size());
    EXPECT_THROW(F.left().X(), Exception<SortedChain<point<>>>);
    EXPECT_THROW(F.right().X(), Exception<SortedChain<point<>>>);
    EXPECT_EQ(&F, &(F << point<double>(0, 0)));
    EXPECT_EQ(1, F.size());
    EXPECT_NO_THROW(F.left().X());
    EXPECT_NO_THROW(F.right().X());
    EXPECT_THROW(F(0.5), Exception<LinearInterpolation<>>);
    EXPECT_EQ(&F, &(F << point<double>(1, 0)));
    EXPECT_EQ(2, F.size());
    EXPECT_NO_THROW(F.left().X());
    EXPECT_NO_THROW(F.right().X());
    EXPECT_EQ(0, F(0.5));
}
#define _EQ(a,b) EXPECT_TRUE(pow(a-b,2)<0.0001)
TEST(LinearInterpolation, SimpleLine)
{
    LinearInterpolation<double> F;
    F << point<double>(0, 0) << point<double>(1, 1);
    EXPECT_EQ(0, F(0));
    EXPECT_EQ(1, F(1));
    _EQ(0.1, F(0.1));
    _EQ(0.2, F(0.2));
    _EQ(0.3, F(0.3));
    _EQ(0.4, F(0.4));
    _EQ(0.5, F(0.5));
    _EQ(0.6, F(0.6));
    _EQ(0.7, F(0.7));
    _EQ(0.8, F(0.8));
    _EQ(0.9, F(0.9));
    F << point<double>(2, 0);
    EXPECT_EQ(0, F(0));
    EXPECT_EQ(1, F(1));
    _EQ(0.1, F(0.1));
    _EQ(0.2, F(0.2));
    _EQ(0.3, F(0.3));
    _EQ(0.4, F(0.4));
    _EQ(0.5, F(0.5));
    _EQ(0.6, F(0.6));
    _EQ(0.7, F(0.7));
    _EQ(0.8, F(0.8));
    _EQ(0.9, F(0.9));
    EXPECT_EQ(0, F(2));
    _EQ(0.1, F(1.9));
    _EQ(0.2, F(1.8));
    _EQ(0.3, F(1.7));
    _EQ(0.4, F(1.6));
    _EQ(0.5, F(1.5));
    _EQ(0.6, F(1.4));
    _EQ(0.7, F(1.3));
    _EQ(0.8, F(1.2));
    _EQ(0.9, F(1.1));
    F << point<double>(0.5, 0.5) << point<double>(1.5, 0.5);
    EXPECT_EQ(0, F(0));
    EXPECT_EQ(1, F(1));
    _EQ(0.1, F(0.1));
    _EQ(0.2, F(0.2));
    _EQ(0.3, F(0.3));
    _EQ(0.4, F(0.4));
    EXPECT_EQ(0.5, F(0.5));
    _EQ(0.6, F(0.6));
    _EQ(0.7, F(0.7));
    _EQ(0.8, F(0.8));
    _EQ(0.9, F(0.9));
    EXPECT_EQ(0, F(2));
    _EQ(0.1, F(1.9));
    _EQ(0.2, F(1.8));
    _EQ(0.3, F(1.7));
    _EQ(0.4, F(1.6));
    EXPECT_EQ(0.5, F(1.5));
    _EQ(0.6, F(1.4));
    _EQ(0.7, F(1.3));
    _EQ(0.8, F(1.2));
    _EQ(0.9, F(1.1));
    for (double x = 0; x < 2; x += 0.1)EXPECT_EQ(F(x), F.func()(x));
}
