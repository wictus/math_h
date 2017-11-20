// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math.h>
#include <random>
#include <math_h/tabledata.h>
using namespace std;
using namespace MathTemplates;
TEST(SortedPoints, size1)
{
    SortedPoints<double> chain;
    EXPECT_EQ(0, chain.size());
    chain << point<double>(0, 1);
    EXPECT_EQ(1, chain.size());
    chain << point<double>(2, 1);
    EXPECT_EQ(2, chain.size());
    chain << point<double>(1, 1);
    EXPECT_EQ(3, chain.size());
    chain << point<double>(3, 1);
    EXPECT_EQ(4, chain.size());
    chain << point<double>(5, 1);
    EXPECT_EQ(5, chain.size());
    chain << point<double>(4, 1);
    EXPECT_EQ(6, chain.size());
    EXPECT_EQ(0, chain.left().X());
    EXPECT_EQ(5, chain.right().X());
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(1, chain[i].Y());
        EXPECT_EQ(i, chain[i].X());
    }
}
TEST(SortedPoints, range)
{
    SortedPoints<> chain=Points<>{{0, 1},{2, 1.1},{3, 1.2},{4, 1.3},{5, 1.4}};
    auto xcut = chain.XRange(0.5, 4.5);
    EXPECT_EQ(3, xcut.size());
    EXPECT_EQ(2, xcut.left().X());
    EXPECT_EQ(4, xcut.right().X());
    auto ycut = chain.YRange(1.15, 1.35).TransponateAndSort();
    EXPECT_EQ(2, ycut.size());
    EXPECT_EQ(1.2, ycut.left().X());
    EXPECT_EQ(1.3, ycut.right().X());
}
#include <math_h/sigma.h>
TEST(SortedPoints, size2)
{
    SortedPoints<value<double>> chain;
    EXPECT_EQ(0, chain.size());
    chain << point<value<double>>(0, 1);
    EXPECT_EQ(1, chain.size());
    chain << point<value<double>>(2, 1);
    EXPECT_EQ(2, chain.size());
    chain << point<value<double>>(1, 1);
    EXPECT_EQ(3, chain.size());
    chain << point<value<double>>(3, 1);
    EXPECT_EQ(4, chain.size());
    chain << point<value<double>>(5, 1);
    EXPECT_EQ(5, chain.size());
    chain << point<value<double>>(4, 1);
    EXPECT_EQ(6, chain.size());
    EXPECT_EQ(0, chain.left().X().val());
    EXPECT_EQ(5, chain.right().X().val());
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(1, chain[i].Y().val());
        EXPECT_EQ(i, chain[i].X().val());
    }
}
