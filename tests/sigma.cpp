// this file is distributed under
// LGPLv3 license
#include <random>
#include <math.h>
#include <gtest/gtest.h>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
TEST(value, error_handling)
{
    EXPECT_NO_THROW(value<double>(1, -0.1));
    EXPECT_NO_THROW(value<double>(1, 0));
    EXPECT_NO_THROW(value<double>(1, 0.1));
    EXPECT_EQ(value<double>(1, 0).uncertainty(), 0);
    EXPECT_EQ(value<double>(1, 0.1).uncertainty(), 0.1);
    EXPECT_FALSE(isfinite(value<double>(1, -0.1).uncertainty()));
    EXPECT_NO_THROW(value<double>(1, -0.1) + value<double>(1, 0.1));
    EXPECT_EQ((value<double>(1, -0.1) + value<double>(1, 0.1)).val(), 2);
    EXPECT_FALSE(isfinite((value<double>(1, -0.1) + value<double>(1, 0.1)).uncertainty()));
}
TEST(value, fromlist)
{
    initializer_list<double> l = {};
    value<double> v;
    EXPECT_THROW(v = l, Exception<value<double>>);
    l = {1.0, 2.0, 3.0};
    EXPECT_THROW(v = l, Exception<value<double>>);
    l = {1.0, 2.0, 3.0, 4.0};
    EXPECT_THROW(v = l, Exception<value<double>>);
    l = {1.0};
    EXPECT_NO_THROW(v = l);
    l = {1.0, 2.0};
    EXPECT_NO_THROW(v = l);

    v = {1.0};
    EXPECT_EQ(1.0, v.val());
    EXPECT_EQ(0.0, v.uncertainty());
    v = {1.0, 2.0};
    EXPECT_EQ(1.0, v.val());
    EXPECT_EQ(2.0, v.uncertainty());
}
TEST(value, base)
{
    mt19937 gen;
    normal_distribution<double> G;
    for (size_t i = 0; i < 10; i++) {
        double x = G(gen);
        value<double> V(x, 0.1);
        EXPECT_EQ(x, V.val());
        EXPECT_EQ(0.1, V.uncertainty());
        EXPECT_EQ(0.1 / x, V.epsilon());
        EXPECT_EQ(x - 0.1, V.min());
        EXPECT_EQ(x + 0.1, V.max());
    }
}
TEST(value, compare1)
{
    value<double> V(0, 0.1);
    EXPECT_EQ(false, V.Contains(-0.2));
    EXPECT_EQ(true, V.Contains(-0.05));
    EXPECT_EQ(true, V.Contains(0));
    EXPECT_EQ(true, V.Contains(0.05));
    EXPECT_EQ(false, V.Contains(0.2));

    EXPECT_EQ(true, V.NotEqual(-0.2));
    EXPECT_EQ(false, V.NotEqual(-0.05));
    EXPECT_EQ(false, V.NotEqual(0));
    EXPECT_EQ(false, V.NotEqual(0.05));
    EXPECT_EQ(true, V.NotEqual(0.2));
}
TEST(value, compare2)
{
    bool v = value<double> {2, 1} .Above({0, 0.9});
    EXPECT_TRUE(v);
    v = value<double> {2, 1} .Above({0, 1.1});
    EXPECT_FALSE(v);
    v = value<double> {0, 1} .Below({2, 0.9});
    EXPECT_TRUE(v);
    v = value<double> {0, 1} .Below({2, 1.1});
    EXPECT_FALSE(v);

    v = value<double> {2, 1} .Above(0.9);
    EXPECT_TRUE(v);
    v = value<double> {2, 1} .Above(1.1);
    EXPECT_FALSE(v);
    v = value<double> {0, 1} .Below(1.1);
    EXPECT_TRUE(v);
    v = value<double> {0, 1} .Below(0.9);
    EXPECT_FALSE(v);
}
#define _EQ(a,b) EXPECT_TRUE(pow(a-b,2)<0.005)
TEST(value, arithmetic_actions1)
{
    mt19937 G;
    uniform_real_distribution<double> val(1, 50), unc(0.1, 10);
    for (size_t cnt = 0; cnt < 1000; cnt++) {
        value<double> A(val(G), unc(G)), B(val(G), unc(G));
        auto sum = A + B;
        EXPECT_EQ(A.val() + B.val(), sum.val());
        _EQ(pow(A.uncertainty(), 2) + pow(B.uncertainty(), 2), pow(sum.uncertainty(), 2));
        auto sub = A - B;
        EXPECT_EQ(A.val() - B.val(), sub.val());
        _EQ(pow(A.uncertainty(), 2) + pow(B.uncertainty(), 2), pow(sub.uncertainty(), 2));
        auto prod = A * B;
        EXPECT_EQ(A.val()*B.val(), prod.val());
        _EQ(pow(A.uncertainty()*B.val(), 2) + pow(A.val()*B.uncertainty(), 2), pow(prod.uncertainty(), 2));
        auto ratio = A / B;
        EXPECT_EQ(A.val() / B.val(), ratio.val());
        _EQ(pow(A.uncertainty() / B.val(), 2) + pow((A.val()*B.uncertainty()) / pow(B.val(), 2), 2), pow(ratio.uncertainty(), 2));
    }
}
TEST(value, arithmetic_actions2)
{
    mt19937 G;
    uniform_real_distribution<double> val(1, 50), unc(0.1, 10);
    for (size_t cnt = 0; cnt < 1000; cnt++) {
        value<double> A(val(G), unc(G));
        double B = val(G);
        auto sum = A + B;
        EXPECT_EQ(A.val() + B, sum.val());
        EXPECT_EQ(A.uncertainty(), sum.uncertainty());
        auto sub = A - B;
        EXPECT_EQ(A.val() - B, sub.val());
        EXPECT_EQ(A.uncertainty(), sub.uncertainty());
        auto prod = A * B;
        EXPECT_EQ(A.val()*B, prod.val());
        EXPECT_EQ(A.uncertainty()*B, prod.uncertainty());
        auto ratio = A / B;
        EXPECT_EQ(A.val() / B, ratio.val());
        EXPECT_EQ(A.uncertainty() / B, ratio.uncertainty());
    }
}
TEST(value, NumCompare)
{
    {
        double v = value<double>(0, 1).NumCompare(0);
        EXPECT_EQ(v, 0);
    }
    {
        double v = value<double>(0, 1).NumCompare(0.5);
        EXPECT_EQ(v, 0.25);
    }
    {
        double v = value<double>(0, 1).NumCompare(1);
        EXPECT_EQ(v, 1);
    }
    {
        double v = value<double>(0, 1).NumCompare(2);
        EXPECT_EQ(v, 4);
    }
    {
        double v = value<double>(2, 1).NumCompare(2);
        EXPECT_EQ(v, 0);
    }
    {
        double v = value<double>(2, 1).NumCompare(1.5);
        EXPECT_EQ(v, 0.25);
    }
    {
        double v = value<double>(2, 1).NumCompare(1);
        EXPECT_EQ(v, 1);
    }
    {
        double v = value<double>(2, 1).NumCompare(0);
        EXPECT_EQ(v, 4);
    }
    {
        double v = value<double>(0, 1).NumCompare({0, 1});
        EXPECT_EQ(v, 0);
    }
    {
        double v = value<double>(0, 1).NumCompare({1, 1});
        EXPECT_EQ(v, 0.25);
    }
    {
        double v = value<double>(0, 1).NumCompare({2, 1});
        EXPECT_EQ(v, 1);
    }
    {
        double v = value<double>(0, 1).NumCompare({4, 1});
        EXPECT_EQ(v, 4);
    }
    {
        double v = value<double>(2, 1).NumCompare({2, 1});
        EXPECT_EQ(v, 0);
    }
    {
        double v = value<double>(2, 1).NumCompare({1, 1});
        EXPECT_EQ(v, 0.25);
    }
    {
        double v = value<double>(2, 1).NumCompare({0, 1});
        EXPECT_EQ(v, 1);
    }
    {
        double v = value<double>(2, 1).NumCompare({ -2, 1});
        EXPECT_EQ(v, 4);
    }
}
TEST(value, func_val)
{
    value<> V(1, 0.1),
    V2([](double x)->double{return 2.0 * x;},V);
    EXPECT_EQ(2 * V.val(), V2.val());
    _EQ(2 * V.uncertainty(), V2.uncertainty());
}
TEST(value, wider)
{
    value<double> V1(1, 0.1), V2(2, 0.1), V3(1, 0.5);
    EXPECT_EQ(V1, V1.make_wider(1));
    EXPECT_EQ(V2, V2.make_wider(1));
    EXPECT_EQ(V3, V3.make_wider(1));
    auto V11 = V1.make_wider(0.5);
    auto V12 = V2.make_wider(0.5);
    auto V13 = V3.make_wider(0.5);
    EXPECT_EQ(V11.val(), V1.val());
    EXPECT_EQ(V12.val(), V2.val());
    EXPECT_EQ(V13.val(), V3.val());
    EXPECT_EQ(V11.uncertainty(), V1.uncertainty() * 0.5);
    EXPECT_EQ(V12.uncertainty(), V2.uncertainty() * 0.5);
    EXPECT_EQ(V13.uncertainty(), V3.uncertainty() * 0.5);
    auto V21 = V1.make_wider(2.5);
    auto V22 = V2.make_wider(2.5);
    auto V23 = V3.make_wider(2.5);
    EXPECT_EQ(V21.val(), V1.val());
    EXPECT_EQ(V22.val(), V2.val());
    EXPECT_EQ(V23.val(), V3.val());
    EXPECT_EQ(V21.uncertainty(), V1.uncertainty() * 2.5);
    EXPECT_EQ(V22.uncertainty(), V2.uncertainty() * 2.5);
    EXPECT_EQ(V23.uncertainty(), V3.uncertainty() * 2.5);
}
TEST(StandardDeviation, Throwing)
{
    StandardDeviation<double> S;
    EXPECT_EQ(0, S.count());
    EXPECT_THROW(S(), Exception<StandardDeviation<double>>);
    EXPECT_EQ(&S, &(S << 0.0));
    EXPECT_EQ(1, S.count());
    EXPECT_THROW(S(), Exception<StandardDeviation<double>>);
    EXPECT_EQ(&S, &(S << 0.0));
    EXPECT_EQ(2, S.count());
    EXPECT_EQ(0, S().val());
    EXPECT_EQ(0, S().uncertainty());
}
TEST(StandardDeviation, Base)
{
    StandardDeviation<double> S;
    S << 0.0;
    EXPECT_EQ(1, S.count());
    EXPECT_THROW(S(), Exception<StandardDeviation<double>>);
    S << 1.0;
    EXPECT_EQ(2, S.count());
    EXPECT_EQ(0.5, S().val());
    EXPECT_EQ(sqrt(0.5), S().uncertainty());
}
TEST(StandardDeviation, Base2)
{
    StandardDeviation<double> S;
    S << 0.0 << 1.0;
    size_t cnt=0;
    for(const auto&x:S)cnt++;
    EXPECT_EQ(2,cnt);
    EXPECT_EQ(2,S.size());
    EXPECT_EQ(2,S.count());
    EXPECT_EQ(0.0,S[0]);
    EXPECT_EQ(1.0,S[1]);
    EXPECT_ANY_THROW(S[2]);
}
#define _EQ2(a,b) EXPECT_TRUE(pow(a-b,2)<0.01)
TEST(StandardDeviation, WithRandomValues)
{
    StandardDeviation<double> S;
    default_random_engine generator;
    normal_distribution<double> distribution(1.0, 3.0);
    for (int i = 0; i < 2000; i++)
        S << distribution(generator);
    _EQ2(1.0, S().val());
    _EQ2(3.0, S().uncertainty());
}
TEST(StandardDeviation, WithRandomValues2)
{
    StandardDeviation<double> S(2);
    default_random_engine generator;
    normal_distribution<double> distribution(1.0, 3.0);
    for (int i = 0; i < 2000; i++)
        S << distribution(generator);
    _EQ2(1.0, S().val());
    _EQ2(6.0, S().uncertainty());
}

TEST(WeightedAverage, Zeros)
{
    WeightedAverage<double> W;
    EXPECT_THROW(W(), Exception<WeightedAverage<double>>);
    EXPECT_THROW(W << value<double>(0, 0), Exception<WeightedAverage<double>>);
    EXPECT_EQ(&W, &(W << value<double>(0, 1)));
    _EQ(0, W().val());
    _EQ(1, W().uncertainty());
    EXPECT_EQ(&W, &(W << value<double>(0, 1)));
    _EQ(0, W().val());
    _EQ(1.0 / sqrt(2.0), W().uncertainty());
    EXPECT_EQ(&W, &(W << value<double>(0, 1)));
    _EQ(0, W().val());
    _EQ(1.0 / sqrt(3.0), W().uncertainty());
}
TEST(WeightedAverage, Ones)
{
    WeightedAverage<double> W;
    EXPECT_THROW(W(), Exception<WeightedAverage<double>>);
    EXPECT_THROW(W << value<double>(1, 0), Exception<WeightedAverage<double>>);
    EXPECT_EQ(&W, &(W << value<double>(1, 1)));
    _EQ(1, W().val());
    _EQ(1, W().uncertainty());
    EXPECT_EQ(&W, &(W << value<double>(1, 1)));
    _EQ(1, W().val());
    _EQ(1.0 / sqrt(2.0), W().uncertainty());
    EXPECT_EQ(&W, &(W << value<double>(1, 1)));
    _EQ(1, W().val());
    _EQ(1.0 / sqrt(3.0), W().uncertainty());
}
TEST(WeightedAverage, Zeros_plus_Ones)
{
    WeightedAverage<double> W;
    EXPECT_THROW(W(), Exception<WeightedAverage<double>>);
    EXPECT_EQ(&W, &(W << value<double>(1, 1)));
    _EQ(1, W().val());
    _EQ(1, W().uncertainty());
    EXPECT_EQ(&W, &(W << value<double>(0, 1)));
    _EQ(0.5, W().val());
    _EQ(1.0 / sqrt(2.0), W().uncertainty());
    EXPECT_EQ(&W, &(W << value<double>(1, 1)));
    _EQ(2.0 / 3.0, W().val());
    _EQ(1.0 / sqrt(3.0), W().uncertainty());
    EXPECT_EQ(&W, &(W << value<double>(0, 1)));
    _EQ(0.5, W().val());
    _EQ(1.0 / sqrt(4.0), W().uncertainty());
}
TEST(CorrelationLinear, simple1)
{
    mt19937 gen;
    normal_distribution<double> G;
    CorrelationLinear<double> A, B, C;
    for (size_t i = 0; i < 10000; i++) {
        double v1 = G(gen), v2 = G(gen);
        EXPECT_EQ(&A, &(A << make_pair(v1, v1)));
        B << make_pair(v1, -v1);
        C << make_pair(v1, v2);
    }
    _EQ(A.R(), 1);
    _EQ(B.R(), -1);
    _EQ(C.R(), 0);
}
TEST(CorrelationLinear, simple2)
{
    mt19937 gen;
    normal_distribution<double> G;
    CorrelationLinear<double> A, B, C;
    for (size_t i = 0; i < 10000; i++) {
        double v1 = G(gen), v2 = G(gen);
        A << make_pair(v1   , v2);
        B << make_pair(v1 * 2., v2);
        C << make_pair(v1   , v2 * 2.);
    }
    _EQ(A.R(), 0);
    _EQ(B.R(), 0);
    _EQ(C.R(), 0);
}
TEST(CorrelationLinear, simple3)
{
    mt19937 gen;
    normal_distribution<double> G;
    CorrelationLinear<double> A, B;
    for (size_t i = 0; i < 10000; i++) {
        double v1 = G(gen), v2 = G(gen);
        A << make_pair(v1, v2);
        B << make_pair(v2, v1);
    }
    EXPECT_EQ(A.R(), B.R());
}
TEST(CorrelationLinear, simple4)
{
    mt19937 gen;
    normal_distribution<double> G;
    CorrelationLinear<double> A, B, C;
    for (size_t i = 0; i < 10000; i++) {
        double v1 = G(gen);
        A << make_pair(v1, v1);
        B << make_pair(v1 * 2., v1);
        C << make_pair(v1, v1 * 3.);
    }
    _EQ(A.R(), 1);
    _EQ(B.R(), 1);
    _EQ(C.R(), 1);
}
TEST(CorrelationLinear, simple5)
{
    mt19937 gen;
    normal_distribution<double> G;
    CorrelationLinear<double> A, B, C;
    for (size_t i = 0; i < 10000; i++) {
        double v1 = G(gen);
        A << make_pair(v1, -v1);
        B << make_pair(-v1 * 2., v1);
        C << make_pair(v1, -v1 * 3.);
    }
    _EQ(A.R(), -1);
    _EQ(B.R(), -1);
    _EQ(C.R(), -1);
}
