// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math.h>
#include <random>
#include <math_h/tabledata.h>
using namespace std;
using namespace MathTemplates;
double TestArray[] = { -2, -1, 0, 1, 2, 3, 4, 5, 6, 7};
TEST(WhereToInsert, BorderConditions)
{
    for (int beg = 0; beg < 10; beg++)
        for (double V = TestArray[beg] - 2; V <= TestArray[beg] + 2; V += 0.5) {
            int index = 0;
            if (beg > 0) {
                index = details::WhereToInsert(beg, beg - 1, TestArray, V);
                EXPECT_EQ(beg, index);
            }
            index = details::WhereToInsert(beg, beg, TestArray, V);
            if (V < TestArray[beg]) {
                EXPECT_EQ(beg, index);
            }
            if (V > TestArray[beg]) {
                EXPECT_EQ(beg + 1, index);
            }
            if (V == TestArray[beg]) {
                EXPECT_TRUE((index == beg) || (index == (beg + 1)));
            }
        }
}
TEST(WhereToInsert, NormalConditions)
{
    for (double x = TestArray[0] - 0.5; x <= TestArray[9] + 0.5; x += 0.5)
        for (int beg = 0; beg < 10; beg++)for (int end = beg + 1; end < 10; end++) {
                int index = details::WhereToInsert(beg, end, TestArray, x);
                if (x < TestArray[beg])EXPECT_EQ(beg, index);
                else if (x > TestArray[end])EXPECT_EQ(end + 1, index);
                else EXPECT_TRUE((x >= TestArray[index - 1]) && (x <= TestArray[index]));
            }
}
TEST(InsertSorted, BasicTest)
{
    Chain<int> X;
    for (int i = 0; i < 50; i++) {
        details::InsertSorted(rand() % 10, X, std_size(X), std_insert(X, int));
        for (int j = 0; j < i; j++)
            EXPECT_TRUE(X[j] <= X[j + 1]);
    }
}
TEST(SortedChain, basetest)
{
    mt19937 engine;
    SortedChain<double> chain;
    uniform_real_distribution<double> get(-10, 10);
    for (size_t i = 0; i < 100; i++)
        chain << get(engine);
    for (size_t i = 1; i < 100; i++)
        EXPECT_TRUE(chain[i] > chain[i - 1]);
}
