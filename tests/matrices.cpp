// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math.h>
#include <math_h/matrices.h>
using namespace std;
using namespace MathTemplates;
typedef Exception<Chain<Chain<>>, 0> invalid_matr;
typedef Exception<Chain<Chain<>>, 1> size_mismatch;
TEST(Matrix, SimpleObjects)
{
    auto A = Unitary<double>(3);
    EXPECT_EQ(A.height(), A.width());
    EXPECT_EQ(3, A.height());
    EXPECT_TRUE(A == A);
    EXPECT_FALSE(A != A);
    EXPECT_TRUE(A.HasSizeAs(A));
    EXPECT_TRUE(A == Unitary<double>(3));
    EXPECT_FALSE(Unitary<double>(3) != A);
    EXPECT_TRUE(Unitary<double>(3) == A);
    EXPECT_FALSE(A != Unitary<double>(3));
    EXPECT_TRUE(A.HasSizeAs(Unitary<double>(3)));
    for (size_t i = 0; i < 3; i++)for (size_t j = 0; j < 3; j++)
            if (i == j)EXPECT_EQ(A(i, j), 1);
            else EXPECT_EQ(A(i, j), 0);
    EXPECT_TRUE(A.HasSizeAs(Zeros<double>(3, 3)));
    EXPECT_FALSE(A == Zeros<double>(3, 3));
    EXPECT_TRUE(Zeros<double>(3, 3).HasSizeAs(A));
    EXPECT_FALSE(Zeros<double>(3, 3) == A);
    EXPECT_EQ(Zeros<double>(4, 2).height(), 4);
    EXPECT_EQ(Zeros<double>(4, 2).width(), 2);
    EXPECT_TRUE(Zeros<double>(4, 2) == MatrixByFormula<double>(4, 2, [](size_t, size_t) {
        return 0.0;
    }));
    auto Z = Zeros<double>(3, 4);
    EXPECT_TRUE(Z == Z);
    EXPECT_FALSE(Z != Z);
    EXPECT_EQ(Z.height(), 3);
    EXPECT_EQ(Z.width(), 4);
    for (size_t i = 0; i < 3; i++)for (size_t j = 0; j < 4; j++)
            EXPECT_EQ(Z(i, j), 0);
}
const function<double(size_t, size_t)> F = [](size_t i, size_t j)
{
    return i * 30 + j;
};
TEST(Matrix, Formula)
{
    auto A = MatrixByFormula<double>(3, 4, F);
    auto B = MatrixData<double>({
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}
    });
    B.Transform([](size_t i, size_t j, const double &) {
        return F(i, j);
    });
    EXPECT_TRUE(A == B);
    for (size_t i = 0; i < 3; i++)for (size_t j = 0; j < 4; j++)
            EXPECT_EQ(A(i, j), F(i, j));
}
TEST(Matrix, Transponate)
{
    auto A = MatrixByFormula<double>(3, 4, F);
    auto B = Transponate(A);
    EXPECT_EQ(A.width(), B.height());
    EXPECT_EQ(B.width(), A.height());
    for (size_t i = 0; i < A.height(); i++)for (size_t j = 0; j < A.width(); j++)
            EXPECT_EQ(A(i, j), B(j, i));
    EXPECT_TRUE(Transponate(A) == B);
}
TEST(Matrix, Multiply)
{
    EXPECT_EQ(Multiply(Zeros<double>(3, 3), Zeros<double>(3, 3)), Zeros<double>(3, 3));
    EXPECT_EQ(Multiply(Zeros<double>(3, 5), Zeros<double>(5, 3)), Zeros<double>(3, 3));
    EXPECT_EQ(Multiply(Zeros<double>(2, 3), Zeros<double>(3, 4)), Zeros<double>(2, 4));
    EXPECT_EQ(Multiply(Unitary<double>(1), Unitary<double>(1)), Unitary<double>(1));
    EXPECT_EQ(Multiply(Unitary<double>(2), Unitary<double>(2)), Unitary<double>(2));
    EXPECT_EQ(Multiply(Unitary<double>(3), Unitary<double>(3)), Unitary<double>(3));
    EXPECT_EQ(Multiply(Unitary<double>(4), Unitary<double>(4)), Unitary<double>(4));
    EXPECT_EQ(Multiply(Unitary<double>(5), Unitary<double>(5)), Unitary<double>(5));
    EXPECT_EQ(MatrixByFormula<double>(4, 4, F), Multiply(MatrixByFormula<double>(4, 4, F), Unitary<double>(4)));
    EXPECT_EQ(MatrixByFormula<double>(4, 4, F), Multiply(Unitary<double>(4), MatrixByFormula<double>(4, 4, F)));
    EXPECT_EQ(Zeros<double>(4, 4), Multiply(MatrixByFormula<double>(4, 4, F), Zeros<double>(4, 4)));
    EXPECT_EQ(Zeros<double>(4, 4), Multiply(Zeros<double>(4, 4), MatrixByFormula<double>(4, 4, F)));
    for (size_t i = 0; i < 5; i++)for (size_t j = 0; j < 5; j++)
            EXPECT_EQ(Multiply(Transponate(RVec<double>(5, i)), RVec<double>(5, j)), (i == j) ? 1.0 : 0.0);
    auto A = MatrixData<double>({
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 0, 1, 2}
    });
    auto B = MatrixData<double>({
        {1, 2},
        {3, 4},
        {5, 6},
        {7, 8}
    });
    auto C = MatrixData<double>({
        {50, 60},
        {114, 140},
        {28, 40}
    });
    EXPECT_EQ(Multiply(A, B), C);
}
TEST(Matrix, Diagonal)
{
    EXPECT_EQ(
        Diagonal<double>({1, 2, 3, 4}),
    MatrixData<double>({
        {1, 0, 0, 0},
        {0, 2, 0, 0},
        {0, 0, 3, 0},
        {0, 0, 0, 4}
    })
    );
}
TEST(Matrix, Permutation)
{
    EXPECT_THROW(Permutation<double>({3, 3}), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(Permutation<double>({0, 3}), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(Permutation<double>({3, 0}), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(Permutation<double>({1, 3}), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(Permutation<double>({2, 2}), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(Permutation<double>({0, 2}), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(Permutation<double>({2, 0}), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(Permutation<double>({1, 2}), Exception<MatrixByFormula<double>>);
    //ToDo: provide such control
    //EXPECT_THROW(Permutation<double>({1,1}),Exception<MatrixByFormula<double>>);
    //EXPECT_THROW(Permutation<double>({0,0}),Exception<MatrixByFormula<double>>);
    //EXPECT_THROW(Permutation<double>({1,1,0}),Exception<MatrixByFormula<double>>);
    //EXPECT_THROW(Permutation<double>({2,1,1}),Exception<MatrixByFormula<double>>);
    //EXPECT_THROW(Permutation<double>({1,0,1}),Exception<MatrixByFormula<double>>);
    EXPECT_EQ(Permutation<double>({0, 1}), MatrixData<double>({{1, 0}, {0, 1}}));
    EXPECT_EQ(Permutation<double>({1, 0}), MatrixData<double>({{0, 1}, {1, 0}}));
    EXPECT_EQ(Permutation<double>({0, 1, 2}), MatrixData<double>({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}));
    EXPECT_EQ(Permutation<double>({1, 0, 2}), MatrixData<double>({{0, 1, 0}, {1, 0, 0}, {0, 0, 1}}));
    EXPECT_EQ(Permutation<double>({2, 0, 1}), MatrixData<double>({{0, 0, 1}, {1, 0, 0}, {0, 1, 0}}));
}
TEST(Matrix, Minor)
{
    auto A = MatrixData<double>({
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    });
    EXPECT_THROW(Minor(A, 3, 3), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(Minor(A, 4, 3), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(Minor(A, 3, 4), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(Minor(A, 4, 4), Exception<MatrixByFormula<double>>);
    EXPECT_EQ(Minor(A, 0, 0), MatrixData<double>({{5, 6}, {8, 9}}));
    EXPECT_EQ(Minor(A, 0, 1), MatrixData<double>({{4, 6}, {7, 9}}));
    EXPECT_EQ(Minor(A, 0, 2), MatrixData<double>({{4, 5}, {7, 8}}));
    EXPECT_EQ(Minor(A, 1, 0), MatrixData<double>({{2, 3}, {8, 9}}));
    EXPECT_EQ(Minor(A, 1, 1), MatrixData<double>({{1, 3}, {7, 9}}));
    EXPECT_EQ(Minor(A, 1, 2), MatrixData<double>({{1, 2}, {7, 8}}));
    EXPECT_EQ(Minor(A, 2, 0), MatrixData<double>({{2, 3}, {5, 6}}));
    EXPECT_EQ(Minor(A, 2, 1), MatrixData<double>({{1, 3}, {4, 6}}));
    EXPECT_EQ(Minor(A, 2, 2), MatrixData<double>({{1, 2}, {4, 5}}));
    EXPECT_EQ(Minor(Minor(A, 0, 0), 1, 1), 5.0);
    EXPECT_EQ(Minor(Minor(A, 0, 0), 0, 0), 9.0);
    EXPECT_EQ(Minor(Minor(A, 2, 1), 1, 1), 1.0);
    EXPECT_EQ(Minor(Minor(A, 2, 1), 0, 0), 6.0);
    EXPECT_EQ(Minor(Minor(A, 2, 1), 0, 1), 4.0);
    EXPECT_EQ(Minor(Minor(A, 1, 1), 1, 1), 1.0);
    EXPECT_THROW(Minor(Minor(Minor(A, 0, 0), 1, 1), 0, 0), Exception<MatrixByFormula<double>>);
}
TEST(Matrix, Determinant)
{
    EXPECT_THROW(Determinant(MatrixData<double>({{1, 2}})), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(Determinant(MatrixData<double>({{1}, {2}})), Exception<MatrixByFormula<double>>);
    EXPECT_EQ(Determinant(MatrixData<double>(5)), 5.0);
    EXPECT_EQ(Determinant(MatrixData<double>({{1, 2}, {3, 4}})), -2);
    EXPECT_EQ(Determinant(MatrixData<double>({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}})), 0);
    EXPECT_EQ(Determinant(Unitary<double>(1)), 1);
    EXPECT_EQ(Determinant(Unitary<double>(2)), 1);
    EXPECT_EQ(Determinant(Unitary<double>(3)), 1);
    EXPECT_EQ(Determinant(Unitary<double>(4)), 1);
    EXPECT_EQ(Determinant(Unitary<double>(5)), 1);
    EXPECT_EQ(Determinant(Zeros<double>(1, 1)), 0);
    EXPECT_EQ(Determinant(Zeros<double>(2, 2)), 0);
    EXPECT_EQ(Determinant(Zeros<double>(3, 3)), 0);
    EXPECT_EQ(Determinant(Zeros<double>(4, 4)), 0);
    EXPECT_EQ(Determinant(Zeros<double>(5, 5)), 0);
}
TEST(Matrix, Inverse)
{
    EXPECT_THROW(CalcInverseMatrix(Zeros<double>(1, 3)), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(CalcInverseMatrix(Zeros<double>(1, 2)), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(CalcInverseMatrix(Zeros<double>(2, 2)), Exception<MatrixByFormula<double>>);
    EXPECT_THROW(CalcInverseMatrix(Zeros<double>(3, 3)), Exception<MatrixByFormula<double>>);

    EXPECT_EQ(CalcInverseMatrix(Unitary<double>(1)), Unitary<double>(1));
    EXPECT_EQ(CalcInverseMatrix(Unitary<double>(2)), Unitary<double>(2));
    EXPECT_EQ(CalcInverseMatrix(Unitary<double>(3)), Unitary<double>(3));
    EXPECT_EQ(CalcInverseMatrix(Unitary<double>(4)), Unitary<double>(4));

    EXPECT_EQ(CalcInverseMatrix(MatrixData<double>(2)), MatrixData<double>(0.5));
    EXPECT_EQ(CalcInverseMatrix(MatrixData<double>(10)), MatrixData<double>(0.1));
    {
        MatrixData<double> A({{1, 2}, {3, 4}});
        EXPECT_EQ(Multiply(A, CalcInverseMatrix(A)), Unitary<double>(A.height()));
        EXPECT_EQ(Multiply(CalcInverseMatrix(A), A), Unitary<double>(A.height()));
    }
    {
        MatrixData<double>A({
            {2, 5, 7},
            {6, 3, 4},
            {5, -2, -3}
        });
        EXPECT_EQ(Multiply(A, CalcInverseMatrix(A)), Unitary<double>(A.height()));
        EXPECT_EQ(Multiply(CalcInverseMatrix(A), A), Unitary<double>(A.height()));
    }
}
TEST(Matrix, Solve)
{
    {
        const MatrixData<double> A({{1, 2, 0}, {3, 4, 4}, {5, 6, 3}}), b({{3}, {7}, {8}}),
        x = Solve(A, b), b2 = Multiply(A, x);
        ASSERT_EQ(b.height(), b2.height());
        for (size_t i = 0; i < b.height(); i++)
            EXPECT_TRUE(pow(b(i, 0) - b2(i, 0), 2) < 0.0000001);
    }{
        const MatrixData<double> A({{2, 0, 2, 0.6}, {3, 3, 4, -2}, {5, 5, 4, 2}, { -1, -2, 3.4, -1}}),
        b({{1}, {2}, {3}, {1}}), x = Solve(A, b), b2 = Multiply(A, x);
        ASSERT_EQ(b.height(), b2.height());
        for (size_t i = 0; i < b.height(); i++)
            EXPECT_TRUE(pow(b(i, 0) - b2(i, 0), 2) < 0.0000001);
    }
}
