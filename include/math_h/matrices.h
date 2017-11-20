// this file is distributed under
// LGPLv3 license
#ifndef ______MATRICES_H_______
#	define ______MATRICES_H_______
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <iostream>
#include <math.h>
#include "error.h"
#include "tabledata.h"
namespace MathTemplates
{
template<class numt = double>
class Matrix
{
public:
    typedef std::function<const numt(const size_t, const size_t)> function;
    virtual ~Matrix() {}
    virtual const size_t height()const = 0;
    virtual const size_t width()const = 0;
protected:
    virtual const numt get_element(const size_t i, const size_t j)const = 0;
public:
    const numt operator()(const size_t i, const size_t j)const
    {
        if (i >= height())throw Exception<Matrix>("Range check error");
        if (j >= width())throw Exception<Matrix>("Range check error");
        return get_element(i, j);
    }
    inline const bool HasSizeAs(const Matrix &other)const
    {
        return (height() == other.height()) && (width() == other.width());
    }
    const bool operator==(const Matrix &other)const
    {
        if (!HasSizeAs(other))return false;
        for (size_t i = 0; i < height(); i++) {
            for (size_t j = 0; j < width(); j++)
                if (operator()(i, j) != other(i, j))
                    return false;
        }
        return true;
    }
    inline const bool operator!=(const Matrix &other)const
    {
        return !operator==(other);
    }
    inline const bool operator==(const numt &c)const
    {
        return (width() == 1) && (height() == 1) && (operator()(0, 0) == c);
    }
    inline const bool operator!=(const numt &c)const
    {
        return !operator==(c);
    }
};
template<typename numt>
inline std::ostream &operator<<(std::ostream &str, const Matrix<numt> &M)
{
    for (size_t i = 0, n = M.height(); i < n; i++) {
        for (size_t j = 0, m = M.width(); j < m; j++)
            str << M(i, j) << "\t";
        str << std::endl;
    }
    return str;
}

template<class numt = double>
class MatrixData: public Matrix<numt>
{
public:
    typedef Chain<Chain<numt>> container;
    typedef std::function<const numt(const size_t, const size_t, const numt &)> transform_function;
    typedef std::function<const numt(const numt &)> transform_func;
private:
    container f_data;
protected:
    virtual const numt get_element(const size_t i, const size_t j)const override
    {
        return f_data[i][j];
    }
public:
    virtual const size_t height()const override
    {
        return f_data.size();
    }
    virtual const size_t width()const override
    {
        return f_data[0].size();
    }
    numt &var_element(const size_t i, const size_t j)
    {
        return f_data[i][j];
    }
    virtual ~MatrixData() {}
    MatrixData(const numt &A)
    {
        f_data.push_back({A});
    }
    MatrixData(const container &A)
    {
        if (A.size() == 0)
            throw Exception<MatrixData, 0>("invalid matrix size");
        if ((A.size() == 1) && (A[0].size() == 0))
            throw Exception<MatrixData, 0>("invalid matrix size");
        const size_t first_row_size = A[0].size();
        if (first_row_size == 0)
            throw Exception<MatrixData, 0>("invalid matrix size");
        for (const auto &row : A) {
            if (row.size() != first_row_size)
                throw Exception<MatrixData, 0>("invalid matrix size");
            f_data.push_back(Chain<numt>());
            const auto row_index = f_data.size() - 1;
            for (const auto &item : row)
                f_data[row_index].push_back(item);
        }
    }
    MatrixData(const Matrix<numt> &source)
    {
        for (size_t i = 0; i < source.height(); i++) {
            f_data.push_back(Chain<numt>());
            for (size_t j = 0; j < source.width(); j++)
                f_data[i].push_back(source(i, j));
        }
    }
    MatrixData &Transform(const transform_function F)
    {
        for (size_t i = 0; i < height(); i++) {
            for (size_t j = 0; j < width(); j++)
                f_data[i][j] = F(i, j, f_data[i][j]);
        }
        return *this;
    }
    MatrixData &Transform(const transform_func F)
    {
        for (size_t i = 0; i < height(); i++) {
            for (size_t j = 0; j < width(); j++)
                f_data[i][j] = F(f_data[i][j]);
        }
        return *this;
    }
};
template<class numt = double>
class MatrixByFormula: public Matrix<numt>
{
private:
    size_t f_N, f_M;
    typename Matrix<numt>::function f_func;
public:
    MatrixByFormula(const size_t N, const size_t M, const typename Matrix<numt>::function F): f_N(N), f_M(M), f_func(F) {}
    virtual ~MatrixByFormula() {}
    virtual const size_t height()const override
    {
        return f_N;
    }
    virtual const size_t width()const override
    {
        return f_M;
    }
protected:
    virtual const numt get_element(const size_t i, const size_t j)const override
    {
        return f_func(i, j);
    }
};

template<class numt>
const MatrixByFormula<numt> Unitary(const size_t N)
{
    return MatrixByFormula<numt>(N, N, [](size_t i, size_t j)->numt {return (i == j) ? 1 : 0;});
}
template<class numt>
const MatrixByFormula<numt> Zeros(const size_t N, const size_t M)
{
    return MatrixByFormula<numt>(N, M, [](size_t, size_t)->numt {return 0;});
}
template<class numt>
const MatrixByFormula<numt> RVec(const size_t N, const size_t i)
{
    if (i >= N)throw Exception<Matrix<numt>>("Invalid reference vector");
    return MatrixByFormula<numt>(N, 1, [i](size_t ii, size_t)->numt {return (ii == i) ? 1 : 0;});
}
template<class numt>
const MatrixByFormula<numt> Diagonal(const Chain<numt> &V)
{
    return MatrixByFormula<numt>(V.size(), V.size(), [V](size_t i, size_t j)->numt {return (i == j) ? V[i] : 0;});
}
template<class numt>
const MatrixByFormula<numt> Permutation(const Chain<size_t> &V)
{
    for (const size_t i : V)if (i >= V.size())
            throw Exception<MatrixByFormula<numt>>("invalid permutation matrix");
    return MatrixByFormula<numt>(V.size(), V.size(), [V](size_t i, size_t j)->numt {return j == V[i] ? 1 : 0;});
}


template<class numt>
const MatrixByFormula<numt> operator*(const Matrix<numt> &source, const numt &k)
{
    return MatrixByFormula<numt>(source.height(), source.width(),
                                 [&source, &k](size_t i, size_t j)->numt {return source(i, j) * k;}
                                );
}
template<class numt>
const MatrixByFormula<numt> operator/(const Matrix<numt> &source, const numt &k)
{
    return MatrixByFormula<numt>(source.height(), source.width(),
                                 [&source, &k](size_t i, size_t j)->numt {return source(i, j) / k;}
                                );
}
template<class numt>
const MatrixByFormula<numt> Transponate(const Matrix<numt> &source)
{
    return MatrixByFormula<numt>(source.width(), source.height(),
                                 [&source](size_t i, size_t j)->numt {return source(j, i);}
                                );
}
template<class numt>
const MatrixByFormula<numt> Minor(const Matrix<numt> &source, const size_t i, const size_t j)
{
    if ((source.height() < 2) || (source.width() < 2) || (i >= source.height()) || (j >= source.width()))
        throw Exception<MatrixByFormula<numt>>("Invalid minor");
    return MatrixByFormula<numt>(source.height() - 1, source.width() - 1,
    [&source, i, j](size_t ii, size_t jj)->numt {
        auto new_ii = ii, new_jj = jj;
        if (ii >= i)new_ii++;
        if (jj >= j)new_jj++;
        return source(new_ii, new_jj);
    }
                                );
}
template<class numt>
const numt Suplement(const Matrix<numt> &source, const size_t i, const size_t j)
{
    switch ((i + j) % 2) {
    case 0:
        return Determinant(Minor(source, i, j));
    case 1:
        return -Determinant(Minor(source, i, j));
    }
    throw;
}
template<class numt>
const MatrixByFormula<numt> Suplements(const Matrix<numt> &source)
{
    return MatrixByFormula<numt>(
               source.height(), source.width(),
    [&source](size_t i, size_t j)->numt {
        return Suplement(source, i, j);
    }
           );
}
template<class numt>
const numt Determinant(const Matrix<numt> &source)
{
    if ((source.height() != source.width()) || (source.height() == 0))
        throw Exception<MatrixByFormula<numt>>("Cannot calculate the determinant");
    if (source.height() == 1)return source(0, 0);
    numt result = 0;
    for (size_t i = 0; i < source.width(); i++) {
        auto a = source(0, i);
        if (a != 0)
            result += a * Suplement(source, 0, i);
    }
    return result;
}
template<class numt>
const MatrixData<numt> CalcInverseMatrix(const Matrix<numt> &source)
{
    numt D = Determinant(source);
    if (D == 0)throw Exception<MatrixByFormula<numt>>("Cannot calculate inverse matrix");
    if (source.height() == 1)
        return MatrixData<numt>(1) / D;
    return Transponate(Suplements(source)) / D;
}
template<class numt>
const MatrixByFormula<numt> Multiply(const Matrix<numt> &A, const Matrix<numt> &B)
{
    if (A.width() != B.height())
        throw Exception<Matrix<numt>>("Matrix Multiplication: size mismatch");
    size_t cycle_length = A.width();
    return MatrixByFormula<numt>(A.height(), B.width(), [&A, &B, cycle_length](size_t i, size_t j)->numt {
        numt result = 0;
        for (size_t k = 0; k < cycle_length; k++)
            result += A(i, k) * B(k, j);
        return result;
    });
}
template<class numt>
struct lup {
    MatrixData<numt> L, U;
    MatrixByFormula<numt>P;
};
template<class numt>
const lup<numt> LUP(const Matrix<numt> &A_)
{
    if (A_.width() != A_.height())
        throw Exception<Matrix<numt>>("LUP decomposition is possible only for squared matrix");
    MatrixData<numt> A = A_;
    const size_t n = A.width();
    Chain<size_t>pi;
    for (size_t i = 0; i < n; i++)pi.push_back(i);
    for (size_t k = 0; k < n; k++) {
        size_t new_k = k;
        numt p = 0;
        for (size_t i = k; i < n; i++) {
            const auto a = pow(A(i, k), 2);
            if (a > p) {
                p = a;
                new_k = i;
            }
        }
        if (0 == p)throw Exception<Matrix<numt>>("LUP: the matrix is singular");
        if (new_k != k) {
            const auto ii = pi[k];
            pi[k] = pi[new_k];
            pi[new_k] = ii;
            for (size_t i = 0; i < n; i++) {
                const auto a = A(k, i);
                A.var_element(k, i) = A(new_k, i);
                A.var_element(new_k, i) = a;
            }
        }
        for (size_t i = k + 1; i < n; i++) {
            A.var_element(i, k) /= A(k, k);
            for (size_t j = k + 1; j < n; j++) {
                A.var_element(i, j) -= A(i, k) * A(k, j);
            }
        }
    }
    MatrixData<numt> L(Zeros<numt>(n, n)), U(Zeros<numt>(n, n));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i > j)L.var_element(i, j) = A(i, j);
            else U.var_element(i, j) = A(i, j);
        }
        L.var_element(i, i) = 1;
    }
    return {.L = L, .U = U, .P = Permutation<numt>(pi)};
}
template<class numt>
const MatrixData<numt> Solve(const Matrix<numt> &A, const Matrix<numt> &B)
{
    if (A.width() != A.height())
        throw Exception<Matrix<numt>>("Solve: parameter 1 must be squared matrix");
    if (B.width() != 1)
        throw Exception<Matrix<numt>>("Solve: parameter 2 must be a vector");
    if (B.height() != A.height())
        throw Exception<Matrix<numt>>("Solve: matrix and vector size mismatch");
    const auto n = A.width();
    const auto D = LUP(A);
    const MatrixData<numt> b = Multiply(D.P, B);
    if (b.width() > 1)
        throw Exception<Matrix<numt>>("");
    Chain<numt> y;
    for (size_t i = 0; i < n; i++) {
        numt v = b(i, 0);
        for (size_t j = 0; j < i; j++)
            v -= y[j] * D.L(i, j);
        y.push_back(v);
    }
    MatrixData<numt> X = Zeros<numt>(n, 1);
    for (size_t ii = n; ii > 0; ii--) {
        const size_t i = ii - 1;
        numt x = y[i];
        for (size_t j = ii; j < n; j++)
            x -= D.U(i, j) * X(j, 0);
        X.var_element(i, 0) = x / D.U(i, i);
    }
    return X;
}
}
#endif
