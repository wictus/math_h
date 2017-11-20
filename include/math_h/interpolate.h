// this file is distributed under
// LGPLv3 license
#ifndef PPRNSANGJGXVGERD
#define PPRNSANGJGXVGERD
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include "functions.h"
#include "integrate.h"
namespace MathTemplates
{
template<class numX = double, class numY = numX>
class LinearInterpolation: public SortedPoints<numX, numY>, public IFunction<const numY, const numX &>
{
public:
    typedef point<numX, numY> Point;
    typedef typename SortedPoints<numX, numY>::Func Func;

public:
    LinearInterpolation() {}
    LinearInterpolation(const Points<numX, numY> &points)
        : SortedPoints<numX, numY>(points) {}
    LinearInterpolation(const SortedChain<point<numX, numY>> &points)
        : SortedPoints<numX, numY>(points) {}
    LinearInterpolation(const Func f, const SortedChain<numX> &chain)
        : SortedPoints<numX, numY>(f, chain) {}
    LinearInterpolation(const Func f, const Chain<numX> &chain)
        : SortedPoints<numX, numY>(f, chain) {}
    LinearInterpolation(const SortedPoints<numX, numY> &source)
        : SortedPoints<numX, numY>(source) {}
    virtual ~LinearInterpolation() {}
    virtual const numY operator()(const numX &x)const override
    {
        auto tbl = [this](size_t i)->const Point& {
            return this->operator[](i);
        };
        const size_t sz = this->size();
        if (x == tbl(0).X())
            return tbl(0).Y();
        if (x == tbl(sz - 1).X())
            return tbl(sz - 1).Y();
        const Point p(x, numY(0));
        const size_t i = details::WhereToInsert<Point>(0, sz - 1, *this, p);
        if ((i == 0) || (i >= sz))
            throw Exception<LinearInterpolation>("Attempt to interpolate outside the given region.");
        const auto Left = tbl(i - 1), Right = tbl(i);
        if ((x - Left.X()) < (Right.X() - x)) {
            numX k = (x - Left.X()) / (Right.X() - Left.X());
            return Left.Y() + (Right.Y() - Left.Y()) * numY(k);
        } else {
            numX k = (Right.X() - x) / (Right.X() - Left.X());
            return Right.Y() + (Left.Y() - Right.Y()) * numY(k);
        }
    }
};
template<class numX = double, class numY = numX>
class IntegratedLinearInterpolation: public SortedPoints<numX, numY>, public IFunction<const numY, const numX &>
{
public:
    typedef typename SortedPoints<numX, numY>::Func Func;
    SortedPoints<numX, numY> raw;
public:
    IntegratedLinearInterpolation(const LinearInterpolation<numX, numY> &source)
        : SortedPoints<numX, numY>(Int_Trapez_Table(source)), raw(source) {}
    IntegratedLinearInterpolation(const IntegratedLinearInterpolation<numX, numY> &source)
        : SortedPoints<numX, numY>(source), raw(source.raw) {}
    virtual ~IntegratedLinearInterpolation() {}
    virtual const numY operator()(const numX &x)const override
    {
        const size_t sz = raw.size();
        if (x == raw.left().X())
            return this->left().Y();
        if (x == raw.right().X())
            return this->right().Y();
        const size_t i = details::WhereToInsert<point<numX, numY>>(0, sz - 1, raw, point<numX, numY>(x, numY(0)));
        if ((i == 0) || (i >= sz))
            throw Exception<IntegratedLinearInterpolation>("Attempt to interpolate outside the given region.");
        const numY Lx = raw[i - 1].X(), Rx = raw[i].X();
        const numY dx = x - Lx, dX = Rx - Lx;
        const numY L0 = raw[i - 1].Y(); //,R0=raw[i].Y();
        const numY L1 = this->operator[](i - 1).Y(), R1 = this->operator[](i).Y();
        const numY k = (R1 - L1 - L0 * dX) / (dX * dX);
        return L1 + L0 * dx + k * dx * dx;
    }
};
template<class numX = double, class numY = numX>
class ReverseIntegratedLinearInterpolation: public SortedPoints<numY, numX>, public IFunction<const numX, const numY &>
{
public:
    typedef typename SortedPoints<numX, numY>::Func Func;
    SortedPoints<numX, numY> raw;
public:
    ReverseIntegratedLinearInterpolation(const LinearInterpolation<numX, numY> &source)
        : SortedPoints<numY, numX>(Int_Trapez_Table_PositiveStrict(source).TransponateAndSort()), raw(source) {}
    ReverseIntegratedLinearInterpolation(const ReverseIntegratedLinearInterpolation<numX, numY> &source)
        : SortedPoints<numY, numX>(source), raw(source.raw) {}
    virtual ~ReverseIntegratedLinearInterpolation() {}
    virtual const numX operator()(const numY &y)const override
    {
        const size_t sz = this->size();
        if (y == this->left().X())
            return this->left().Y();
        if (y == this->right().X())
            return this->right().Y();
        const size_t i = details::WhereToInsert<point<numY, numX>>(0, sz - 1, *this, point<numY, numX>(y, numX(0)));
        if ((i == 0) || (i >= sz))
            throw Exception<ReverseIntegratedLinearInterpolation>("Attempt to interpolate outside the given region.");
        const numY Lx = raw[i - 1].X(), Rx = raw[i].X(), dX = Rx - Lx;
        const numY L0 = raw[i - 1].Y(), R0 = raw[i].Y();
        const numY L1 = this->operator[](i - 1).X(), R1 = this->operator[](i).X();
        const numY dy = y - L1, dY = R1 - L1;
        if (L0 != R0) {
            const numY k = (R1 - L1 - L0 * dX) / (dX * dX);
            return Lx + (sqrt(L0 * L0 + numX(4) * k * dy) - L0) / (numX(2) * k);
        } else {
            return Lx + dy * dX / dY;
        }
    }
};

template<class numtX = double, class numtY = numtX, class numtZ = numtY>
class BiLinearInterpolation: public BiSortedPoints<numtX, numtY, numtZ>, public IFunction<numtZ, const numtX &, const numtY &>
{
public:
    BiLinearInterpolation(const Chain<numtX> &X, const Chain<numtY> &Y)
        : BiSortedPoints<numtX, numtY, numtZ>(X, Y) {}
    BiLinearInterpolation(const SortedChain<numtX> &X, const SortedChain<numtY> &Y)
        : BiSortedPoints<numtX, numtY, numtZ>(X, Y) {}
    BiLinearInterpolation()
        : BiLinearInterpolation({}, {}) {}
    BiLinearInterpolation(const BiSortedPoints<numtX, numtY, numtZ> &source)
        : BiSortedPoints<numtX, numtY, numtZ>(source) {}
    virtual ~BiLinearInterpolation() {}
    virtual numtZ operator()(const numtX &x, const numtY &y)const override
    {
        int szx = this->X().size();
        int i = details::WhereToInsert<numtX>(0, szx - 1, this->X(), x);
        if ((i <= 0) || (i >= szx))
            throw Exception<BiLinearInterpolation>("Attempt to interpolate outside the given region (x).");
        int szy = this->Y().size();
        int j = details::WhereToInsert<numtY>(0, szy - 1, this->Y(), y);
        if ((j <= 0) || (j >= szy))
            throw Exception<BiLinearInterpolation>("Attempt to interpolate outside the given region (y).");
        numtZ U = numtZ(this->X()[i] - this->X()[i - 1]);
        numtZ T = numtZ(this->Y()[j] - this->Y()[j - 1]);
        numtZ u = numtZ(x - this->X()[i - 1]) / U;
        numtZ t = numtZ(y - this->Y()[j - 1]) / T;
        numtZ _u = numtZ(this->X()[i] - x) / U;
        numtZ _t = numtZ(this->Y()[j] - y) / T;
        return this->operator[](i - 1)[j - 1] * _u * _t
               + this->operator[](i)[j - 1] * u * _t
               + this->operator[](i - 1)[j] * _u * t
               + this->operator[](i)[j] * u * t;
    }
};
};
#endif
