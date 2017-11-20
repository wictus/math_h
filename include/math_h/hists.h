// this file is distributed under
// LGPLv3 license
#ifndef ____HISTS_H_____
#	define ____HISTS_H_____
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include "tabledata.h"
#include "sigma.h"
#include "error.h"
namespace MathTemplates
{
template<class numt>
const Chain<value<numt>> BinsByStep(const numt from, const numt step, const numt to)
{
    if (0 >= step)throw Exception<Chain<value<numt>>>("wrong bin width");
    if (to <= from)throw Exception<Chain<value<numt>>>("wrong range");
    numt delta = step / numt(2);
    Chain<value<numt>> res;
    for (numt x = from + delta; x < to; x += step)
        res.push_back(value<numt>(x, delta));
    return res;
}
template<class numt>
const Chain<value<numt>> BinsByCount(const size_t count, const numt from, const numt to)
{
    if (0 == count)throw Exception<Chain<value<numt>>>("wrong bins count");
    return BinsByStep(from, (to - from) / numt(count), to);
}

template<class numtX = double, class numtY = numtX>
class hist: public SortedPoints<value<numtX>, value<numtY>>
{
public:
    typedef point<value<numtX>, value<numtY>> Point;
    hist() {}
    hist(const SortedChain<value<numtX>> &data)
    {
        for (const auto &v : data)this->append_item_from_sorted(Point(v,numtY(0)));
    }
    hist(const SortedChain<point<numtX, numtY>> &data)
    {
        for (const auto &P : data)this->append_item_from_sorted(Point(P));
    }
    hist(const SortedChain<Point> &data)
    {
        for (const auto &P : data)this->append_item_from_sorted(P);
    }
    hist(const typename SortedPoints<value<numtX>, value<numtY>>::Func f, const SortedChain<value<numtX>> &chain)
        : SortedPoints<value<numtX>, value<numtY>>(f, chain) {}
    hist(const SortedPoints<value<numtX>, value<numtY>> &source)
        : SortedPoints<value<numtX>, value<numtY>>(source) {}
    hist(const Points<numtX, numtY> &data)
    {
        for (const auto &P : data)this->operator<<(Point(P));
    }
    hist(const Chain<Point> &data)
    {
        for (const auto &P : data)this->operator<<(P);
    }
    hist Clone()const
    {
        return hist(*this);
    }
    virtual ~hist() {}
    hist &operator=(const SortedPoints<value<numtX>, value<numtY>> &points)
    {
        SortedPoints<value<numtX>, value<numtY>>::operator=(points);
        return *this;
    }
    hist CloneEmptyBins()const
    {
        Chain<Point> initer;
        for (const Point &P : *this)initer.push_back(P);
        return hist(initer);
    }
    value<numtY> TotalSum()const
    {
        value<numtY> res = 0;
        for (const Point &P : *this)res += P.Y();
        return res;
    }
    const SortedPoints<numtX, numtY> toLine()const
    {
        SortedPoints<numtX, numtY> res;
        for (const Point &P : *this)
            res << point<numtX, numtY>(P.X().val(), P.Y().val());
        return res;
    }
    const hist Scale(const size_t sc_x)const
    {
        SortedChain<value<numtX>> new_x, sorted_x;
        for (const auto &item : *this)
            sorted_x << item.X();
        for (size_t i = sc_x - 1, n = sorted_x.size(); i < n; i += sc_x) {
            auto min = sorted_x[i + 1 - sc_x].min();
            auto max = sorted_x[i].max();
            new_x << value<numtX>((max + min) / numtX(2), (max - min) / numtX(2));
        }
        hist res(new_x);
        for (size_t i = 0; i < new_x.size(); i++) {
            numtY v = 0;
            for (size_t ii = 0; ii < sc_x; ii++)
                v += this->operator[](i * sc_x + ii).Y().val();
            res.Bin(i) = value<numtY>::std_error(v);
        }
        return res;
    }
    hist &imbibe(const hist &second)
    {
        for (int i = 0, n = this->size(); i < n; i++) {
            if (this->operator[](i).X() == second[i].X()) {
                this->Bin(i) = value<numtY>::std_error(this->operator[](i).Y().val() + second[i].Y().val());
            } else
                throw Exception<hist>("Cannot imbibe histogram. bins differ");
        }
        return *this;
    }
};
template<class numtX = double, class numtY = numtX, class numtZ = numtY>
class hist2d: public BiSortedPoints<value<numtX>, value<numtY>, value<numtZ>>
{
public:
    hist2d(const SortedChain<value<numtX>> &X, const SortedChain<value<numtY>> &Y)
        : BiSortedPoints<value<numtX>, value<numtY>, value<numtZ>>(X, Y) {}
    hist2d(): hist2d({}, {}) {}
    hist2d(const BiSortedPoints<value<numtX>, value<numtY>, value<numtZ>> &source)
        : BiSortedPoints<value<numtX>, value<numtY>, value<numtZ>>(source) {}
    hist2d Clone()const
    {
        return hist2d(*this);
    }

    virtual ~hist2d() {}

    const hist2d Scale(const size_t sc_x, const size_t sc_y)const
    {
        SortedChain<value<numtX>> new_x;
        for (size_t i = sc_x - 1, n = this->X().size(); i < n; i += sc_x) {
            auto min = this->X()[i + 1 - sc_x].min();
            auto max = this->X()[i].max();
            new_x << (value<numtX>((max + min) / numtX(2), (max - min) / numtX(2)));
        }
        SortedChain<value<numtY>> new_y;
        for (size_t i = sc_y - 1, n = this->Y().size(); i < n; i += sc_y) {
            auto min = this->Y()[i + 1 - sc_y].min();
            auto max = this->Y()[i].max();
            new_y << (value<numtY>((max + min) / numtY(2), (max - min) / numtY(2)));
        }
        hist2d res(new_x, new_y);
        for (size_t i = 0; i < new_x.size(); i++)for (size_t j = 0; j < new_y.size(); j++) {
                numtZ v = 0;
                for (size_t ii = 0; ii < sc_x; ii++) {
                    for (size_t jj = 0; jj < sc_y; jj++)
                        v += this->operator[](i * sc_x + ii)[j * sc_y + jj].val();
                }
                res.Bin(i, j) = value<numtZ>::std_error(v);
            }
        return res;
    }
    hist2d &imbibe(const hist2d &second)
    {
        if ((this->X().size() != second.X().size()) || (this->Y().size() != second.Y().size()))
            throw Exception<hist2d>("cannot imbibe second histogram: bins differ");
        for (int i = 0, n = this->size(); i < n; i++)for (int j = 0, m = this->operator[](i).size(); j < m; j++)
                this->Bin(i, j) = value<numtY>::std_error(this->operator[](i)[j].val() + second[i][j].val());
        return *this;
    }
};

template<class numtX = double, class numtY = numtX>
class Distribution1D: public hist<numtX, numtY>
{
private:
    unsigned long long counter;
public:
    Distribution1D(const Chain<value<numtX>> &data): hist<numtX, numtY>(data)
    {
        this->FillWithValues(value<numtY>::std_error(0));
        counter = 0;
    }
    Distribution1D(const SortedChain<value<numtX>> &data): hist<numtX, numtY>(data)
    {
        this->FillWithValues(value<numtY>::std_error(0));
        counter = 0;
    }
    Distribution1D &Fill(const numtX &v)
    {
        counter++;
        for (size_t i = 0, n = this->size(); i < n; i++) {
            if (this->Bin(i).X().Contains(v))
                this->Bin(i) = value<numtY>::std_error(this->Bin(i).Y().val() + numtY(1));
        }
        return *this;
    }
    const unsigned long long Entries()const
    {
        return counter;
    }
};
template<class numtX = double, class numtY = numtX, class numtZ = numtY>
class Distribution2D: public hist2d<numtX, numtY, numtZ>
{
private:
    unsigned long long counter;
    void init()
    {
        counter = 0;
        this->FullCycleVar([](const value<numtX> &, const value<numtY> &, value<numtZ> &z) {
            z = value<numtZ>::std_error(0);
        });
    }
public:
    Distribution2D(const Chain<value<numtX>> &X, const Chain<value<numtY>> &Y)
        : hist2d<numtX, numtY, numtZ>(X, Y)
    {
        init();
    }
    Distribution2D(const SortedChain<value<numtX>> &X, const SortedChain<value<numtY>> &Y)
        : hist2d<numtX, numtY, numtZ>(X, Y)
    {
        init();
    }
    Distribution2D &Fill(const numtX &x, const numtY &y)
    {
        counter++;
        for (size_t i = 0, I = this->X().size(); i < I; i++)if (this->X()[i].Contains(x)) {
                for (size_t j = 0, J = this->Y().size(); j < J; j++)if (this->Y()[j].Contains(y))
                        this->Bin(i, j) = value<numtZ>::std_error(this->operator[](i)[j].val() + numtZ(1));
            }
        return *this;
    }
    const unsigned long long Entries()const
    {
        return counter;
    }
};
};
#endif
