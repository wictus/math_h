// this file is distributed under
// LGPLv3 license
#ifndef ____table_data_H_____
#	define ____table_data_H_____
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <vector>
#include <functional>
#include "error.h"
namespace MathTemplates
{
template<class comparable = double>
using Chain=std::vector<comparable>;
namespace details
{
template<class comparable, class indexer = Chain<comparable>>
size_t  WhereToInsert(const size_t from, const size_t to, const indexer &X, const comparable &x)
{
    if (from > to) return from;
    size_t beg = from, end = to;
    if (x > X[end]) return end + 1;
    if (x < X[beg]) return beg;
    while (1 < (end - beg)) {
        size_t mid = (beg + end) / 2;
        if (x < X[mid]) end = mid;
        else if (x > X[mid]) beg = mid;
        else return mid;
    }
    return end;
}
template<class comparable, class indexer, class Size, class Insert>
void InsertSorted(const comparable &x, indexer &X, const Size size, const Insert insert)
{
    if (size() == 0) insert(0, x);
    else insert(WhereToInsert(0, size() - 1, X, x), x);
}
template<class comparable, class indexer, class Size, class Insert>
void InsertSorted(const comparable &&x, indexer &X, const Size size, const Insert insert)
{
    InsertSorted(x, X, size, insert);
}
}
#define std_size(vector) [&vector](){return (vector).size();}
#define std_insert(vector,type) [&vector](int pos,type x){(vector).insert((vector).begin()+pos,x);}
#define field_size(vector)  [this](){return (vector).size();}
#define field_insert(vector,type)  [this](int pos,type x){(vector).insert((vector).begin()+pos,x);}
template<class comparable = double>
class SortedChain
{
private:
    Chain<comparable> data;
public:
    const Chain<comparable>&operator()()const{return data;}
    SortedChain() {}
    SortedChain(const SortedChain &points)
    {
        for (const auto &p : points.data)
            data.push_back(p);
    }
    const SortedChain IndexRange(const size_t from, const size_t count)const
    {
        SortedChain result;
        size_t final = from + count;
        for (size_t i = from; (i < final) && (i < data.size()); i++)
            result.data.push_back(data[i]);
        return result;
    }
    const SortedChain SelectByCondition(const std::function<bool(const comparable &)> condition)const
    {
        SortedChain result;
        for (size_t i = 0; i < data.size(); i++)
            if (condition(data[i]))
                result.data.push_back(data[i]);
        return result;
    }
    SortedChain &operator<<(const comparable &p)
    {
        details::InsertSorted(p, data, field_size(data), field_insert(data, comparable));
        return *this;
    }
    SortedChain(const std::vector<comparable> &points)
    {
        for (const auto &p : points)
            operator<<(p);
    }
    virtual ~SortedChain() {}
    SortedChain &operator=(const SortedChain &points)
    {
        data.clear();
        for (const auto &p : points.data)
            data.push_back(p);
        return *this;
    }
    void clear()
    {
        data.clear();
    }

    const size_t size()const
    {
        return data.size();
    }
    const comparable &operator[](const size_t i)const
    {
        if (size() <= i)
            throw Exception<SortedChain>("Range check error");
        return data[i];
    }
    const comparable &left()const
    {
        if (size() < 1)
            throw Exception<SortedChain>("Attempt to obtain empty properties.");
        return data[0];
    }
    const comparable &right()const
    {
        if (size() < 1)
            throw Exception<SortedChain>("Attempt to obtain empty properties.");
        return data[size() - 1];
    }
    typedef typename std::vector<comparable>::const_iterator const_iterator;
    const_iterator begin()const
    {
        return data.begin();
    }
    const_iterator cbegin()const
    {
        return data.cbegin();
    }
    const_iterator end() const
    {
        return data.end();
    }
    const_iterator cend() const
    {
        return data.cend();
    }
protected:
    comparable &accessBin(const size_t i)
    {
        if (data.size() <= i)
            throw Exception<SortedChain>("range check error");
        return data[i];
    }
    void append_item_from_sorted(const comparable &c)
    {
        data.push_back(c);
    }
};
template<class numX>
const Chain<numX> ChainWithStep(const numX &from, const numX &step, const numX &to)
{
    if (from >= to)throw Exception<Chain<numX>>("wrong binning ranges");
    if (step <= 0)throw Exception<Chain<numX>>("wrong binning step");
    Chain<numX> res;
    for (numX x = from; x <= to; x += step)res.push_back(x);
    return res;
}
template<class numX>
const Chain<numX> ChainWithCount(const size_t cont, const numX &from, const numX &to)
{
    if (from >= to)throw Exception<Chain<numX>>("wrong binning ranges");
    if (0 == cont)throw Exception<Chain<numX>>("wrong bins count");
    numX step = (to - from) / numX(cont);
    Chain<numX> res;
    for (numX x = from; x <= to; x += step)res.push_back(x);
    return res;
}

template<class numtX = double, class numtY = numtX>
class point
{
private:
    numtX x;
    numtY y;
public:
    const numtX &X()const
    {
        return x;
    }
    const numtY &Y()const
    {
        return y;
    }
    point(const numtX &pos): x(pos), y(numtY(0)) {}
    point(const numtX &pos, const numtY &val): x(pos), y(val) {}
    template<class numtX2, class numtY2>
    point(const point<numtX2,numtY2> &source): x(source.X()), y(source.Y()) {}
    template<class numt = numtY>
    point(const std::initializer_list<numt> &source)
    {
        if (source.size() == 0)
            throw Exception<point>("wrong initialization of point from emply list");
        if (source.size() > 2)
            throw Exception<point>("wrong initialization of value from list with more than two numbers");
        Chain<numt> v;
        for (const numt &x : source)v.push_back(x);
        x = numtX(v[0]);
        y = numtY(v[1]);
    }
    template<class numtX2 = numtX, class numtY2 = numtY>
    point &operator=(const point<numtX2,numtY2> &source)
    {
        x = source.X();
        y = source.Y();
	return *this;
    }
    point &operator=(const numtY &s)
    {
        y = s;
	return *this;
    }
    bool operator<(const point &b)const
    {
        return x < b.x;
    }
    bool operator>(const point &b)const
    {
        return x > b.x;
    }
    const point operator+(const numtY &val)const
    {
        return {X(), Y() + val};
    }
    const point operator-(const numtY &val)const
    {
        return {X(), Y() - val};
    }
    const point operator*(const numtY &val)const
    {
        return {X(), Y() * val};
    }
    const point operator/(const numtY &val)const
    {
        return {X(), Y() / val};
    }
    const point operator+(const point &other)const
    {
        if (other.X() == X())
            return operator+(other.Y());
        else
            throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
    }
    const point operator-(const point &other)const
    {
        if (other.X() == X())
            return operator-(other.Y());
        else
            throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
    }
    const point operator*(const point &other)const
    {
        if (other.X() == X())
            return operator*(other.Y());
        else
            throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
    }
    const point operator/(const point &other)const
    {
        if (other.X() == X())
            return operator/(other.Y());
        else
            throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
    }
};
template<class numtX = double, class numtY = numtX>
inline const point<numtX,numtY>make_point(const numtX&x,const numtY&y){return point<numtX,numtY>(x,y);}
template<class numtX = double, class numtY = numtX, class numtZ = numtY>
class point3d
{
private:
    numtX x;
    numtY y;
    numtZ z;
public:
    point3d(const numtX &_x, const numtY &_y, const numtZ &_z)
        : x(_x), y(_y), z(_z) {}
    virtual ~point3d() {}
    const numtX &X()const
    {
        return x;
    }
    const numtY &Y()const
    {
        return y;
    }
    const numtZ &Z()const
    {
        return z;
    }
    template<class numtX2, class numtY2,class numtZ2>
    point3d(const point3d<numtX2,numtY2,numtZ2>&source)
	: x(source.X()), y(source.Y()), z(source.Z()){}
};
template<class numtX = double, class numtY = numtX,class numtZ=numtY>
inline const point3d<numtX,numtY,numtZ>make_point(const numtX&x,const numtY&y,const numtZ&z){
    return point3d<numtX,numtY,numtZ>(x,y,z);
}


template<class numX = double, class numY = numX>
using Points=Chain<point<numX, numY>>;
template<class numX = double, class numY = numX>
class SortedPoints: public SortedChain<point<numX, numY>>
{
public:
    typedef std::function<numY(const numX &)> Func;
    SortedPoints() {}
    SortedPoints(const Points<numX, numY> &chain)
    {
        for (const auto &x : chain)
            SortedChain<point<numX, numY>>::operator<<(x);
    }
    SortedPoints(const SortedChain<point<numX, numY>> &chain)
    {
        for (const auto &x : chain)
            SortedChain<point<numX, numY>>::append_item_from_sorted(x);
    }
    SortedPoints(const Func f, const SortedChain<numX> &chain)
    {
        for (numX x : chain)
            SortedChain<point<numX, numY>>::append_item_from_sorted(point<numX, numY>(x, f(x)));
    }
    virtual ~SortedPoints() {}
    SortedPoints Clone()const
    {
        return SortedPoints(*this);
    }
protected:
    point<numX, numY> &Bin(const size_t i)
    {
        return SortedChain<point<numX, numY>>::accessBin(i);
    }
public:
    const  Points<numY, numX> Transponate()const
    {
        Points<numY, numX> res;
        for (const auto &p : *this)
            res.push_back(point<numX, numY>(p.Y(), p.X()));
        return res;
    }
    const  SortedChain<point<numY, numX>> TransponateAndSort()const
    {
        SortedChain<point<numY, numX>> res;
        for (const auto &p : *this)
            res << point<numX, numY>(p.Y(), p.X());
        return res;
    }
    const SortedPoints XRange(const numX &from, const numX &to)const
    {
        return SortedChain<point<numX, numY>>::SelectByCondition([from, to](const point<numX, numY> &P) {
            return (P.X() >= from) && (P.X() <= to);
        });
    }
    const SortedPoints YRange(const numY &from, const numY &to)const
    {
        return SortedChain<point<numX, numY>>::SelectByCondition([from, to](const point<numX, numY> &P) {
            return (P.Y() >= from) && (P.Y() <= to);
        });
    }
    const SortedPoints XExclude(const numX &from, const numX &to)const
    {
        return SortedChain<point<numX, numY>>::SelectByCondition([from, to](const point<numX, numY> &P) {
            return (P.X() < from) || (P.X() > to);
        });
    }
    const SortedPoints YExclude(const numY &from, const numY &to)const
    {
        return SortedChain<point<numX, numY>>::SelectByCondition([from, to](const point<numX, numY> &P) {
            return (P.Y() < from) || (P.Y() > to);
        });
    }
    SortedPoints &FillWithValues(const numY &v)
    {
        for (size_t i = 0, n = this->size(); i < n; i++)
            Bin(i) = v;
        return *this;
    }
    SortedPoints &Transform(const std::function<numY(const numX &, const numY &)> &F)
    {
        for (size_t i = 0, n = this->size(); i < n; i++)
            Bin(i) = F(Bin(i).X(), Bin(i).Y());
        return *this;
    }
    SortedPoints &operator+=(const SortedPoints &second)
    {
	if(second.size()!=this->size())
	    throw Exception<SortedPoints>("Cannot perform arithmetic operation for two sets of points with different size");
        for (size_t i = 0, n = this->size(); i < n; i++) {
	    Bin(i)=Bin(i)+second[i];
        }
        return *this;
    }
    SortedPoints &operator+=(const Func f)
    {
        for (size_t i = 0, n = this->size(); i < n; i++)
            Bin(i)=Bin(i)+f(Bin(i).X());
        return *this;
    }
    SortedPoints &operator+=(const numY &c)
    {
        for (size_t i = 0, n = this->size(); i < n; i++)
            Bin(i)=Bin(i)+ c;
        return *this;
    }

    SortedPoints &operator-=(const SortedPoints &second)
    {
	if(second.size()!=this->size())
	    throw Exception<SortedPoints>("Cannot perform arithmetic operation for two sets of points with different size");
        for (size_t i = 0, n = this->size(); i < n; i++) {
	    Bin(i)=Bin(i)-second[i];
        }
        return *this;
    }
    SortedPoints &operator-=(const Func f)
    {
        for (size_t i = 0, n = this->size(); i < n; i++)
            Bin(i)=Bin(i)-f(Bin(i).X());
        return *this;
    }
    SortedPoints &operator-=(const numY &c)
    {
        for (size_t i = 0, n = this->size(); i < n; i++)
            Bin(i)=Bin(i)-c;
        return *this;
    }

    SortedPoints &operator*=(const SortedPoints &second)
    {
	if(second.size()!=this->size())
	    throw Exception<SortedPoints>("Cannot perform arithmetic operation for two sets of points with different size");
        for (size_t i = 0, n = this->size(); i < n; i++) {
	    Bin(i)=Bin(i)*second[i];
        }
        return *this;
    }
    SortedPoints &operator*=(const Func f)
    {
        for (size_t i = 0, n = this->size(); i < n; i++)
            Bin(i)=Bin(i)*f(Bin(i).X());
        return *this;
    }
    SortedPoints &operator*=(const numY &c)
    {
        for (size_t i = 0, n = this->size(); i < n; i++)
            Bin(i)=Bin(i)*c;
        return *this;
    }

    SortedPoints &operator/=(const SortedPoints &second)
    {
	if(second.size()!=this->size())
	    throw Exception<SortedPoints>("Cannot perform arithmetic operation for two sets of points with different size");
        for (size_t i = 0, n = this->size(); i < n; i++) {
	    Bin(i)=Bin(i)/second[i];
        }
        return *this;
    }
    SortedPoints &operator/=(const Func f)
    {
        for (size_t i = 0, n = this->size(); i < n; i++)
            Bin(i)=Bin(i)/f(Bin(i).X());
        return *this;
    }
    SortedPoints &operator/=(const numY &c)
    {
        for (size_t i = 0, n = this->size(); i < n; i++)
            Bin(i)=Bin(i)/c;
        return *this;
    }

    const SortedPoints operator+(const SortedPoints &other)const
    {
        return SortedPoints(*this) += other;
    }
    const SortedPoints operator-(const SortedPoints &other)const
    {
        return SortedPoints(*this) -= other;
    }
    const SortedPoints operator*(const SortedPoints &other)const
    {
        return SortedPoints(*this) *= other;
    }
    const SortedPoints operator/(const SortedPoints &other)const
    {
        return SortedPoints(*this) /= other;
    }

    const SortedPoints operator+(const numY &other)const
    {
        return SortedPoints(*this) += other;
    }
    const SortedPoints operator-(const numY &other)const
    {
        return SortedPoints(*this) -= other;
    }
    const SortedPoints operator*(const numY &other)const
    {
        return SortedPoints(*this) *= other;
    }
    const SortedPoints operator/(const numY &other)const
    {
        return SortedPoints(*this) /= other;
    }

    const SortedPoints operator+(const Func other)const
    {
        return SortedPoints(*this) += other;
    }
    const SortedPoints operator-(const Func other)const
    {
        return SortedPoints(*this) -= other;
    }
    const SortedPoints operator*(const Func other)const
    {
        return SortedPoints(*this) *= other;
    }
    const SortedPoints operator/(const Func other)const
    {
        return SortedPoints(*this) /= other;
    }

};
template<class numtX = double, class numtY = numtX, class numtZ = numtY>
class BiSortedPoints
{
private:
    SortedChain<numtX> m_x_axis;
    SortedChain<numtY> m_y_axis;
    Chain<Chain<numtZ>> m_data;
    void init()
    {
        m_data.clear();
        for (size_t i = 0, I = m_x_axis.size(); i < I; i++) {
            m_data.push_back(Chain<numtZ>());
            for (size_t j = 0, J = m_y_axis.size(); j < J; j++)
                m_data[m_data.size() - 1].push_back(numtZ(0));
        }
    }
public:
    const SortedChain<numtX> &X()const
    {
        return m_x_axis;
    }
    const SortedChain<numtY> &Y()const
    {
        return m_y_axis;
    }
    BiSortedPoints(const Chain<numtX> &X, const Chain<numtY> &Y)
        : m_x_axis(X), m_y_axis(Y)
    {
        init();
    }
    BiSortedPoints(const SortedChain<numtX> &X, const SortedChain<numtY> &Y)
        : m_x_axis(X), m_y_axis(Y)
    {
        init();
    }
    BiSortedPoints(): BiSortedPoints({}, {}) {}
    BiSortedPoints(const BiSortedPoints &source): m_x_axis(source.X()), m_y_axis(source.Y())
    {
        for (size_t i = 0, I = source.m_data.size(); i < I; i++) {
            m_data.push_back(Chain<numtZ>());
            for (const auto &item : source.m_data[i])
                m_data[i].push_back(item);
        }
    }
    BiSortedPoints Clone()const
    {
        return BiSortedPoints(*this);
    }
    virtual ~BiSortedPoints() {}
    typedef typename Chain<Chain<numtZ>>::const_iterator const_iterator;
    const_iterator begin()const
    {
        return m_data.cbegin();
    }
    const_iterator cbegin()const
    {
        return m_data.cbegin();
    }
    const_iterator end() const
    {
        return m_data.cend();
    }
    const_iterator cend() const
    {
        return m_data.cend();
    }
    size_t size()const
    {
        return m_data.size();
    }
    const Chain<numtZ> &operator[](const size_t i)const
    {
        if (size() <= i)throw Exception<BiSortedPoints>("X-range check error");
        return m_data[i];
    }
    numtZ &Bin(const size_t i, const size_t j)
    {
        if (size() <= i)throw Exception<BiSortedPoints>("X-range check error " + std::to_string(i) + ":" + std::to_string(j));
        if (m_data[i].size() <= j)throw Exception<BiSortedPoints>("Y-range check error " + std::to_string(i) + ":" + std::to_string(j));
        return m_data[i][j];
    }
    const SortedPoints<numtX, numtZ>CutX(const size_t j)const
    {
        if (m_y_axis.size() <= j)throw Exception<BiSortedPoints>("range check error");
        SortedPoints<numtX, numtZ> res;
        for (size_t i = 0, I = m_x_axis.size(); i < I; i++) {
            res << point<numtX, numtZ>(m_x_axis[i], m_data[i][j]);
        }
        return res;
    }
    const SortedPoints<numtY, numtZ>CutY(const size_t i)const
    {
        if (m_x_axis.size() <= i)throw Exception<BiSortedPoints>("range check error");
        SortedPoints<numtY, numtZ> res;
        for (size_t j = 0, J = m_y_axis.size(); j < J; j++) {
            res << point<numtY, numtZ>(m_y_axis[j], m_data[i][j]);
        }
        return res;
    }
    const point3d<numtX, numtY, numtZ> operator()(const size_t i, const size_t j)const
    {
        if (size() <= i)throw Exception<BiSortedPoints>("range check error");
        if (m_y_axis.size() <= j)throw Exception<BiSortedPoints>("range check error");
        return point3d<numtX, numtY, numtZ>(m_x_axis[i], m_y_axis[j], m_data[i][j]);
    }
    void FullCycle(const std::function<void(const point3d<numtX, numtY, numtZ>&)>f)const
    {
        for (size_t i = 0, I = m_x_axis.size(); i < I; i++)
            for (size_t j = 0, J = m_y_axis.size(); j < J; j++) {
                point3d<numtX, numtY, numtZ> P(m_x_axis[i], m_y_axis[j], m_data[i][j]);
                f(P);
            }
    }
    void FullCycle(const std::function<void(const numtX &, const numtY &, const numtZ &)>f)const
    {
        for (size_t i = 0, I = m_x_axis.size(); i < I; i++)
            for (size_t j = 0, J = m_y_axis.size(); j < J; j++)
                f(m_x_axis[i], m_y_axis[j], m_data[i][j]);
    }
    BiSortedPoints &FullCycleVar(const std::function<void(const numtX &, const numtY &, numtZ &)>f)
    {
        for (size_t i = 0, I = m_x_axis.size(); i < I; i++) {
            for (size_t j = 0, J = m_y_axis.size(); j < J; j++)
                f(m_x_axis[i], m_y_axis[j], m_data[i][j]);
        }
        return *this;
    }
    BiSortedPoints &FullCycleVar(const std::function<void(numtZ &)>f)
    {
        for (size_t i = 0, I = m_x_axis.size(); i < I; i++) {
            for (size_t j = 0, J = m_y_axis.size(); j < J; j++)
                f(m_data[i][j]);
        }
        return *this;
    }
};
};
#endif
