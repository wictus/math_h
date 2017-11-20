// this file is distributed under
// LGPLv3 license
#ifndef YJIOGPKSIIYJVKND
#define YJIOGPKSIIYJVKND
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <iostream>
#include <memory>
#include <math.h>
#include "error.h"
#include "tabledata.h"
namespace MathTemplates
{
template<typename numt = double>
class value
{
private:
    numt Value, Error;
    struct cache {
        numt epsilon, min, max;
        cache(const numt &e, const numt &b, const numt &a)
            : epsilon(e), min(b), max(a) {}
    };
    std::shared_ptr<cache> f_cache;
    inline void invalidate()
    {
        if (Error < 0)
            Error = INFINITY;
        f_cache = nullptr;
    }
    void calc()const
    {
        if (!f_cache) {
            const_cast<value &>(*this).f_cache = std::make_shared<cache>(
                    Error / Value, Value - Error, Value + Error
                                                 );
        }
    }
public:
    virtual ~value() {}
    inline const numt &val()const
    {
        return Value;
    }
    inline const numt &uncertainty()const
    {
        return Error;
    }
    value(): Value(numt(0)), Error(numt(0))
    {
        invalidate();
    }
    value(const numt v): Value(v), Error(numt(0))
    {
        invalidate();
    }
    value(const numt v, const numt err): Value(v), Error(err)
    {
        invalidate();
    }
    value(const std::function<numt(const numt&)>F,const value&source)
    :Value(F(source.val())),Error(sqrt((pow(F(source.min()) - Value, 2) + pow(F(source.max()) - Value, 2)) / numt(2)))
    {
	invalidate();
    }

    value(const std::initializer_list<numt> &source)
    {
        if (source.size() == 0)
            throw Exception<value>("wrong initialization of value from emply list");
        if (source.size() > 2)
            throw Exception<value>("wrong initialization of value from list with more than two numbers");
        Chain<numt> v;
        for (const numt &x : source)v.push_back(x);
        Value = v[0];
        if (v.size() == 1)Error = numt(0);
        else Error = v[1];
        invalidate();
    }
    static const value std_error(const numt &v)
    {
        if (v < 0)throw Exception<value>("Cannot calculate std error for negative value");
        auto res = value(v, sqrt(v));
        if (res.Error < numt(1))res.Error = numt(1);
        return res;
    }
    inline static const value interval(const numt &a, const numt &b)
    {
        if (b < a)throw Exception<value>("Bad interval");
        return value((b + a) / numt(2), (b - a) / numt(2));
    }
    template<class numt2>
    value(const value<numt2> &source): Value(source.val()), Error(source.uncertainty())
    {
        invalidate();
    }
    value &operator=(const value &source)
    {
        Value = source.Value;
        Error = source.Error;
        invalidate();
        return *this;
    }
    const value make_wider(const numt &scale)const
    {
        return value(Value, scale * Error);
    }
    const numt &epsilon()const
    {
        calc();
        return f_cache->epsilon;
    }
    const numt &min()const
    {
        calc();
        return f_cache->min;
    }
    const numt &max()const
    {
        calc();
        return f_cache->max;
    }
    //Physical comparing of magnitudes with uncertainties
    const bool Contains(const numt &x)const
    {
        return (x >= min()) && (x <= max());
    }
    const bool Contains(const value &x)const
    {
        return (x.max() >= min()) && (x.min() <= max());
    }
    const bool NotEqual(const numt &x)const
    {
        return (x < min()) || (x > max());
    }
    const bool NotEqual(const value &x)const
    {
        return (x.max() < min()) || (x.min() > max());
    }
    const bool Below(const numt &x)const
    {
        return max() < x;
    }
    const bool Below(const value &x)const
    {
        return max() < x.min();
    }
    const bool Above(const numt &x)const
    {
        return min() > x;
    }
    const bool Above(const value &x)const
    {
        return min() > x.max();
    }
    //chi-square-like numeric comparing of magnitudes
    const numt NumCompare(const numt &x)const
    {
        return pow((Value - x) / Error, 2);
    }
    const numt NumCompare(const value &x)const
    {
        return pow((Value - x.Value) / (Error + x.Error), 2);
    }
    //Inheriting number-like comparing
    inline const bool operator<(const value &other)const
    {
        return Value < other.Value;
    }
    inline const bool operator>(const value &other)const
    {
        return Value > other.Value;
    }
    inline const bool operator==(const value &other)const
    {
        return Value == other.Value;
    }
    inline const bool operator>=(const value &other)const
    {
        return Value >= other.Value;
    }
    inline const bool operator<=(const value &other)const
    {
        return Value <= other.Value;
    }
    //arithmetic actions
    value &operator+=(const value &other)
    {
        Error = sqrt(pow(Error, 2) + pow(other.Error, 2));
        Value += other.Value;
        invalidate();
        return *this;
    }
    value &operator-=(const value &other)
    {
        Error = sqrt(pow(Error, 2) + pow(other.Error, 2));
        Value -= other.Value;
        invalidate();
        return *this;
    }
    value &operator*=(const value &other)
    {
        Error = sqrt(pow(Error * other.Value, 2) + pow(other.Error * Value, 2));
        Value *= other.Value;
        invalidate();
        return *this;
    }
    value &operator/=(const value &other)
    {
        Error = sqrt(pow(Error / other.Value, 2)
                     + pow(other.Error * Value / pow(other.Value, 2), 2));
        Value /= other.Value;
        invalidate();
        return *this;
    }
    inline const value operator+(const value &other)const
    {
        return value(*this) += other;
    }
    inline const value operator+(const value&&other)const
    {
        return value(*this) += other;
    }
    inline const value operator-(const value &other)const
    {
        return value(*this) -= other;
    }
    inline const value operator-(const value&&other)const
    {
        return value(*this) -= other;
    }
    inline const value operator*(const value &other)const
    {
        return value(*this) *= other;
    }
    inline const value operator*(const value&&other)const
    {
        return value(*this) *= other;
    }
    inline const value operator/(const value &other)const
    {
        return value(*this) /= other;
    }
    inline const value operator/(const value&&other)const
    {
        return value(*this) /= other;
    }
};
template<class numt>
inline const value<numt> std_error(const numt &v)
{
    return value<numt>::std_error(v);
}
template<typename numt>
inline std::istream &operator>>(std::istream &str, value<numt> &P)
{
    numt v, u;
    str >> v >> u;
    P = {v, u};
    return str;
}
template<typename numt>
inline std::ostream &operator<<(std::ostream &str, const value<numt> &P)
{
    return str << P.val() << " " << P.uncertainty();
}

template<typename numt = double>
class StandardDeviation
{
private:
    Chain<numt> m_list;
    numt m_sum;
    std::shared_ptr<value<numt>>m_cache;
    numt m_scale;
public:
    typedef typename Chain<numt>::const_iterator const_iterator;
    const_iterator begin()const{return m_list.begin();}
    const_iterator end() const{return m_list.end();}
    const numt&operator[](const size_t index){
	if(index>=m_list.size())throw Exception<StandardDeviation>("range check error");
	return m_list[index];
    }
    inline const size_t size()const
    {
        return m_list.size();
    }

    StandardDeviation(const numt &scale = 1)
    {
        if (scale <= 0)throw Exception<StandardDeviation>("Uncertainty scaling factor must be greater than zero");
        m_scale = scale;
        m_sum = 0;
    }
    virtual ~StandardDeviation() {}
    StandardDeviation &operator<<(const numt &x)
    {
        m_list.push_back(x);
        m_sum += x;
        m_cache = nullptr;
        return *this;
    }
    inline const size_t count()const
    {
        return m_list.size();
    }
    inline const numt &scaling_factor()const
    {
        return m_scale;
    }
    const value<numt> &operator()()const
    {
        using namespace std;
        if (!m_cache) {
            size_t sz = m_list.size();
            if (sz <= 1)
                throw Exception<StandardDeviation>("No data to check. for sigma needed at least two elements.");
            numt average = m_sum / sz;
            numt m_sigsqr = 0;
            for (auto value : m_list)
                m_sigsqr += pow(value - average, 2);
            m_sigsqr /= sz - 1;
            const_cast<StandardDeviation &>(*this).m_cache =
                make_shared<value<numt>>(average, sqrt(m_sigsqr) * m_scale);
        }
        return *m_cache;
    }
};
template<typename numt = double>
class WeightedAverage
{
private:
    numt Sum;
    numt Wnorm;
    std::shared_ptr<value<numt>>m_cache;
public:
    WeightedAverage()
    {
        Sum = 0;
        Wnorm = 0;
    }
    virtual ~WeightedAverage() {}
    WeightedAverage &operator<<(const value<numt> &X)
    {
        if (X.uncertainty() == 0)
            throw Exception<WeightedAverage>("Cannot add value with zero error");
        numt w = 1.0 / pow(X.uncertainty(), 2);
        Sum += w * X.val();
        Wnorm += w;
        m_cache = nullptr;
        return *this;
    }
    const value<numt> &operator()()const
    {
        using namespace std;
        if (!m_cache) {
            if (Wnorm <= 0)throw Exception<WeightedAverage>("Attempt to check empty data");
            const_cast<WeightedAverage &>(*this).m_cache = make_shared<value<numt>>(Sum / Wnorm, 1.0 / sqrt(Wnorm));
        }
        return *m_cache;
    }
};
template<typename numt = double>
class CorrelationLinear
{
private:
    StandardDeviation<numt> _X, _Y, _XY;
public:
    CorrelationLinear(const numt &scale = 1):
        _X(scale), _Y(scale), _XY(scale) {}
    virtual ~CorrelationLinear() {}
    CorrelationLinear &operator<<(const std::pair<numt, numt> &P)
    {
        _X << P.first;
        _Y << P.second;
        _XY << P.first *P.second;
        return *this;
    }
    const numt Covariance()const
    {
        return _XY().val() - _X().val() * _Y().val();
    }
    const numt R()const
    {
        return Covariance() / (_X().uncertainty() * _Y().uncertainty());
    }
    const StandardDeviation<numt>&X()const{return _X;}
    const StandardDeviation<numt>&Y()const{return _Y;}
};
};
#endif
