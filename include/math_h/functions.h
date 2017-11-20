// this file is distributed under
// LGPLv3 license
#ifndef EMEFWAYNIGJGENCP
#define EMEFWAYNIGJGENCP
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <math.h>
#include <functional>
#include <memory>
namespace MathTemplates
{
template<class result, typename... Args>
class IFunction
{
public:
    virtual ~IFunction() {}
    virtual result operator()(Args...)const = 0;
    inline const std::function<result(Args...)> func()const
    {
        return [this](Args... args)->result {return this->operator()(args...);};
    }
};
template<class result, typename... Args>
class StdFunctionWrap: public IFunction<result, Args...>
{
private:
    std::function<result(Args...)> m_func;
public:
    StdFunctionWrap(const std::function<result(Args...)>f)
    {
        m_func = f;
    }
    StdFunctionWrap(const StdFunctionWrap &source)
    {
        m_func = source.m_func;
    }
    virtual ~StdFunctionWrap() {}
    virtual result operator()(Args...args)const
    {
	return m_func(args...);
    };
};
template<class result, typename... Args>
inline const StdFunctionWrap<result, Args...>
Function(const std::function<result(Args...)>f)
{
    return StdFunctionWrap<result, Args...>(f);
}
template<class result, typename... Args>
inline const std::shared_ptr<StdFunctionWrap<result, Args...>>
        PFunction(const std::function<result(Args...)>f)
{
    return std::make_shared<StdFunctionWrap<result, Args...>>(f);
}
template<class numt = double>
inline const numt PI()
{
    return 3.1415926;
}
template<class numt = double>
numt Gaussian(const numt &x, const numt &X_max, const numt &sigma)
{
    return exp(-pow((X_max - x) / sigma, 2) / 2) / (sigma * sqrt(2 * PI<numt>()));
}
template<class numt = double>
numt Lorentzian(const numt &x, const numt &X_max, const numt &gamma)
{
    return gamma / (PI<numt>() * (pow(x - X_max, 2) + pow(gamma, 2)));
}
template<class numt = double>
numt BreitWigner(const numt &x, const numt &pos, const numt &gamma)
{
    return Lorentzian<numt>(x, pos, gamma / numt(2));
}
template<class numt = double>
numt Novosibirsk(const numt &x, const numt &pos, const numt &sigma, const numt &asym)
{
    if (pow(asym / sigma, 2) <= 0.000001)
        return Gaussian<numt>(x, pos, sigma);
    numt Lsqlog4 = asym * sqrt(log(4));
    numt qb = sinh(Lsqlog4) / Lsqlog4;
    numt q4 = qb * asym * (x - pos) / sigma + 1;
    if (q4 <= 0)return 0;
    numt slq4 = pow(log(q4), 2);
    numt k = sigma * sqrt(2 * PI<numt>() * (pow(asym, 2) + 1.0));
    return exp(-slq4 / (pow(asym, 2) * 2) + pow(asym, 2)) / k;
}
template<class numt = double>
numt FermiFunc(const numt &x, const numt &X_border, const numt &diffuse)
{
    return 1.0 / (1.0 + exp((x - X_border) / diffuse));
}
template<class numt, class indexer>
numt Polynom(const numt &x, const indexer p, const unsigned int P, const int index_offset = 0)
{
    numt res = 0;
    numt c = 1;
    for (unsigned int i = 0; i <= P; i++) {
        res += c * p[index_offset + i];
        if (i < P)c *= x;
    }
    return res;
}
template<unsigned int P, class numt, class indexer, int index_offset = 0>
inline numt Polynom(const numt &x, const indexer p)
{
    return Polynom<numt, indexer>(x, p, P, index_offset);
}
};
#endif
