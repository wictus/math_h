// this file is distributed under
// LGPLv3 license
#ifndef ___________VECTORS_H_____
#	define ___________VECTORS_H_____
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <tuple>
#include <math.h>
#include <type_traits>
#include "error.h"
#include "randomfunc.h"
#if __cplusplus>201700L
#define ____optimized_version_of_vectors_h_____
#else
#warning compiler does not support "if constexpr(...)". c++>=17 is needed. classes from vectors.h will work slower
#endif
namespace MathTemplates
{
template<size_t size = 3, class numt = double>class Vector;
template<size_t size = 3, class linetype = Vector<3, double>>class VectorTransformation;
template<size_t size = 3, class numt = double>class Direction;

template<class numt>
class Vector<1, numt>
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class VectorTransformation;
public:
    enum {Dimensions = 1};
    typedef numt NumberType;
    typedef Direction<1, numt> DType;
private:
    numt m_x;
protected:
    inline const numt &___last_component()const
    {
        return m_x;
    }
    inline Vector(const numt &x): m_x(x) {}
public:
    virtual ~Vector() {}
    inline const std::tuple<numt> to_tuple()const
    {
        return std::make_tuple(m_x);
    }
    template<class numt2>
    inline Vector(const Vector<Dimensions, numt2> &source): m_x(source.___last_component()) {}
    template<class... Args>
    inline Vector(const std::tuple<Args...> &v): m_x(std::get<0>(v)) {}
    inline static const Vector zero()
    {
        return Vector(std::make_tuple(numt(0)));
    }
#ifdef ____optimized_version_of_vectors_h_____
    template<size_t index>
    inline static const Vector basis_vector()
    {
	static_assert(index ==1,"dimension index is out of range");
        return Vector(std::make_tuple(numt(1)));
    }
    template<size_t index>
    inline const numt &component()const
    {
	static_assert(index == 1,"dimension index is out of range");
        return m_x;
    }
#else
    template<size_t index>
    static const Vector basis_vector()
    {
	static_assert(index>0,"dimension index is out of range");
	if(index>1)throw Exception<Vector>("dimension index is out of range");
        return Vector(std::make_tuple(numt(1)));
    }
    template<size_t index>
    const numt &component()const
    {
	static_assert(index>0,"dimension index is out of range");
	if(index>1)throw Exception<Vector>("dimension index is out of range");
        return m_x;
    }
#endif
    inline const numt &x()const
    {
        return m_x;
    }
    Vector &operator=(const Vector &source)
    {
        m_x = source.m_x;
        return *this;
    }
    Vector &operator+=(const Vector &second)
    {
        m_x += second.m_x;
        return *this;
    }
    const Vector operator+(const Vector &second)const
    {
        return Vector(std::make_tuple(m_x + second.m_x));
    }
    Vector &operator-=(const Vector &second)
    {
        m_x -= second.m_x;
        return *this;
    }
    const Vector operator-(const Vector &second)const
    {
        return Vector(std::make_tuple(m_x - second.m_x));
    }

    Vector &operator*=(const numt &second)
    {
        m_x *= second;
        return *this;
    }
    const Vector operator*(const numt &second)const
    {
        return Vector(std::make_tuple(m_x * second));
    }
    Vector &operator/=(const numt &second)
    {
        m_x /= second;
        return *this;
    }
    const Vector operator/(const numt &second)const
    {
        return Vector(std::make_tuple(m_x / second));
    }
    const numt operator*(const Vector &second)const
    {
        return m_x * second.m_x;
    }
    inline const numt M_sqr()const
    {
        return operator*(*this);
    }
    inline const numt M()const
    {
        return sqrt(M_sqr());
    }
    const bool operator==(const Vector &second)const
    {
        return (m_x == second.m_x);
    }
    const bool CloseTo(const Vector &second, const numt &epsilon)const
    {
        return operator-(second).M() < epsilon;
    }
};
template<size_t size, class numt>
class Vector
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class VectorTransformation;
public:
    enum {Dimensions = size};
    typedef numt NumberType;
    typedef Vector < size - 1, numt > VectorN;
    typedef Direction<size, numt> DType;
private:
    VectorN m_other;
    numt m_x;
    inline Vector(const VectorN &other, const numt &x): m_other(other), m_x(x) {}
protected:
    inline const VectorN &___recursive()const
    {
        return m_other;
    }
    inline const numt &___last_component()const
    {
        return m_x;
    }
public:
    virtual ~Vector() {}
    inline auto to_tuple()const->decltype(std::tuple_cat(m_other.to_tuple(), std::make_tuple(m_x)))
    {
        return std::tuple_cat(m_other.to_tuple(), std::make_tuple(m_x));
    }
    template<class numt2>
    inline Vector(const Vector<Dimensions, numt2> &source): m_other(source.___recursive()), m_x(source.___last_component()) {}
    template<class... Args>
    inline Vector(const std::tuple<Args...> &v): m_other(v), m_x(std::get < Dimensions - 1 > (v)) {}
    inline static const Vector zero()
    {
        return Vector(VectorN::zero(), numt(0));
    }
#ifdef ____optimized_version_of_vectors_h_____
    template<size_t index>
    inline static const Vector basis_vector()
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=Dimensions,"dimension index is out of range");
        if constexpr(index == Dimensions) return Vector(VectorN::zero(), numt(1));
	else return Vector(VectorN::template basis_vector<index>(), numt(0));
    }
    template<size_t index>
    inline const numt &component()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=Dimensions,"dimension index is out of range");
        if constexpr(index == Dimensions)return m_x;
	else return m_other.template component<index>();
    }
#else
    template<size_t index>
    static const Vector basis_vector()
    {
	static_assert(index > 0,"dimension index is out of range");
	if(index>Dimensions)throw Exception<Vector>("dimension index is out of range");
        if(index == Dimensions) return Vector(VectorN::zero(), numt(1));
	else return Vector(VectorN::template basis_vector<index>(), numt(0));
    }
    template<size_t index>
    const numt &component()const
    {
	static_assert(index > 0,"dimension index is out of range");
	if(index>Dimensions)throw Exception<Vector>("dimension index is out of range");
        if(index == Dimensions)return m_x;
	else return m_other.template component<index>();
    }
#endif
    inline const numt &x()const
    {
        return component<1>();
    }
    inline const numt &y()const
    {
        return component<2>();
    }
    inline const numt &z()const
    {
        return component<3>();
    }
    Vector &operator+=(const Vector &second)
    {
        m_x += second.m_x;
        m_other += second.m_other;
        return *this;
    }
    const Vector operator+(const Vector &second)const
    {
        return Vector(m_other + second.m_other, m_x + second.m_x);
    }
    Vector &operator-=(const Vector &second)
    {
        m_x -= second.m_x;
        m_other -= second.m_other;
        return *this;
    }
    const Vector operator-(const Vector &second)const
    {
        return Vector(m_other - second.m_other, m_x - second.m_x);
    }

    Vector &operator*=(const numt &second)
    {
        m_x *= second;
        m_other *= second.m_other;
        return *this;
    }
    const Vector operator*(const numt &second)const
    {
        return Vector(m_other * second, m_x * second);
    }
    Vector &operator/=(const numt &second)
    {
        m_x /= second;
        m_other /= second.m_other;
        return *this;
    }
    const Vector operator/(const numt &second)const
    {
        return Vector(m_other / second, m_x / second);
    }
    const numt operator*(const Vector &second)const
    {
        return (m_other * second.m_other) + (m_x * second.m_x);
    }
    inline const numt M_sqr()const
    {
        return operator*(*this);
    }
    inline const numt M()const
    {
        return sqrt(M_sqr());
    }
    const bool operator==(const Vector &second)const
    {
        return (m_x == second.m_x) && (m_other == second.m_other);
    }
    const bool CloseTo(const Vector &second, const numt &epsilon)const
    {
        return operator-(second).M() < epsilon;
    }
};
template<class numt, class... Args>
inline const Vector < sizeof...(Args) + 1, numt > desCartes(const numt &x, Args... args)
{
    return Vector < sizeof...(Args) + 1, numt > (std::make_tuple(x, args...));
}
template<size_t i, size_t d, class numt = double>
inline const Vector<d, numt> axis()
{
    return Vector<d, numt>::template basis_vector<i>();
}
template<class numt = double>
inline const Vector<2, numt> x()
{
    return Vector<2, numt>::template basis_vector<1>();
}
template<class numt = double>
inline const Vector<2, numt> y()
{
    return Vector<2, numt>::template basis_vector<2>();
}
template<class numt = double>
inline const Vector<2, numt> zero()
{
    return Vector<2, numt>::zero();
}
template<class numt = double>
inline const Vector<3, numt> X()
{
    return Vector<3, numt>::template basis_vector<1>();
}
template<class numt = double>
inline const Vector<3, numt> Y()
{
    return Vector<3, numt>::template basis_vector<2>();
}
template<class numt = double>
inline const Vector<3, numt> Z()
{
    return Vector<3, numt>::template basis_vector<3>();
}
template<class numt = double>
inline const Vector<3, numt> Zero()
{
    return Vector<3, numt>::zero();
}
template<size_t i, class numt = double>
inline const Vector<i, numt> operator-(const Vector<i, numt> &V)
{
    return V * numt(-1);
}





template<class linetype>
class VectorTransformation<1, linetype>
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class VectorTransformation;
public:
    enum {DimensionsFinal = 1};
    enum {DimensionsInitial = linetype::Dimensions};
    typedef typename linetype::NumberType NumberType;
    typedef linetype VIType;
    typedef Vector<1, NumberType> VFType;
    typedef Vector < DimensionsInitial + 1, NumberType > VIWType;
    typedef VectorTransformation<DimensionsFinal, VIWType> PlusOneColumn;
    typedef Vector < DimensionsInitial - 1, NumberType > VINType;
    typedef VectorTransformation<DimensionsFinal, VINType> MinusOneColumn;
private:
    VIType m_line;
protected:
    inline VectorTransformation(const VIType &line): m_line(line) {}
    inline const MinusOneColumn ___minus_one_column()const
    {
        return MinusOneColumn(m_line.___recursive());
    }
    inline const VFType ___last_column()const
    {
        return desCartes(m_line.___last_component());
    }
    inline const PlusOneColumn ___add_column(const VFType &col)const
    {
        return VIWType(m_line, col.___last_component());
    }
    inline const VIType &___last_line()const
    {
        return m_line;
    }
    inline const VectorTransformation &AddColumns()const
    {
        return *this;
    }
public:
    virtual ~VectorTransformation() {}
#ifdef ____optimized_version_of_vectors_h_____
    template<size_t index, size_t jindex>
    inline const NumberType &element()const
    {
	static_assert(index ==1 ,"dimension index is out of range");
        return m_line.template component<jindex>();
    }
#else
    template<size_t index, size_t jindex>
    const NumberType&element()const
    {
	static_assert(index>0, "dimension index is out of range");
	if(index>1)throw Exception<VectorTransformation>("dimension index out of range");
	return m_line.template component<jindex>();
    }
#endif
    const VFType operator*(const VIType &v)const
    {
        return desCartes(m_line * v);
    }
    const VectorTransformation<DimensionsFinal, Vector<1, NumberType>>
            operator*(const VectorTransformation<DimensionsInitial, Vector<1, NumberType>> &B)const
    {
        return desCartes(m_line * B.___last_column());
    }
    template<size_t third_size>
    const VectorTransformation<DimensionsFinal, Vector<third_size, NumberType>>
            operator*(const VectorTransformation<DimensionsInitial, Vector<third_size, NumberType>> &B)const
    {
        return operator*(B.___minus_one_column()).___add_column(operator*(B.___last_column()));
    }
    template<class otherlinetype>
    inline VectorTransformation(const VectorTransformation<DimensionsFinal, otherlinetype> &source): m_line(source.___last_line()) {}
    template<class... Args>
    inline VectorTransformation(const std::tuple<Args...> &v): m_line(std::get < DimensionsFinal - 1 > (v)) {}
    inline VectorTransformation(const VFType &A, const VIType &B): m_line(B *A.___last_component()) {}
    const bool operator==(const VectorTransformation &B)const
    {
        return m_line == B.m_line;
    }
    const VectorTransformation operator*(const NumberType &v)const
    {
        return VectorTransformation(m_line * v);
    }
    const VectorTransformation operator/(const NumberType &v)const
    {
        return VectorTransformation(m_line / v);
    }
    const VectorTransformation operator+(const VectorTransformation &B)const
    {
        return VectorTransformation(m_line + B.m_line);
    }
    const VectorTransformation operator-(const VectorTransformation &B)const
    {
        return VectorTransformation(m_line - B.m_line);
    }
    static inline const VectorTransformation zero()
    {
        return VectorTransformation(VIType::zero());
    }
    static inline const VectorTransformation one()
    {
        return VectorTransformation(VIType::template basis_vector<DimensionsFinal>());
    }
    template<size_t x, size_t y>
    static inline const VectorTransformation RotationInPlane(const NumberType &angle)
    {
        return VectorTransformation(
                   (x == DimensionsFinal) ? ((VIType::template basis_vector<x>() * cos(angle)) - (VIType::template basis_vector<y>() * sin(angle))) :
                   (y == DimensionsFinal) ? ((VIType::template basis_vector<y>() * cos(angle)) + (VIType::template basis_vector<x>() * sin(angle))) :
                   VIType::template basis_vector<DimensionsFinal>()
                                                );
    }
    template<class... Args>
    inline const VectorTransformation < DimensionsFinal, Vector < DimensionsInitial + 1 + sizeof...(Args), NumberType >> AddColumns(const VFType &col, Args...args)const
    {
        return ___add_column(col).AddColumns(args...);
    }
};
template<size_t sizef, class linetype>
class VectorTransformation
{
    template<size_t sizeff, class n>friend class Vector;
    template<size_t sizeff, class n>friend class Direction;
    template<size_t sizeff, class n>friend class VectorTransformation;
public:
    enum {DimensionsFinal = sizef};
    enum {DimensionsInitial = linetype::Dimensions};
    typedef typename linetype::NumberType NumberType;
    typedef linetype VIType;
    typedef Vector<sizef, NumberType> VFType;
    typedef VectorTransformation < sizef - 1, linetype > MinorTransformation;
    typedef Vector < DimensionsInitial + 1, NumberType > VIWType;
    typedef VectorTransformation<DimensionsFinal, VIWType> PlusOneColumn;
    typedef Vector < DimensionsInitial - 1, NumberType > VINType;
    typedef VectorTransformation<DimensionsFinal, VINType> MinusOneColumn;
private:
    MinorTransformation m_minor;
    VIType m_line;
protected:
    inline VectorTransformation(const MinorTransformation &minor, const VIType &line): m_minor(minor), m_line(line) {}
    inline const MinusOneColumn ___minus_one_column()const
    {
        const auto new_minor = m_minor.___minus_one_column();
        const auto new_line = m_line.___recursive();
        return MinusOneColumn(new_minor, new_line);
    }
    inline const VFType ___last_column()const
    {
        return VFType(m_minor.___last_column(), m_line.___last_component());
    }
    inline const PlusOneColumn ___add_column(const VFType &col)const
    {
        const auto new_minor = m_minor.___add_column(col.___recursive());
        const auto new_line = VIWType(m_line, col.___last_component());
        return PlusOneColumn(new_minor, new_line);
    }
    inline const VIType &___last_line()const
    {
        return m_line;
    }
    inline const MinorTransformation &___recursive()const
    {
        return m_minor;
    }
    inline const VectorTransformation &AddColumns()const
    {
        return *this;
    }
public:
    virtual ~VectorTransformation() {}
#ifdef ____optimized_version_of_vectors_h_____
    template<size_t index, size_t jindex>
    inline const NumberType &element()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=DimensionsFinal,"dimension index is out of range");
	if constexpr(index==DimensionsFinal) return m_line.template component<jindex>();
	else return m_minor.template element<index, jindex>();
    }
#else
    template<size_t index, size_t jindex>
    const NumberType &element()const
    {
	static_assert(index > 0,"dimension index is out of range");
	if(index>DimensionsFinal)throw Exception<VectorTransformation>("dimension index is out of range");
	if(index==DimensionsFinal) return m_line.template component<jindex>();
	else return m_minor.template element<index, jindex>();
    }
#endif
    const VFType operator*(const VIType &v)const
    {
        return VFType(m_minor * v, m_line * v);
    }
    const VectorTransformation<DimensionsFinal, Vector<1, NumberType>>
            operator*(const VectorTransformation<DimensionsInitial, Vector<1, NumberType>> &B)const
    {
        const auto P = m_minor * B;
        const auto C = desCartes(m_line * B.___last_column());
        return VectorTransformation<DimensionsFinal, Vector<1, NumberType>>(P, C);
    }
    template<size_t third_size>
    const VectorTransformation<DimensionsFinal, Vector<third_size, NumberType>>
            operator*(const VectorTransformation<DimensionsInitial, Vector<third_size, NumberType>> &B)const
    {
        const auto P = operator*(B.___minus_one_column());
        const VFType C = operator*(B.___last_column());
        return P.___add_column(C);
    }
    template<class otherlinetype>
    inline VectorTransformation(const VectorTransformation<DimensionsFinal, otherlinetype> &source): m_minor(source.___recursive()), m_line(source.___last_line()) {}
    template<class... Args>
    inline VectorTransformation(const std::tuple<Args...> &v): m_minor(v), m_line(std::get < DimensionsFinal - 1 > (v)) {}
    inline VectorTransformation(const VFType &A, const VIType &B): m_minor(A.___recursive(), B), m_line(B *A.___last_component()) {}
    const bool operator==(const VectorTransformation &B)const
    {
        return (m_line == B.m_line) && (m_minor == B.m_minor);
    }
    const VectorTransformation operator*(const NumberType &v)const
    {
        return VectorTransformation(m_minor * v, m_line * v);
    }
    const VectorTransformation operator/(const NumberType &v)const
    {
        return VectorTransformation(m_minor / v, m_line / v);
    }
    const VectorTransformation operator+(const VectorTransformation &B)const
    {
        return VectorTransformation(m_minor + B.m_minor, m_line + B.m_line);
    }
    const VectorTransformation operator-(const VectorTransformation &B)const
    {
        return VectorTransformation(m_minor - B.m_minor, m_line - B.m_line);
    }
    static inline const VectorTransformation zero()
    {
        return VectorTransformation(MinorTransformation::zero(), VIType::zero());
    }
    static inline const VectorTransformation one()
    {
        return VectorTransformation(MinorTransformation::one(), VIType::template basis_vector<DimensionsFinal>());
    }
    template<size_t x, size_t y>
    static inline const VectorTransformation RotationInPlane(const NumberType &angle)
    {
        return VectorTransformation(
                   MinorTransformation::template RotationInPlane<x, y>(angle),
                   (x == DimensionsFinal) ? ((VIType::template basis_vector<x>() * cos(angle)) - (VIType::template basis_vector<y>() * sin(angle))) :
                   (y == DimensionsFinal) ? ((VIType::template basis_vector<y>() * cos(angle)) + (VIType::template basis_vector<x>() * sin(angle))) :
                   VIType::template basis_vector<DimensionsFinal>()
                                                );
    }
    template<class... Args>
    inline const VectorTransformation < DimensionsFinal, Vector < DimensionsInitial + 1 + sizeof...(Args), NumberType >> AddColumns(const VFType &col, Args...args)const
    {
        return ___add_column(col).AddColumns(args...);
    }
};
template<class linetype, class... Args>
inline const VectorTransformation < sizeof...(Args) + 1, linetype > lines(const linetype &x, Args... args)
{
    return VectorTransformation < sizeof...(Args) + 1, linetype > (std::make_tuple(x, args...));
}
template<class numt, class... Args>
inline const VectorTransformation < 1, Vector < sizeof...(Args) + 1, numt >> line(const numt &x, Args... args)
{
    return lines(desCartes(x, args...));
}
template<class numt, class...Args>
inline const VectorTransformation < 1 + sizeof...(Args), Vector<1, numt >> column(const numt &x, Args...args)
{
    return VectorTransformation < 1 + sizeof...(Args), Vector<1, numt >> (std::make_tuple(x, args...));
}
template<class VecT, class...Args>
inline const VectorTransformation < VecT::Dimensions, Vector < 1 + sizeof...(Args), typename VecT::NumberType >> columns(const VecT &x, Args...args)
{
    return VectorTransformation<VecT::Dimensions, Vector<1, typename VecT::NumberType>>(x.to_tuple()).AddColumns(args...);
}
template<size_t s, class numt = double>
inline const VectorTransformation<s, Vector<s, numt>> ZERO()
{
    return VectorTransformation<s, Vector<s, numt>>::zero();
}
template<size_t s, class numt = double>
inline const VectorTransformation<s, Vector<s, numt>> ONE()
{
    return VectorTransformation<s, Vector<s, numt>>::one();
}
template<size_t size, class VIType>
inline const VectorTransformation<size, VIType> TensorProduct(const Vector<size, typename VIType::NumberType> &A, const VIType &B)
{
    return VectorTransformation<size, VIType>(A, B);
}
#define AC(n) (A.template component<n>())
#define ZeRo numt(0)
template<class numt = double>
inline const VectorTransformation<1, Vector<2, numt>> SkewM(const Vector<2, numt> &A)
{
    return lines(desCartes(-AC(2), AC(1)));
}
template<class numt = double>
inline const VectorTransformation<3, Vector<3, numt>> SkewM(const Vector<3, numt> &A)
{
    return lines(
               desCartes(ZeRo, -AC(3), AC(2)),
               desCartes(AC(3),  ZeRo, -AC(1)),
               desCartes(-AC(2), AC(1),  ZeRo)
           );
}
template<class numt = double>
inline const VectorTransformation<7, Vector<7, numt>> SkewM(const Vector<7, numt> &A)
{
    return lines(
               desCartes(ZeRo, -AC(4), -AC(7), AC(2), -AC(6), AC(5), AC(3)),
               desCartes(AC(4),  ZeRo, -AC(5), -AC(1), AC(3), -AC(7), AC(6)),
               desCartes(AC(7), AC(5),  ZeRo, -AC(6), -AC(2), AC(4), -AC(1)),
               desCartes(-AC(2), AC(1), AC(6),  ZeRo, -AC(7), -AC(3), AC(5)),
               desCartes(AC(6), -AC(3), AC(2), AC(7),  ZeRo, -AC(1), -AC(4)),
               desCartes(-AC(5), AC(7), -AC(4), AC(3), AC(1),  ZeRo, -AC(2)),
               desCartes(-AC(3), -AC(6), AC(1), -AC(5), AC(4), AC(2), ZeRo)
           );
}
#undef ZeRo
#undef AC
template<class A, class B>
inline auto operator^(const A &a, const B &b)->decltype(SkewM(a)*b)
{
    return SkewM(a) * b;
}








template<class numt>
class Direction<1, numt>
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class VectorTransformation;
private:
    bool sign;
    static const numt PHI(RANDOM &r)
    {
        static const RandomUniform<numt> res(-PI<numt>(), +PI<numt>());
        return res(r);
    }
protected:
    inline const bool &___sign()const
    {
        return sign;
    }
public:
    enum {Dimensions = 1};
    enum {Thetas = 0};
    typedef numt NumberType;
    typedef Vector<1, numt> VType;
    virtual ~Direction() {}
    inline Direction(RANDOM &RG): sign(PHI(RG) >= 0) {}
    template<class numt2>
    inline Direction(const Direction<Dimensions, numt2> &source): sign(source.___sign()) {}
    inline Direction(const VType &v): sign(v.x() >= 0) {}
    template<class... Args>
    inline Direction(const std::tuple<Args...> &args): sign(std::get<0>(args)) {}
    inline const VType operator*(const numt &rho)const
    {
        if (sign)return desCartes(rho);
        else return desCartes(-rho);
    }
    inline const bool operator==(const Direction &second)const
    {
        return (sign == second.sign);
    }
    inline const numt dir()const
    {
        return sign ? numt(1) : numt(-1);
    }
    inline const VectorTransformation<Dimensions, VType> Rotations()const
    {
        return sign ? lines(desCartes(numt(1))) : lines(desCartes(numt(-1)));
    }
    inline const VectorTransformation<Dimensions, VType> AntiRotations()const
    {
        return sign ? lines(desCartes(numt(1))) : lines(desCartes(numt(-1)));
    }
};

template<class numt>
class Direction<2, numt>
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class VectorTransformation;
public:
    enum {Dimensions = 2};
    enum {Thetas = 0};
    typedef numt NumberType;
    typedef Vector<2, numt> VType;
private:
    numt m_phi;
    void ___phi_norm()
    {
        while (m_phi > PI<numt>())m_phi -= PI<numt>() * numt(2);
        while (m_phi < -PI<numt>())m_phi += PI<numt>() * numt(2);
    }
    static const numt PHI(RANDOM &r)
    {
        static const RandomUniform<numt> res(-PI<numt>(), +PI<numt>());
        return res(r);
    }
    inline const numt &___last_component()const
    {
        return m_phi;
    }
public:
    virtual ~Direction() {}
    template<class... Args>
    inline Direction(const std::tuple<Args...> &v): m_phi(std::get<0>(v))
    {
        ___phi_norm();
    }
    inline Direction(RANDOM &RG): m_phi(PHI(RG)) {}
    template<class numt2>
    inline Direction(const Direction<Dimensions, numt2> &source): m_phi(source.___last_component())
    {
        ___phi_norm();
    }
    inline const numt &phi()const
    {
        return m_phi;
    }
    Direction(const VType &V): m_phi(atan2(V.y(), V.x())) {}
    const VType operator*(const numt &rho)const
    {
        return desCartes(cos(m_phi), sin(m_phi)) * rho;
    }
    const bool operator==(const Direction &second)const
    {
        return (m_phi == second.m.phi);
    }
    const VectorTransformation<Dimensions, VType> Rotations()const
    {
        return VectorTransformation<Dimensions, VType>::template RotationInPlane<1, 2>(m_phi);
    }
    const VectorTransformation<Dimensions, VType> AntiRotations()const
    {
        return VectorTransformation<Dimensions, VType>::template RotationInPlane<1, 2>(-m_phi);
    }
};
template<class numt>
class Direction<3, numt>
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class VectorTransformation;
public:
    enum {Dimensions = 3};
    enum {Thetas = 1};
    typedef numt NumberType;
    typedef Vector<3, numt> VType;
    typedef Direction<2, numt> DirectionN;
private:
    DirectionN m_ld;
    numt m_theta;
    void ___theta_check()const
    {
        if (m_theta < numt(0))throw Exception<Direction>("wrong theta value <0");
        if (m_theta > PI<numt>())throw Exception<Direction>("wrong theta value >pi");
    }
    static const numt CTHETA(RANDOM &r)
    {
        static const RandomUniform<numt> res(numt(-1), numt(+1));
        return res(r);
    }
protected:
    inline const DirectionN &___recursive()const
    {
        return m_ld;
    }
    inline const numt &___last_component()const
    {
        return m_theta;
    }
public:
    virtual ~Direction() {}
    template<class numt2>
    inline Direction(const Direction<Dimensions, numt2> &source): m_ld(source.___recursive()), m_theta(source.___last_component())
    {
        ___theta_check();
    }
    template<class... Args>
    inline Direction(const std::tuple<Args...> &v): m_ld(v), m_theta(std::get<1>(v))
    {
        ___theta_check();
    }
    inline Direction(RANDOM &RG): m_ld(RG), m_theta(acos(CTHETA(RG))) {}
    inline const numt &phi()const
    {
        return m_ld.phi();
    }
#ifdef ____optimized_version_of_vectors_h_____
    template<size_t index = 1>
    inline const numt &th()const
    {
	static_assert(index ==1,"dimension index is out of range");
        return m_theta;
    }
#else
    template<size_t index = 1>
    const numt &th()const
    {
	static_assert(index >0,"dimension index is out of range");
	if(index>1)throw Exception<Direction>("dimension index is out of range");
        return m_theta;
    }
#endif
    Direction(const VType &V): m_ld(V.___recursive()), m_theta(acos(V.template component<Dimensions>() / V.M())) {}
    const VType operator*(const numt &rho)const
    {
        return VType(m_ld * (rho * sin(m_theta)), rho * cos(m_theta));
    }
    const bool operator==(const Direction &second)const
    {
        return (m_theta == second.m_theta) && (m_ld == second.m_ld);
    }
    const VectorTransformation<Dimensions, VType> Rotations()const
    {
        return
            VectorTransformation<Dimensions, VType>(
                m_ld.Rotations().___add_column(DirectionN::VType::zero()),
                VType::template basis_vector<Dimensions>())
            * VectorTransformation<Dimensions, VType>::template RotationInPlane<3, 1>(m_theta);
    }
    const VectorTransformation<Dimensions, VType> AntiRotations()const
    {
        return
            VectorTransformation<Dimensions, VType>::template RotationInPlane<3, 1>(-m_theta) *
        VectorTransformation<Dimensions, VType>(
            m_ld.AntiRotations().___add_column(DirectionN::VType::zero()),
            VType::template basis_vector<Dimensions>());
    }
};
template<size_t size, class numt>
class Direction
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class VectorTransformation;
public:
    enum {Dimensions = size};
    enum {Thetas = size - 2};
    typedef numt NumberType;
    typedef Vector<size, numt> VType;
    typedef Direction < size - 1, numt > DirectionN;
private:
    DirectionN m_ld;
    numt m_theta;
    void ___theta_check()const
    {
        if (m_theta < numt(0))throw Exception<Direction>("wrong theta value <0");
        if (m_theta > PI<numt>())throw Exception<Direction>("wrong theta value >pi");
    }
    static const numt CTHETA(RANDOM &r)
    {
        static const RandomUniform<numt> res(numt(-1), numt(+1));
        return res(r);
    }
protected:
    inline const DirectionN &___recursive()const
    {
        return m_ld;
    }
    inline const numt &___last_component()const
    {
        return m_theta;
    }
public:
    virtual ~Direction() {}
    template<class numt2>
    inline Direction(const Direction<Dimensions, numt2> &source): m_ld(source.___recursive()), m_theta(source.___last_component())
    {
        ___theta_check();
    }
    template<class... Args>
    inline Direction(const std::tuple<Args...> &v): m_ld(v), m_theta(std::get<Thetas>(v))
    {
        ___theta_check();
    }
    inline Direction(RANDOM &RG): m_ld(RG), m_theta(acos(CTHETA(RG))) {}
    inline const numt &phi()const
    {
        return m_ld.phi();
    }
#ifdef ____optimized_version_of_vectors_h_____
    template<size_t index>
    inline const numt &th()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index <= Thetas,"dimension index is out of range");
        if constexpr(index == Thetas) return m_theta;
	else return m_ld.template th<index>();
    }
#else
    template<size_t index>
    const numt &th()const
    {
	static_assert(index > 0,"dimension index is out of range");
	if(index > Thetas)throw Exception<Direction>("dimension index is out of range");
        if(index == Thetas) return m_theta;
	else return m_ld.template th<index>();
    }
#endif
    Direction(const VType &V): m_ld(V.___recursive()), m_theta(acos(V.template component<Dimensions>() / V.M())) {}
    const VType operator*(const numt &rho)const
    {
        return VType(m_ld * (rho * sin(m_theta)), rho * cos(m_theta));
    }
    const VectorTransformation<Dimensions, VType> Rotations()const
    {
        return
            VectorTransformation<Dimensions, VType>(
                m_ld.Rotations().___add_column(DirectionN::VType::zero()),
                VType::template basis_vector<Dimensions>())
            * VectorTransformation<Dimensions, VType>::template RotationInPlane < Dimensions, Dimensions - 1 > (m_theta);
    }
    const VectorTransformation<Dimensions, VType> AntiRotations()const
    {
        return
            VectorTransformation<Dimensions, VType>::template RotationInPlane < Dimensions, Dimensions - 1 > (-m_theta) *
        VectorTransformation<Dimensions, VType>(
            m_ld.AntiRotations().___add_column(DirectionN::VType::zero()),
            VType::template basis_vector<Dimensions>());
    }
};
template<class numt = double>
inline const Direction< 1, numt > direction()
{
    return Direction < 1, numt > (std::make_tuple(true));
}
template<class numt = double, class... Args>
inline const Direction < 2 + sizeof...(Args), numt > direction(const numt &phi, Args... other)
{
    return Direction < 2 + sizeof...(Args), numt > (std::make_tuple(phi, other...));
}
template<size_t size, class numt = double>
inline const Direction< size, numt > direction(const Vector<size, numt> &V)
{
    return Direction < size, numt > (V);
}
template<size_t size = 3, class numt = double, class RG = RANDOM>
inline const Direction<size, numt> randomIsotropic(RG &generator)
{
    return Direction<size, numt>(generator);
}
template<class VType>
struct VectorDecomposition {
    VType tau;
    VType n;
};
template<class VType>
inline const VectorDecomposition<VType> decompose_by_direction(const VType &source, const typename VType::DType &dir)
{
    typedef typename VType::NumberType numt;
    const auto t = dir * ((dir * numt(1)) * source);
    return {.tau = t, .n = source - t};
}
template<class VType>
inline const VectorDecomposition<VType> decompose_by_plane_normale(const VType &source, const typename VType::DType &pn)
{
    typedef typename VType::NumberType numt;
    const auto t = pn * ((pn * numt(1)) * source);
    return {.tau = source - t, .n = t};
}
template<class numt = double>
inline const VectorTransformation<2, Vector<2, numt>> Rotation(const numt &theta)
{
    const numt cost = cos(theta), sint = sin(theta);
    return lines(
               desCartes(cost, -sint),
               desCartes(sint, cost)
           );
}
template<class numt = double>
inline const VectorTransformation<3, Vector<3, numt>> Rotation(const Direction<3, numt> &axis, const numt &theta)
{
    const auto n = axis * numt(1);
    const numt cost = cos(theta), sint = sin(theta), one = 1;
    return lines(
               desCartes(cost + (one - cost) * n.x() * n.x(),	(one - cost) * n.x() * n.y() - sint * n.z(),	(one - cost) * n.x() * n.z() + sint * n.y()),
               desCartes((one - cost) * n.y() * n.x() + sint * n.z(),	    cost + (one - cost) * n.y() * n.y(),	(one - cost) * n.y() * n.z() - sint * n.x()),
               desCartes((one - cost) * n.z() * n.x() - sint * n.y(),	(one - cost) * n.z() * n.y() + sint * n.x(),	    cost + (one - cost) * n.z() * n.z())
           );
}












template<class numt = double, class Space = Vector<3, numt>>
class LorentzVector
{
public:
    typedef Space SpaceVectorType;
    typedef numt TimeCoordinateType;
private:
    TimeCoordinateType m_time;
    SpaceVectorType m_space;
public:
    virtual ~LorentzVector() {}
    inline LorentzVector(const numt &t, const Space &S): m_time(t), m_space(S) {}
    inline const TimeCoordinateType &E()const{return m_time;}
    inline const SpaceVectorType &P()const{return m_space;}
    template<class numt2, class Space2>
    inline LorentzVector(const LorentzVector<numt2, Space2> &source): m_time(source.E()), m_space(source.P()) {}
    LorentzVector &operator=(const LorentzVector &source)
    {
        m_space = source.m_space;
        m_time = source.m_time;
        return *this;
    }
    const bool operator==(const LorentzVector &second)const
    {
        return (m_time == second.m_time) && (m_space == second.m_space);
    }
    LorentzVector &operator+=(const LorentzVector &second)
    {
        m_time += second.m_time;
        m_space += second.m_space;
        return *this;
    }
    const LorentzVector operator+(const LorentzVector &second)const
    {
        return LorentzVector(m_time + second.m_time, m_space + second.m_space);
    }
    LorentzVector &operator-=(const LorentzVector &second)
    {
        m_time -= second.m_time;
        m_space -= second.m_space;
        return *this;
    }
    const LorentzVector operator-(const LorentzVector &second)const
    {
        return LorentzVector(m_time - second.m_time, m_space - second.m_space);
    }
    const numt operator*(const LorentzVector &second)const
    {
        return (m_time * second.m_time) - (m_space * second.m_space);
    }
    inline const numt M_sqr()const
    {
        return operator*(*this);
    }
    inline const numt M()const
    {
        return sqrt(M_sqr());
    }
    inline const numt Ekin()const
    {
	return E()-M();
    }

    inline static const LorentzVector zero()
    {
        return LorentzVector(numt(0), Space::zero());
    }
    template<class...Args>
    inline const LorentzVector Rotate(Args...args)const
    {
        return LorentzVector(E(), MathTemplates::Rotation(args...) * P());
    }
    const LorentzVector Transform(const Space &Beta)const
    {
        const numt beta = Beta.M();
        if (beta == 0.0)return *this;
        if (beta >= numt(1))throw Exception<LorentzVector>("Bad Lorentz transformation");
        const auto bn = Beta / beta;
        const numt gamma = numt(1) / sqrt(numt(1) - beta * beta);
        const auto ST = ONE<Space::Dimensions>() + TensorProduct(bn, bn) * (gamma - numt(1));
        const auto TT = -Beta * gamma;
        return LorentzVector((E() * gamma) + (P() * TT), (ST * P()) + (TT * E()));
    }
    const Space Beta()const
    {
        return P() / E();
    }
};
template<class numt = double, class Space = Vector<3, numt>>
inline const LorentzVector<numt, Space> lorentzVector(const numt &t, const Space &s)
{
    return LorentzVector<numt, Space>(t, s);
}
template<class numt = double, class Space = Vector<3, numt>>
inline const LorentzVector<numt, Space> lorentz_byPM(const Space &s, const numt &l4)
{
    return LorentzVector<numt, Space>(sqrt(s.M_sqr() + l4 * l4), s);
}
template<class numt = double, class Space = Vector<3, numt>>
inline const LorentzVector<numt, Space> lorentz_byEM(const numt &t, const numt &l4, const typename Space::DType &dir)
{
    numt Sp = sqrt(t * t - l4 * l4);
    return LorentzVector<numt, Space>(t, dir * Sp);
}
template<class numt = double, class Space = Vector<3, numt>>
inline const LorentzVector<numt, Space> lorentz_byEM(const numt &t, const numt &l4, const Space &Dir)
{
    return lorentz_byEM(t, l4, direction(Dir));
}
template<class numt = double, class Space = Vector<3, numt>>
inline const LorentzVector<numt, Space> lorentz_byEkM(const numt &e, const numt &l4, const typename Space::DType &dir)
{
    return lorentz_byEM(e+l4, l4, dir);
}
template<class numt = double, class Space = Vector<3, numt>>
inline const LorentzVector<numt, Space> lorentz_byEkM(const numt &e, const numt &l4, const Space &Dir)
{
    return lorentz_byEM(e+l4, l4, direction(Dir));
}
template<size_t size, class numt>
const std::pair<LorentzVector<numt, Vector<size, numt>>, LorentzVector<numt, Vector<size, numt>>>
        binaryDecay(const numt &IM, const numt &m1, const numt &m2, const Direction<size, numt> &dir)
{
    if (m1 < 0)throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Negative mass1 error");
    if (m2 < 0)throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Negative mass2 error");
    if (IM < (m1 + m2))throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Invariant mass of decaying system is less then masses of products");
    const auto E1 = (IM * IM - m2 * m2 + m1 * m1) / (IM * numt(2));
    const auto p = sqrt(E1 * E1 - m1 * m1);
    const auto P = dir * p;
    return std::make_pair(lorentz_byPM(P, m1), lorentz_byPM(-P, m2));
}

};
#undef ____optimized_version_of_vectors_h_____
#endif
