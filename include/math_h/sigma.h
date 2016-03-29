// this file is distributed under 
// MIT license
#ifndef YJIOGPKSIIYJVKND
#define YJIOGPKSIIYJVKND
#include <iostream>
#include <list>
#include <vector>
#include <functional>
#include <math.h>
#include "error.h"
namespace MathTemplates{
	template<typename numt>
	class value{
	private:
		numt Value,Error;
	public:
		value():Value(numt(0)),Error(numt(0)){}
		value(const numt v):Value(v),Error(numt(0)){}
		value(const numt v,const numt err):Value(v),Error(err){
			if(Error<0)throw Exception<value>("Error cannot be negative");
		}
		static value std_error(const numt v){
			if(v<0)throw Exception<value>("Cannot calculate std error for negative value");
			auto res=value(v,sqrt(v));
			if(res.Error<numt(1))res.Error=numt(1);
			return res;
		}
		value(const value&source):Value(source.Value),Error(source.Error){}
		value&operator=(const value&source){
			Value=source.Value;
			Error=source.Error;
			return *this;
		}
		
		numt val()const{return Value;}
		numt delta()const{return Error;}
		numt epsilon()const{return Error/Value;}
		numt min()const{return Value-Error;}
		numt max()const{return Value+Error;}
		bool contains(const numt& x)const{return (x>=min())&&(x<=max());}
		bool contains(const value&x)const{return (x.max()>=min())&&(x.min()<=max());}
		bool contains(const numt&&x)const{return contains(x);}
		bool contains(const value&&x)const{return contains(x);}
		
		value&operator+=(const value&other){
			Error=sqrt(pow(Error,2)+pow(other.Error,2));
			Value+=other.Value;
			return *this;
		}
		value&operator+=(const value&&other){
			return operator+=(other);
		}
		value&operator-=(const value&other){
			Error=sqrt(pow(Error,2)+pow(other.Error,2));
			Value-=other.Value;
			return *this;
		}
		value&operator-=(const value&&other){
			return operator-=(other);
		}
		value&operator*=(const value&other){
			Error=sqrt(pow(Error*other.Value,2)+pow(other.Error*Value,2));
			Value*=other.Value;
			return *this;
		}
		value&operator*=(const value&&other){
			return operator*=(other);
		}
		value&operator/=(const value&other){
			Error=sqrt(pow(Error/other.Value,2)+pow(other.Error*Value/pow(other.Value,2),2));
			Value/=other.Value;
			return *this;
		}
		value&operator/=(const value&&other){
			return operator/=(other);
		}
		
		const value operator+(const value&other)const{return value(*this)+=other;}
		const value operator+(const value&&other)const{return value(*this)+=other;}
		const value operator-(const value&other)const{return value(*this)-=other;}
		const value operator-(const value&&other)const{return value(*this)-=other;}
		const value operator*(const value&other)const{return value(*this)*=other;}
		const value operator*(const value&&other)const{return value(*this)*=other;}
		const value operator/(const value&other)const{return value(*this)/=other;}
		const value operator/(const value&&other)const{return value(*this)/=other;}
		const value func(const std::function<numt(const numt&)>F)const{
			numt V=F(val());
			return value(V,sqrt( (pow(F(min())-V,2)+pow(F(max())-V,2)) / numt(2) ));
		}
		
		bool operator<(const value&other)const{return Value<other.Value;}
		bool operator<(const value&&other)const{return Value<other.Value;}
		bool operator>(const value&other)const{return Value>other.Value;}
		bool operator>(const value&&other)const{return Value>other.Value;}
		bool operator==(const value&other)const{return Value==other.Value;}
		bool operator==(const value&&other)const{return Value==other.Value;}
		bool operator>=(const value&other)const{return Value>=other.Value;}
		bool operator>=(const value&&other)const{return Value>=other.Value;}
		bool operator<=(const value&other)const{return Value<=other.Value;}
		bool operator<=(const value&&other)const{return Value<=other.Value;}
	};
	template<typename numt>
	std::ostream&operator<<(std::ostream&str,const value<numt>&P){
		return str<<P.val()<<"+/-"<<P.delta();
	}
	template<typename numt>
	std::ostream&operator<<(std::ostream&str,const value<numt>&&P){return str<<P;}
	template<typename numt>
	const value<numt> func_value(const std::function<numt(const numt&)>F,const value<numt>&X){return X.func(F);}
	template<typename numt>
	inline const value<numt> func_value(const std::function<numt(const numt&)> F,const value<numt>&&X){return X.func(F);}
	template<typename numt>
	const std::function<value<numt>(const value<numt>&)> value_func(const std::function<numt(const numt&)>F){
		return [F](const value<numt>&X)->value<numt>{return X.func(F);};
	}
	
	template<typename numt>
	class StandardDeviation{
	private:
		std::list<numt> m_list;
		numt m_sum;
		value<numt>*m_cache;
		numt m_scale;
	public:
		StandardDeviation(const numt scale=1){
			if(scale<=0)throw Exception<StandardDeviation>("Uncertainty scaling factor must be greater than zero");
			m_scale=scale;
			m_sum=0;
			m_cache=new value<numt>(INFINITY,INFINITY);
		}
		virtual ~StandardDeviation(){
			delete m_cache;
		}
		StandardDeviation &operator<<(const numt x){
			m_list.push_back(x);
			m_sum+=x;
			(*m_cache)=value<numt>(INFINITY,INFINITY);
			return *this;
		}
		const int count()const{return m_list.size();}
		const numt scaling_factor()const{return m_scale;}
		const value<numt>&operator()()const{
			using namespace std;
			if(isfinite(m_cache->val())){
				return *m_cache;
			}else{
				int sz=m_list.size();
				if(sz<=1)
					throw Exception<StandardDeviation>("No data to check. for sigma needed at least two elements.");
				numt average=m_sum/sz;
				numt m_sigsqr=0;
				for(auto value:m_list)
					m_sigsqr+=pow(value-average,2);
				m_sigsqr/=sz-1;
				(*m_cache)=value<numt>(average,sqrt(m_sigsqr)*m_scale);
				return *m_cache;
			}
		}
	};
	template<typename numt>
	class WeightedAverageCalculator{
	private:
		numt Sum;
		numt Wnorm;
		value<numt>*m_cache;
	public:
		WeightedAverageCalculator(){
			Sum=0;Wnorm=0;
			m_cache=new value<numt>(INFINITY,INFINITY);
		}
		virtual ~WeightedAverageCalculator(){delete m_cache;}
		WeightedAverageCalculator &operator<<(const value<numt>&X){
			if(X.delta()<=0)
				throw Exception<WeightedAverageCalculator>("Cannot add value with zero error");
			numt w=1.0/pow(X.delta(),2);
			Sum+=w*X.val();
			Wnorm+=w;
			(*m_cache)=value<numt>(INFINITY,INFINITY);
			return *this;
		}
		WeightedAverageCalculator &operator<<(const value<numt>&&X){
			return operator<<(X);
		}
		const value<numt>&operator()()const{
			using namespace std;
			if(isfinite(m_cache->val())){
				return *m_cache;
			}else{
				if(Wnorm<=0)throw Exception<WeightedAverageCalculator>("Attempt to check empty data");
				(*m_cache)=value<numt>(Sum/Wnorm,1.0/sqrt(Wnorm));
				return *m_cache;
			}
		}
	};
};
#endif
