// this file is distributed under 
// MIT license
#ifndef ____STRUCTURES_H_____
#	define ____STRUCTURES_H_____
#include <vector>
#include <functional>
#include "sigma.h"
#include "error.h"
namespace MathTemplates{
	using namespace std;
	template<class numtX,class numtY=numtX>class point{
	private:
		numtX x;
		numtY y;
	public:
		point(const numtX&pos):x(pos){}
		point(const numtX&&pos):x(pos){}
		point(const numtX&pos,const numtY&val):x(pos),y(val){}
		point(const numtX&&pos,const numtY&val):x(pos),y(val){}
		point(const numtX&pos,const numtY&&val):x(pos),y(val){}
		point(const numtX&&pos,const numtY&&val):x(pos),y(val){}
		point(const point&source):x(source.x),y(source.y){}
		const numtX&X()const{return x;}
		const numtY&Y()const{return y;}
	protected:
		numtX&__X(){return x;}
		numtY&__Y(){return y;}
	};
	template<class numtX,class numtY=numtX> class point_editable_x:public point<numtX,numtY>{
	public:
		point_editable_x(const point<numtX,numtY>&source):point<numtX,numtY>(source){}
		point_editable_x(const point<numtX,numtY>&&source):point<numtX,numtY>(source){}
		numtX&varX(){return point<numtX,numtY>::__X();}
	};
	template<class numtX,class numtY=numtX> class point_editable_y:public point<numtX,numtY>{
	public:
		point_editable_y(const point<numtX,numtY>&source):point<numtX,numtY>(source){}
		point_editable_y(const point<numtX,numtY>&&source):point<numtX,numtY>(source){}
		numtY&varY(){return point<numtX,numtY>::__Y();}
	};
	template<class numtX,class numtY=numtX> class point_editable:public point<numtX,numtY>{
	public:
		point_editable(const point<numtX,numtY>&source):point<numtX,numtY>(source){}
		point_editable(const point<numtX,numtY>&&source):point<numtX,numtY>(source){}
		numtX&varX(){return point<numtX,numtY>::__X();}
		numtY&varY(){return point<numtX,numtY>::__Y();}
	};
	template<class numX, class numY=numX>
	bool operator<(const point<numX,numY>&a,const point<numX,numY>&b){
		return a.X()<b.X();
	}
	template<class numX, class numY=numX>
	bool operator>(const point<numX,numY>&a,const point<numX,numY>&b){
		return a.X()>b.X();
	}
	
	namespace details{
		template<class comparable, class indexer=vector<comparable>>
		int  WhereToInsert(const int from,const int to,const indexer&X,const comparable&x){
			if(from>to) return from;
			int beg=from,end=to;
			if(x>X[end]) return end+1;
			if(x<X[beg]) return beg;
			while(1<(end-beg)){
				int mid=(beg+end)/2;
				if(x<X[mid]) end=mid;
				else
					if(x>X[mid]) beg=mid;
					else return mid;
			}
			return end;
		}
		template<class comparable, class indexer=vector<comparable>>
		int Search(const int from,const int to,const indexer&X,const comparable&x){
			if(from>to) return from-1;
			int beg=from,end=to;
			if(x>X[end]) return from-1;
			if(x<X[beg]) return from-1;
			while(1<(end-beg)){
				int mid=(beg+end)/2;
				if(x<X[mid]) end=mid;
				else
					if(x>X[mid]) beg=mid;
					else return mid;
			}
			if(X[beg]==x)return beg;
			if(X[end]==x)return end;
			return from-1;
		}
	}
	template<class comparable,class indexer, class Size, class Insert>
	void InsertSorted(const comparable&x,indexer&X,const Size size,const Insert insert){
		insert(details::WhereToInsert(0,size()-1,X,x),x);
	}
	template<class comparable,class indexer, class Size, class Insert>
	void InsertSorted(const comparable&&x,indexer&X,const Size size,const Insert insert){
		InsertSorted(x,X,size,insert);
	}
	#define std_size(vector) [&vector](){return (vector).size();}
	#define std_insert(vector,type) [&vector](int pos,type x){(vector).insert((vector).begin()+pos,x);}
	#define field_size(vector)  [this](){return (vector).size();}
	#define field_insert(vector,type)  [this](int pos,type x){(vector).insert((vector).begin()+pos,x);}
	
	template<class numX>
	const vector<numX> ChainWithStep(const numX from,const numX step,const numX to){
		if(from>=to)throw Exception<vector<numX>>("wrong binning ranges");
		if(step<=0)throw Exception<vector<numX>>("wrong binning step");
		vector<numX> res;
		for(numX x=from;x<=to;x+=step)res.push_back(x);
		return res;
	}
	template<class numX>
	const vector<numX> ChainWithCount(const size_t cont,const numX from,const numX to){
		if(from>=to)throw Exception<vector<numX>>("wrong binning ranges");
		if(0==cont)throw Exception<vector<numX>>("wrong bins count");
		numX step=(to-from)/numX(cont);
		vector<numX> res;
		for(numX x=from;x<=to;x+=step)res.push_back(x);
		return res;
	}
	
	
	template<class numX, class numY=numX>
	class SortedPoints{
	private:
		vector<point_editable_y<numX,numY>> data;
	public:
		SortedPoints(){}
		SortedPoints&operator<<(const point<numX,numY>&p){
			typedef point_editable_y<numX,numY> P;
			InsertSorted(P(p),data,field_size(data),field_insert(data,P));
			return *this;
		}
		SortedPoints&operator<<(const point<numX,numY>&&p){
			return operator<<(p);
		}
		SortedPoints(const initializer_list<point<numX,numY>>&points){
			for(const auto&p:points)operator<<(p);
		}
		SortedPoints(const initializer_list<point<numX,numY>>&&points):SortedPoints(points){}
		SortedPoints(const SortedPoints&points){
			for(const auto&p:points.data)operator<<(p);
		}
		SortedPoints(const SortedPoints&&points):SortedPoints(points){}
		SortedPoints(const function<numY(numX)> f,const vector<numX>&chain){
			for(numX x:chain)operator<<(point<numX,numY>(x,f(x)));
		}
		SortedPoints(const function<numY(numX)> f,const vector<numX>&&chain):SortedPoints(f,chain){}
		SortedPoints(const function<numY(numX)> f,const initializer_list<numX>&&chain){
			for(numX x:chain)operator<<(point<numX,numY>(x,f(x)));
		}
		virtual ~SortedPoints(){}

		int size()const{return data.size();}
		const point<numX,numY>&operator[](const int i)const{
			if(size()<=i)
				throw Exception<SortedPoints>("Range check error");
			return data[i];
		}
		typedef typename vector<point_editable_y<numX,numY>>::const_iterator const_iterator;
		const_iterator begin()const{return data.begin();}
		const_iterator cbegin()const{return data.cbegin();}
		const_iterator end() const{return data.end();}
		const_iterator cend() const{return data.cend();}
		const point<numX,numY>&left()const{
			if(size()<1)
				throw Exception<SortedPoints>("Attempt to obtain empty properties.");
			return data[0];
		}
		const point<numX,numY>&right()const{
			if(size()<1)
				throw Exception<SortedPoints>("Attempt to obtain empty properties.");
			return data[size()-1];
		}
		numX min()const{return left().X();}
		numX max()const{return right().X();}
		//Arithmetic actions
		const SortedPoints<numY,numX> Transponate()const{
			SortedPoints<numY,numX> res;
			for(const point_editable_y<numX,numY>&p:data)
				res<<point<numX,numY>(p.Y(),p.X());
			return res;
		}
		SortedPoints &operator+=(const numY value){
			for(point_editable_y<numX,numY>&point:data)
				point.varY()+=value;
			return *this;
		}
		SortedPoints &operator-=(const numY value){
			for(point_editable_y<numX,numY>&point:data)
				point.varY()-=value;
			return *this;
		}
		SortedPoints &operator*=(const numY value){
			for(point_editable_y<numX,numY>&point:data)
				point.varY()*=value;
			return *this;
		}
		SortedPoints &operator/=(const numY value){
			for(point_editable_y<numX,numY>&point:data)
				point.varY()/=value;
			return *this;
		}
		SortedPoints &transform(const std::function<numY(numY)>F){
			for(point_editable_y<numX,numY>&point:data)
				point.varY()=F(point.Y());
			return *this;
		}
		SortedPoints &transform(const std::function<numY(numX,numY)>F){
			for(point_editable_y<numX,numY>&point:data)
				point.varY()=F(point.X(),point.Y());
			return *this;
		}
	};
	
	template<class numtX,class numtY=numtX>class hist{
	public:
		typedef point_editable_y<value<numtX>,value<numtY>> Point;
	private:
		vector<Point> m_data;
	public:
		//Copying
		hist(){}
		hist(const initializer_list<value<numtX>>&data){
			for(const auto& v:data)m_data.push_back(Point(v));
		}
		hist(const initializer_list<value<numtX>>&&data):hist(data){}
		hist(const vector<value<numtX>>&data){
			for(const auto& v:data)m_data.push_back(Point(v));
		}
		hist(const vector<value<numtX>>&&data):hist(data){}
		hist(const initializer_list<point<numtX,numtY>>&data){
			for(const auto& P:data)m_data.push_back(Point(P));
		}
		hist(const initializer_list<point<numtX,numtY>>&&data):hist(data){}
		hist(const vector<point<numtX,numtY>>&data){
			for(const auto& P:data)m_data.push_back(Point(P));
		}
		hist(const vector<point<numtX,numtY>>&&data):hist(data){}
		hist(const initializer_list<Point>&data){
			for(const auto& P:data)m_data.push_back(P);
		}
		hist(const initializer_list<Point>&&data):hist(data){}
		hist(const vector<Point>&data){
			for(const auto& P:data)m_data.push_back(P);
		}
		hist(const vector<Point>&&data):hist(data){}
		hist(const hist&source){
			for(const Point& P:source.m_data)m_data.push_back(P);
		}
		virtual ~hist(){}
		hist&operator=(const hist& source){
			m_data.clear();
			for(const Point& P:source.m_data)m_data.push_back(P);
			return *this;
		}
		//Iterating
		typedef typename vector<Point>::const_iterator const_iterator;
		const_iterator begin()const{return m_data.cbegin();}
		const_iterator cbegin()const{return m_data.cbegin();}
		const_iterator end() const{return m_data.cend();}
		const_iterator cend() const{return m_data.cend();}
		size_t size()const{return m_data.size();}
		const Point&operator[](const size_t i)const{
			if(m_data.size()<=i)
				throw Exception<hist>("range check error");
			return m_data[i];
		}
		Point&Bin(const size_t i){
			if(m_data.size()<=i)
				throw Exception<hist>("range check error");
			return m_data[i];
		}
		//Simple transform
		hist CloneEmptyBins()const{
			vector<Point> initer;
			for(const Point&P:m_data)initer.push_back(P);
			return hist(initer);
		}
		numtY Total()const{
			numtY res=0;
			for(const Point&P:m_data)res+=P.Y().val();
			return res;
		}
		value<numtY> TotalSum()const{
			value<numtY> res(0,0);
			for(const Point&P:m_data)res+=P.Y();
			return res;
		}
		SortedPoints<numtX,numtY> Line()const{
			SortedPoints<numtX,numtY> res;
			for(const Point&P:m_data)
				res<<make_pair(P.X().val(),P.Y().val());
			return res;
		}
		
		//Arithmetic operations
		hist&operator+=(const hist& second){
			for(size_t i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()+=second[i].Y();
				}else
					throw Exception<hist>("Cannot add histogram. bins differ");
			}
			return *this;
		}
		hist&operator+=(const function<numtY(numtX)>f){
			for(size_t i=0,n=size();i<n;i++)
				m_data[i].varY()+=value<numtY>(f(m_data[i].X().val()),0);
			return *this;
		}
		hist&operator+=(const value<numtY>&c){
			for(size_t i=0,n=size();i<n;i++)m_data[i].varY()+=c;
			return *this;
		}
		hist&operator+=(const value<numtY>&&c){return operator+=(c);}
		
		hist&operator-=(const hist& second){
			for(size_t i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()-=second[i].Y();
				}else
					throw Exception<hist>("Cannot add histogram. bins differ");
			}
			return *this;
		}
		hist&operator-=(const function<numtY(numtX)>f){
			for(size_t i=0,n=size();i<n;i++)
				m_data[i].varY()-=value<numtY>(f(m_data[i].X().val()),0);
			return *this;
		}
		hist&operator-=(const value<numtY>&c){
			for(size_t i=0,n=size();i<n;i++)m_data[i].varY()-=c;
			return *this;
		}
		hist&operator-=(const value<numtY>&&c){return operator-=(c);}
		
		hist&operator*=(const hist& second){
			for(size_t i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()*=second[i].Y();
				}else
					throw Exception<hist>("Cannot add histogram. bins differ");
			}
			return *this;
		}
		hist&operator*=(const function<numtY(numtX)>f){
			for(size_t i=0,n=size();i<n;i++)
				m_data[i].varY()*=value<numtY>(f(m_data[i].X().val()),0);
			return *this;
		}
		hist&operator*=(const value<numtY>&c){
			for(size_t i=0,n=size();i<n;i++)m_data[i].varY()*=c;
			return *this;
		}
		hist&operator*=(const value<numtY>&&c){return operator*=(c);}
		
		hist&operator/=(const hist& second){
			for(size_t i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()*=second[i].Y();
				}else
					throw Exception<hist>("Cannot add histogram. bins differ");
			}
			return *this;
		}
		hist&operator/=(const function<numtY(numtX)>f){
			for(size_t i=0,n=size();i<n;i++)
				m_data[i].varY()/=value<numtY>(f(m_data[i].X().val()),0);
			return *this;
		}
		hist&operator/=(const value<numtY>&c){
			for(size_t i=0,n=size();i<n;i++)m_data[i].varY()/=c;
			return *this;
		}
		hist&operator/=(const value<numtY>&&c){return operator/=(c);}
		
		//Advanced transrormations
		hist&FillWithValues(const value<numtY>&v){
			for(Point&P:m_data)P.varY()=v;
			return *this;
		}
		hist&FillWithValues(const value<numtY>&&v){
			return FillWithValues(v);
		}
		hist Scale(const size_t sc_x)const{
			//uncertanties are set to standard sqrt
			vector<value<numtX>> new_x,sorted_x;
			for(const auto&item:m_data)
				InsertSorted(item.X(),sorted_x,std_size(sorted_x),std_insert(sorted_x,value<numtX>));
			for(size_t i=sc_x-1,n=sorted_x.size();i<n;i+=sc_x){
				auto min=sorted_x[i+1-sc_x].min();
				auto max=sorted_x[i].max();
				new_x.push_back(value<numtX>((max+min)/numtX(2),(max-min)/numtX(2)));
			}
			hist res(new_x);
			for(size_t i=0;i<new_x.size();i++){
				numtY v=0;
				for(size_t ii=0;ii<sc_x;ii++)
					v+=m_data[i*sc_x+ii].Y().val();
				res.Bin(i).varY()=value<numtY>(v,sqrt(v));
			}
			return res;
		}
		hist&imbibe(const hist& second){
			//Sum of histograms. the uncertanties are set to standard sqrt
			for(int i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()=value<numtY>(m_data[i].Y().val()+second[i].Y().val());
				}else
					throw Exception<hist>("Cannot imbibe histogram. bins differ");
			}
			return *this;
		}
	};
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,const hist<numtX,numtY>&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,const hist<numtX,numtY>&&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&&b){auto res=a;res+=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,const function<numtY(numtX)>b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&&a,const function<numtY(numtX)>b){auto res=a;res+=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,const value<numtY>&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&&a,const value<numtY>&b){return a+b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,const value<double>&&b){return a+b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&&a,const value<double>&&b){return a+b;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,const hist<numtX,numtY>&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,const hist<numtX,numtY>&&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&&b){auto res=a;res-=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,const function<numtY(numtX)>b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&&a,const function<numtY(numtX)>b){auto res=a;res-=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,const value<numtY>&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&&a,const value<numtY>&b){return a-b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,const value<numtY>&&b){return a-b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&&a,const value<numtY>&&b){return a-b;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,const hist<numtX,numtY>&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,const hist<numtX,numtY>&&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&&b){auto res=a;res*=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,const function<numtY(numtX)>b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&&a,const function<numtY(numtX)>b){auto res=a;res*=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,const value<numtY>&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&&a,const value<numtY>&b){return a*b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,const value<numtY>&&b){return a*b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&&a,const value<numtY>&&b){return a*b;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,const hist<numtX,numtY>&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,const hist<numtX,numtY>&&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&&b){auto res=a;res/=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,const function<numtY(numtX)>b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&&a,const function<numtY(numtX)>b){auto res=a;res/=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,const value<numtY>&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&&a,const value<numtY>&b){return a/b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,const value<numtY>&&b){return a/b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&&a,const value<numtY>&&b){return a/b;}
	
	
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class point3d{
	private:
		numtX x;
		numtY y;
		numtZ z;
	public: 
		point3d(const numtX&_x,const numtY&_y,const numtZ&_z)
		:x(_x),y(_y),z(_z){}
		point3d(const numtX&&_x,const numtY&&_y,const numtZ&&_z)
		:x(_x),y(_y),z(_z){}
		virtual ~point3d(){}
		const numtX&X()const{return x;}
		const numtY&Y()const{return y;}
		const numtZ&Z()const{return z;}
	};
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class hist2d{
	private:
		vector<value<numtX>> m_x_axis;
		vector<value<numtY>> m_y_axis;
		vector<vector<value<numtZ>>> m_data;
		void init(){
			m_data.clear();
			for(size_t i=0,I=m_x_axis.size();i<I;i++){
				m_data.push_back(vector<value<numtZ>>());
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					m_data[m_data.size()-1].push_back(value<numtZ>(0));
			}
		}
	public:
		hist2d(const initializer_list<value<numtX>>&X,const initializer_list<value<numtY>>&Y){
			for(const auto&x:X)m_x_axis.push_back(x);
			for(const auto&y:Y)m_y_axis.push_back(y);
			init();
		}
		hist2d(const initializer_list<value<numtX>>&&X,const initializer_list<value<numtY>>&&Y):hist2d(X,Y){}
		hist2d(const vector<value<numtX>>&X,const vector<value<numtY>>&Y){
			for(const auto&x:X)m_x_axis.push_back(x);
			for(const auto&y:Y)m_y_axis.push_back(y);
			init();
		}
		hist2d(const vector<value<numtX>>&&X,const vector<value<numtY>>&&Y):hist2d(X,Y){}
		hist2d():hist2d({},{}){}
		hist2d(const hist2d&source){
			for(const auto&item:source.m_x_axis)m_x_axis.push_back(item);
			for(const auto&item:source.m_y_axis)m_y_axis.push_back(item);
			for(size_t i=0,I=source.m_data.size();i<I;i++){
				m_data.push_back(vector<value<numtZ>>());
				for(const auto&item:source.m_data[i])
					m_data[i].push_back(item);
			}
		}
		virtual ~hist2d(){}
		typedef typename vector<vector<value<numtZ>>>::const_iterator const_iterator;
		const_iterator begin()const{return m_data.cbegin();}
		const_iterator cbegin()const{return m_data.cbegin();}
		const_iterator end() const{return m_data.cend();}
		const_iterator cend() const{return m_data.cend();}
		size_t size()const{return m_data.size();}
		const vector<value<numtZ>>&operator[](const size_t i)const{
			if(size()<=i)throw Exception<hist2d>("range check error");
			return m_data[i];
		}
		value<numtZ>&Bin(const size_t i,const size_t j){
			if(size()<=i)throw Exception<hist2d>("range check error");
			if(m_data[i].size()<=j)throw Exception<hist2d>("range check error");
			return m_data[i][j];
		}
		const vector<value<numtX>>&X()const{return m_x_axis;}
		const vector<value<numtY>>&Y()const{return m_y_axis;}
		
		hist2d Scale(const size_t sc_x,const size_t sc_y)const{
			//uncertanties are set to standard sqrt
			vector<value<numtX>> new_x,sorted_x;
			for(const auto&item:X())
				InsertSorted(item,sorted_x,std_size(sorted_x),std_insert(sorted_x,value<numtX>));
			for(size_t i=sc_x-1,n=sorted_x.size();i<n;i+=sc_x){
				auto min=sorted_x[i+1-sc_x].min();
				auto max=sorted_x[i].max();
				new_x.push_back(value<numtX>((max+min)/numtX(2),(max-min)/numtX(2)));
			}
			vector<value<numtY>> new_y,sorted_y;
			for(const auto&item:Y())
				InsertSorted(item,sorted_y,std_size(sorted_y),std_insert(sorted_y,value<numtY>));
			for(size_t i=sc_y-1,n=sorted_y.size();i<n;i+=sc_y){
				auto min=sorted_y[i+1-sc_y].min();
				auto max=sorted_y[i].max();
				new_y.push_back(value<numtY>((max+min)/numtY(2),(max-min)/numtY(2)));
			}
			hist2d res(new_x,new_y);
			for(size_t i=0;i<new_x.size();i++)for(size_t j=0;j<new_y.size();j++){
				numtZ v=0;
				for(size_t ii=0;ii<sc_x;ii++)
					for(size_t jj=0;jj<sc_y;jj++)
						v+=m_data[i*sc_x+ii][j*sc_y+jj].val();
					res.Bin(i,j)=value<numtZ>(v,sqrt(v));
			}
			return res;
		}
		hist2d&imbibe(const hist2d& second){
			//Sum of histograms. the uncertanties are set to standard sqrt
			if((X().size()!=second.X().size())||(Y().size()!=second.Y().size()))
				throw Exception<hist2d>("cannot imbibe second histogram: bins differ");
			for(int i=0,n=size();i<n;i++)for(int j=0,m=m_data[i].size();j<m;j++)
				m_data[i][j]=value<numtY>(m_data[i][j].val()+second[i][j].val());
			return *this;
		}
		
		const point3d<value<numtX>,value<numtY>,value<numtZ>> operator()(const size_t i,const size_t j)const{
			if(size()<=i)throw Exception<hist2d>("range check error");
			if(m_y_axis.size()<=j)throw Exception<hist2d>("range check error");
			return point3d<value<numtX>,value<numtY>,value<numtZ>>(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
		void FullCycle(const function<void(const point3d<value<numtX>,value<numtY>,value<numtZ>>&)>f)const{
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++){
					point3d<value<numtX>,value<numtY>,value<numtZ>> P(m_x_axis[i],m_y_axis[j],m_data[i][j]);
					f(P);
				}
		}
		void FullCycle(const function<void(const value<numtX>&,const value<numtY>&,const value<numtZ>&)>f)const{
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					f(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
		void FullCycleVar(const function<void(const value<numtX>&,const value<numtY>&,value<numtZ>&)>f){
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					f(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
	};
	
	
	template<class numtX,class numtY=numtX>
	class Distribution1D:public hist<numtX,numtY>{
	public:
		Distribution1D(const initializer_list<value<numtX>>&data):hist<numtX,numtY>(data){hist<numtX,numtY>::FillWithValues(value<numtY>(0));}
		Distribution1D(const initializer_list<value<numtX>>&&data):hist<numtX,numtY>(data){hist<numtX,numtY>::FillWithValues(value<numtY>(0));}
		Distribution1D(const vector<value<numtX>>&data):hist<numtX,numtY>(data){hist<numtX,numtY>::FillWithValues(value<numtY>(0));}
		Distribution1D(const vector<value<numtX>>&&data):hist<numtX,numtY>(data){hist<numtX,numtY>::FillWithValues(value<numtY>(0));}
		Distribution1D&operator<<(const numtX v){
			for(size_t i=0,n=hist<double>::size();i<n;i++)
				if(hist<double>::Bin(i).X().contains(v))
					hist<double>::Bin(i).varY()=value<numtY>(hist<double>::Bin(i).Y().val()+numtY(1));
				return *this;
		}
	};
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class Distribution2D:public hist2d<numtX,numtY,numtZ>{
	public:
		Distribution2D(const initializer_list<value<numtX>>&X,const initializer_list<value<numtY>>&Y):hist2d<numtX,numtY,numtZ>(X,Y){}
		Distribution2D(const initializer_list<value<numtX>>&&X,initializer_list<value<numtY>>&&Y):hist2d<numtX,numtY,numtZ>(X,Y){}
		Distribution2D(const vector<value<numtX>>&X,const vector<value<numtY>>&Y):hist2d<numtX,numtY,numtZ>(X,Y){}
		Distribution2D(const vector<value<numtX>>&&X,vector<value<numtY>>&&Y):hist2d<numtX,numtY,numtZ>(X,Y){}
		Distribution2D&operator<<(const pair<numtX,numtY>&p){
			for(size_t i=0,I=hist2d<numtX,numtY,numtZ>::X().size();i<I;i++)if(hist2d<numtX,numtY,numtZ>::X()[i].contains(p.first))
				for(size_t j=0,J=hist2d<numtX,numtY,numtZ>::Y().size();j<J;j++)if(hist2d<numtX,numtY,numtZ>::Y()[j].contains(p.second))
					hist2d<numtX,numtY,numtZ>::Bin(i,j)=value<numtZ>(hist2d<numtX,numtY,numtZ>::operator[](i)[j].val()+numtZ(1));
				return *this;
		}
		Distribution2D&operator<<(const pair<numtX,numtY>&&p){
			return operator<<(p);
		}
	};
	
};

#endif