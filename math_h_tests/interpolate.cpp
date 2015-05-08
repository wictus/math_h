#include <gtest/gtest.h>
#include <vector>
#include <math.h>
#include <interpolate.h>
using namespace std;
double TestArray[]={-2,-1,0,1,2,3,4,5,6,7};
TEST(WhereToInsert,BorderConditions){
	for(int beg=0;beg<10;beg++)
		for(double V=TestArray[beg]-2;V<=TestArray[beg]+2;V+=0.5){
			int index=WhereToInsert(beg,beg-1,TestArray,V);
			EXPECT_EQ(beg,index);
			index=WhereToInsert(beg,beg,TestArray,V);
			if(V<TestArray[beg])EXPECT_EQ(beg,index);
			if(V>TestArray[beg])EXPECT_EQ(beg+1,index);
			if(V==TestArray[beg])ASSERT_TRUE((index==beg)||(index==(beg+1)));
		}
}
TEST(WhereToInsert,NormalConditions){
	for(double x=TestArray[0]-0.5;x<=TestArray[9]+0.5;x+=0.5)
		for(int beg=0;beg<10;beg++)for(int end=beg+1;end<10;end++){
			int index=WhereToInsert(beg,end,TestArray,x);
			if(x<TestArray[beg])EXPECT_EQ(beg,index);
			else if(x>TestArray[end])EXPECT_EQ(end+1,index);
			else ASSERT_TRUE((x>=TestArray[index-1])&&(x<=TestArray[index]));
		}
}
TEST(InsertSorted,BasicTest){
	vector<int> X;
	for(int i=0;i<50;i++){
		InsertSorted(rand()%10,X,std_size(X),std_insert(X,int));
		for(int j=0;j<i;j++)
			ASSERT_TRUE(X[j]<=X[j+1]);
	}
}
TEST(LinearInterpolation,Create){
	LinearInterpolation<double> F;
	EXPECT_EQ(0,F.size());
	ASSERT_ANY_THROW(F.min());
	ASSERT_ANY_THROW(F.max());
	EXPECT_EQ(&F,&(F<<make_pair(0,0)));
	EXPECT_EQ(1,F.size());
	ASSERT_NO_THROW(F.min());
	ASSERT_NO_THROW(F.max());
	ASSERT_ANY_THROW(F(0.5));
	EXPECT_EQ(&F,&(F<<make_pair(1,0)));
	EXPECT_EQ(2,F.size());
	ASSERT_NO_THROW(F.min());
	ASSERT_NO_THROW(F.max());
	EXPECT_EQ(0,F(0.5));
}
#define _EQ(a,b) EXPECT_TRUE(pow(a-b,2)<0.0001)
TEST(LinearInterpolation,SimpleLine){
	LinearInterpolation<double> F;
	F<<make_pair(0,0)<<make_pair(1,1);
	EXPECT_EQ(0,F(0));
	EXPECT_EQ(1,F(1));
	_EQ(0.1,F(0.1));
	_EQ(0.2,F(0.2));
	_EQ(0.3,F(0.3));
	_EQ(0.4,F(0.4));
	_EQ(0.5,F(0.5));
	_EQ(0.6,F(0.6));
	_EQ(0.7,F(0.7));
	_EQ(0.8,F(0.8));
	_EQ(0.9,F(0.9));
	F<<make_pair(2,0);
	EXPECT_EQ(0,F(0));
	EXPECT_EQ(1,F(1));
	_EQ(0.1,F(0.1));
	_EQ(0.2,F(0.2));
	_EQ(0.3,F(0.3));
	_EQ(0.4,F(0.4));
	_EQ(0.5,F(0.5));
	_EQ(0.6,F(0.6));
	_EQ(0.7,F(0.7));
	_EQ(0.8,F(0.8));
	_EQ(0.9,F(0.9));
	EXPECT_EQ(0,F(2));
	_EQ(0.1,F(1.9));
	_EQ(0.2,F(1.8));
	_EQ(0.3,F(1.7));
	_EQ(0.4,F(1.6));
	_EQ(0.5,F(1.5));
	_EQ(0.6,F(1.4));
	_EQ(0.7,F(1.3));
	_EQ(0.8,F(1.2));
	_EQ(0.9,F(1.1));
	F<<make_pair(0.5,0.5)<<make_pair(1.5,0.5);
	EXPECT_EQ(0,F(0));
	EXPECT_EQ(1,F(1));
	_EQ(0.1,F(0.1));
	_EQ(0.2,F(0.2));
	_EQ(0.3,F(0.3));
	_EQ(0.4,F(0.4));
	EXPECT_EQ(0.5,F(0.5));
	_EQ(0.6,F(0.6));
	_EQ(0.7,F(0.7));
	_EQ(0.8,F(0.8));
	_EQ(0.9,F(0.9));
	EXPECT_EQ(0,F(2));
	_EQ(0.1,F(1.9));
	_EQ(0.2,F(1.8));
	_EQ(0.3,F(1.7));
	_EQ(0.4,F(1.6));
	EXPECT_EQ(0.5,F(1.5));
	_EQ(0.6,F(1.4));
	_EQ(0.7,F(1.3));
	_EQ(0.8,F(1.2));
	_EQ(0.9,F(1.1));
}
