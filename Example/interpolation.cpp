// this file is distributed under 
// MIT license
#include <math_h/interpolate.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	LinearInterpolation<double> test{
		point<double>(0,0),point<double>(1,1),point<double>(2,4),point<double>(3,9)
	};
	test<<point<double>(-1,1)<<point<double>(-2,4)<<point<double>(-3,9);
	
	auto x_chain_for_displaying=ChainWithStep(-3.0,0.1,3.0);
	Plotter::Instance().SetOutput(".","interpolation");
	Plot<double> output;
	output.Line(SortedPoints<double>(test.func(),x_chain_for_displaying),"y (interpolated)");
	output.Line(test*2.0,"y*2")<<"set key on";
	
	BiLinearInterpolation<double> test2({0.0,1.0,2.0,3.0,4.0},{0.0,1.0,2.0,3.0,4.0});
	test2.FullCycleVar([](const double&x,const double&y,double&z){
		z=sin(sqrt(pow(x,2)+pow(y,2))/2.0);
	});
	BiSortedPoints<double> for_plot(ChainWithStep(0.0,0.1,4.0),ChainWithStep(0.0,0.1,4.0));
	for_plot.FullCycleVar([&test2](const double&x,const double&y,double&z){
		z=test2(x,y);
	});
	PlotHist2d<double>(normal).Surface(for_plot);
	return 0;
}