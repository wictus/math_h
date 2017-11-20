#include <math_h/interpolate.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
//This is an example how to use
// linear interpolation class
int main()
{
    //create the table of points for interpolation
    LinearInterpolation<> test=Points<>{{0, 0}, {1, 1}, {2, 4}, {3, 9}};
    //you can add more points later
    test << make_point(-1, 1) << make_point(-2, 4) << make_point(-3, 9);
    //create chain of values with lesser step
    SortedPoints<> plot_data;
    for(double x=-3.;x<=3.;x+=0.1){
	const auto y=test(x);//that's how we interpolate
	plot_data<<make_point(x,y);
    }
    Plot("interpolation1").Points(plot_data);
    return 0;
}
