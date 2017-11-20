#include <iostream>
#include <math_h/vectors.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
//This is an example how to use vectors.h that provides
//a simple Monte Carlo simulation of particle scattering
//The vectors must be given in unit system that assumes c=1
//and equal units for energy, momentum and mass
int main()
{
    RANDOM RG;
    //Particles masses
    const double M1 = 0.75, M2 = 1.0;
    //momentum distribution
    const RandomGauss<> P(0.3,0.01);

    PlotDistr2D<>
    P_vs_P("P-P", BinsByCount(100, 0.0, 0.5), BinsByCount(100, 0.0, 0.5), "plot_PP"),
    E_vs_E("Ek-Ek", BinsByCount(100, 0.0, 0.1), BinsByCount(100, 0.0, 0.1), "plot_EE"),
    Th_vs_Th("Th-Th", BinsByCount(90, 0.0, 180.), BinsByCount(90, 0.0, 90.), "plot_ThTh");
    for (size_t i = 0; i < 100000; i++) {
	//Create Lorentz vectors for input particles
        const auto Pr0 = lorentz_byPM(Z() * P(RG), M1);
        const auto Pt0 = lorentz_byPM(Zero(), M2);
        const auto Total = Pr0 + Pt0;
        //Let them scatter isotropically in CM
        const auto finalCM = binaryDecay(Total.M(), M1, M2, randomIsotropic<3>(RG));
        //Converting to lab system.
        const auto final1 = finalCM.first.Transform(-Total.Beta());
        const auto final2 = finalCM.second.Transform(-Total.Beta());
        //Gather kinematic statistics
        Th_vs_Th.Fill(direction(final1.P()).th() * 180. / PI(), direction(final2.P()).th() * 180. / PI());
        P_vs_P.Fill(final1.P().M(), final2.P().M());
        E_vs_E.Fill(final1.Ekin(), final2.Ekin());
    }
}

