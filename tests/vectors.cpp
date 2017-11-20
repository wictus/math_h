// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/vectors.h>
#include <math_h/hists.h>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
const double epsilon = 0.0000000001;
TEST(Vector, scalar_prod2)
{
    EXPECT_EQ(0, x()*y());
    EXPECT_EQ(0, y()*x());
    EXPECT_EQ(1, x()*x());
    EXPECT_EQ(1, y()*y());

    EXPECT_EQ(0, x() * (-y()));
    EXPECT_EQ(0, y() * (-x()));
    EXPECT_EQ(0, (-x())*y());
    EXPECT_EQ(0, (-y())*x());

    EXPECT_EQ(-1, x() * (-x()));
    EXPECT_EQ(-1, y() * (-y()));
    EXPECT_EQ(-1, (-x())*x());
    EXPECT_EQ(-1, (-y())*y());
}
TEST(Vector, scalar_prod3)
{
    EXPECT_EQ(0, X()*Y());
    EXPECT_EQ(0, Y()*Z());
    EXPECT_EQ(0, X()*Z());
    EXPECT_EQ(1, X()*X());
    EXPECT_EQ(1, Y()*Y());
    EXPECT_EQ(1, Z()*Z());

    EXPECT_EQ(0, X() * (-Y()));
    EXPECT_EQ(0, (-X())*Z());

    EXPECT_EQ(-1, (-X())*X());
    EXPECT_EQ(-1, Y() * (-Y()));
    EXPECT_EQ(1, (-Z()) * (-Z()));
}
TEST(Vector, pseudoscalar_prod)
{
    EXPECT_EQ(x()^y(), desCartes(1.));
    EXPECT_EQ(y()^x(), desCartes(-1.));
    EXPECT_EQ(x()^x(), desCartes(0.));
    EXPECT_EQ(y()^y(), desCartes(0.));
    EXPECT_EQ((-x())^x(), desCartes(0.));
    EXPECT_EQ((-y())^y(), desCartes(0.));
}
TEST(Vector, pseudoscalar_prod_algebra)
{
    RANDOM RG;
    RandomUniform<> M(-50, 50);
    for (size_t i = 0; i < 50; i++) {
        const auto
        A = randomIsotropic<2>(RG) * M(RG),
        B = randomIsotropic<2>(RG) * M(RG),
        C = randomIsotropic<2>(RG) * M(RG);
        const auto a = M(RG);
        EXPECT_TRUE((A ^ A).M() < epsilon);
        EXPECT_TRUE((B ^ B).M() < epsilon);
        EXPECT_TRUE(((A ^ B) + (B ^ A)).M() < epsilon);
        EXPECT_TRUE(((A ^ B)*a - ((A * a)^B)).M() < epsilon);
        EXPECT_TRUE(((A ^ B)*a - (A ^ (B * a))).M() < epsilon);
        EXPECT_TRUE((((A + B)^C) - ((A ^ C) + (B ^ C))).M() < epsilon);
    }
}
TEST(Vector, vector_prod3_basis)
{
    EXPECT_EQ(X()^Y(), Z());
    EXPECT_EQ(X()^X(), Zero());
    EXPECT_EQ(Y()^X(), -Z());
    EXPECT_EQ(X() ^ (-X()), Zero());

    EXPECT_EQ(1, (X()^Y())*Z());
    EXPECT_EQ(1, (Y()^Z())*X());
    EXPECT_EQ(1, (Z()^X())*Y());

}
TEST(Vector, vector_prod3_algebra)
{
    RANDOM RG;
    RandomUniform<> M(-50, 50);
    for (size_t i = 0; i < 50; i++) {
        const auto
        A = randomIsotropic<3>(RG) * M(RG),
        B = randomIsotropic<3>(RG) * M(RG),
        C = randomIsotropic<3>(RG) * M(RG);
        const auto a = M(RG);
        EXPECT_TRUE((A ^ A).M() < epsilon);
        EXPECT_TRUE((B ^ B).M() < epsilon);
        EXPECT_TRUE((C ^ C).M() < epsilon);
        EXPECT_TRUE(((A ^ B) + (B ^ A)).M() < epsilon);
        EXPECT_TRUE(((A ^ C) + (C ^ A)).M() < epsilon);
        EXPECT_TRUE(((C ^ B) + (B ^ C)).M() < epsilon);
        EXPECT_TRUE(((A ^ B)*a - ((A * a)^B)).M() < epsilon);
        EXPECT_TRUE(((A ^ B)*a - (A ^ (B * a))).M() < epsilon);
        EXPECT_TRUE((((A + B)^C) - ((A ^ C) + (B ^ C))).M() < epsilon);
        EXPECT_TRUE((((A ^ B)^C) + ((C ^ A)^B) + ((B ^ C)^A)).M() < epsilon);
        EXPECT_TRUE(((A ^ (B ^ C)) - (B * (A * C)) + (C * (A * B))).M() < epsilon);
        EXPECT_TRUE(abs(((A ^ B)*C) - (A * (B ^ C))) < epsilon);
    }
}
TEST(Vector, vector_prod7_algebra)
{
    RANDOM RG;
    RandomUniform<> M(-50, 50);
    for (size_t i = 0; i < 50; i++) {
        const auto
        A = randomIsotropic<7>(RG) * M(RG),
        B = randomIsotropic<7>(RG) * M(RG),
        C = randomIsotropic<7>(RG) * M(RG);
        const auto a = M(RG);
        EXPECT_TRUE((A ^ A).M() < epsilon);
        EXPECT_TRUE((B ^ B).M() < epsilon);
        EXPECT_TRUE((C ^ C).M() < epsilon);
        EXPECT_TRUE(((A ^ B) + (B ^ A)).M() < epsilon);
        EXPECT_TRUE(((A ^ C) + (C ^ A)).M() < epsilon);
        EXPECT_TRUE(((C ^ B) + (B ^ C)).M() < epsilon);
        EXPECT_TRUE(((A ^ B)*a - ((A * a)^B)).M() < epsilon);
        EXPECT_TRUE(((A ^ B)*a - (A ^ (B * a))).M() < epsilon);
        EXPECT_TRUE((((A + B)^C) - ((A ^ C) + (B ^ C))).M() < epsilon);
        EXPECT_TRUE(abs(((A ^ B)*C) - (A * (B ^ C))) < epsilon);
    }
}
TEST(Vector, base1d)
{
    RANDOM RG;
    RandomUniform<> X(-50., 50);
    for (size_t i = 0; i < 50; i++) {
        const auto x = X(RG);
        const auto V = desCartes(x);
        EXPECT_EQ(x, V.x());
        EXPECT_EQ(x, V.component<1>());
    }
}
TEST(Vector, base2d)
{
    RANDOM RG;
    RandomUniform<> X(-50., 50);
    for (size_t i = 0; i < 50; i++) {
        const auto x = X(RG), y = X(RG);
        const auto V = desCartes(x, y);
        EXPECT_EQ(x, V.x());
        EXPECT_EQ(x, V.component<1>());
        EXPECT_EQ(y, V.y());
        EXPECT_EQ(y, V.component<2>());
    }
}
TEST(Vector, base3d)
{
    RANDOM RG;
    RandomUniform<> X(-50., 50);
    for (size_t i = 0; i < 50; i++) {
        const auto x = X(RG), y = X(RG), z = X(RG);
        const auto V = desCartes(x, y, z);
        EXPECT_EQ(x, V.x());
        EXPECT_EQ(x, V.component<1>());
        EXPECT_EQ(y, V.y());
        EXPECT_EQ(y, V.component<2>());
        EXPECT_EQ(z, V.z());
        EXPECT_EQ(z, V.component<3>());
    }
}
TEST(Vector, base4d)
{
    RANDOM RG;
    RandomUniform<> X(-50., 50);
    for (size_t i = 0; i < 50; i++) {
        const auto x = X(RG), y = X(RG), z = X(RG), zz = X(RG);
        const auto V = desCartes(x, y, z, zz);
        EXPECT_EQ(x, V.x());
        EXPECT_EQ(x, V.component<1>());
        EXPECT_EQ(y, V.y());
        EXPECT_EQ(y, V.component<2>());
        EXPECT_EQ(z, V.z());
        EXPECT_EQ(z, V.component<3>());
        EXPECT_EQ(zz, V.component<4>());
    }
}

TEST(Direction, base2d)
{
    RANDOM RG;
    RandomUniform<> Phi(-PI(), PI());
    for (size_t i = 0; i < 50; i++) {
        const auto phi = Phi(RG);
        const auto D = direction(phi);
        EXPECT_EQ(D.phi(), phi);
    }
}
TEST(Direction, base3d)
{
    RANDOM RG;
    RandomUniform<> Phi(-PI(), PI());
    RandomUniform<> Theta(0, PI());
    for (size_t i = 0; i < 50; i++) {
        const auto phi = Phi(RG);
        const auto theta = Theta(RG);
        const auto D = direction(phi, theta);
        EXPECT_EQ(D.phi(), phi);
        EXPECT_EQ(D.th(), theta);
    }
}
TEST(Direction, base4d)
{
    RANDOM RG;
    RandomUniform<> Phi(-PI(), PI());
    RandomUniform<> Theta(0, PI<>());
    for (size_t i = 0; i < 50; i++) {
        const auto phi = Phi(RG);
        const auto theta1 = Theta(RG);
        const auto theta2 = Theta(RG);
        const auto D = direction(phi, theta1, theta2);
        EXPECT_EQ(D.phi(), phi);
        EXPECT_EQ(D.th<1>(), theta1);
        EXPECT_EQ(D.th<2>(), theta2);
    }
}
TEST(Direction, base5d)
{
    RANDOM RG;
    RandomUniform<> Phi(-PI(), PI());
    RandomUniform<> Theta(0, PI());
    for (size_t i = 0; i < 50; i++) {
        const auto phi = Phi(RG);
        const auto theta1 = Theta(RG);
        const auto theta2 = Theta(RG);
        const auto theta3 = Theta(RG);
        const auto D = direction(phi, theta1, theta2, theta3);
        EXPECT_EQ(D.phi(), phi);
        EXPECT_EQ(D.th<1>(), theta1);
        EXPECT_EQ(D.th<2>(), theta2);
        EXPECT_EQ(D.th<3>(), theta3);
    }
}

TEST(Vector, Isotropic1)
{
    RANDOM RG;
    const auto c = BinsByCount(2, -2.0, 2.0);
    Distribution1D<> D(c);
    for (size_t i = 0; i < 1000000; i++) {
        D.Fill(randomIsotropic<1>(RG).dir());
    }
    const auto x = D.TotalSum().val() / D.size();
    for (const auto &p : D)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
}
TEST(Vector, Isotropic2)
{
    RANDOM RG;
    const auto c = BinsByCount(10, -PI(), PI());
    Distribution1D<> Phi(c);
    for (size_t i = 0; i < 1000000; i++) {
        Phi.Fill(randomIsotropic<2>(RG).phi());
    }
    const auto x = Phi.TotalSum().val() / Phi.size();
    for (const auto &p : Phi)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
}
TEST(Vector, Isotropic3)
{
    RANDOM RG;
    const auto c = BinsByStep(-1.0, 0.1, 1.0);
    Distribution1D<> X(c), Y(c), Z(c);
    for (size_t i = 0; i < 1000000; i++) {
        const auto V = randomIsotropic<3>(RG) * 1.0;
        X.Fill(V.x());
        Y.Fill(V.y());
        Z.Fill(V.z());
    }
    const auto x = X.TotalSum().val() / X.size();
    for (const auto &p : X)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
    for (const auto &p : Y)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
    for (const auto &p : Z)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
}
TEST(Vector, Rotation2)
{
    RANDOM RG;
    RandomUniform<> Phi(0, PI() * 2.0);
    for (size_t i = 0; i < 50; i++) {
        const auto I = randomIsotropic<2>(RG) * 1.0;
        const auto ang = Phi(RG);
        const auto F = Rotation(ang) * I;
        EXPECT_TRUE(abs(I.M() - F.M()) < epsilon);
    }
}

TEST(Vector, Rotation3)
{
    RANDOM RG;
    RandomUniform<> M(0, 5.0), Phi(-PI(), PI());
    for (size_t i = 0; i < 50; i++) {
        const auto I = randomIsotropic<3>(RG) * M(RG);
        const auto ang = Phi(RG);
        const auto dir = randomIsotropic<3>(RG);
        const auto F = Rotation(dir, ang) * I;
        EXPECT_TRUE(abs((dir * 1.0) * I - (dir * 1.0) * F) < epsilon);
    }
}

TEST(Vector, decompose2)
{
    RANDOM RG;
    RandomUniform<> M(0, 5.0);
    for (size_t i = 0; i < 50; i++) {
        const auto chosen_direction = randomIsotropic<2>(RG);
        const auto vector_to_decompose = randomIsotropic<2>(RG) * M(RG);
        const auto decomposition = decompose_by_direction(vector_to_decompose, chosen_direction);
        EXPECT_TRUE(abs(decomposition.n * (chosen_direction * 1.0)) < epsilon);
        EXPECT_TRUE(vector_to_decompose.CloseTo(decomposition.tau + decomposition.n, epsilon));
    }
}
TEST(Vector, decompose3)
{
    RANDOM RG;
    RandomUniform<> M(0, 5.0);
    for (size_t i = 0; i < 50; i++) {
        const auto chosen_direction = randomIsotropic<3>(RG);
        const auto vector_to_decompose = randomIsotropic<3>(RG) * M(RG);
        const auto decomposition = decompose_by_direction(vector_to_decompose, chosen_direction);
        EXPECT_TRUE(abs(decomposition.n * (chosen_direction * 1.0)) < epsilon);
        EXPECT_TRUE(vector_to_decompose.CloseTo(decomposition.tau + decomposition.n, epsilon));
    }
}
TEST(Vector, decompose4)
{
    RANDOM RG;
    RandomUniform<> M(0, 5.0);
    for (size_t i = 0; i < 50; i++) {
        const auto chosen_direction = randomIsotropic<4>(RG);
        const auto vector_to_decompose = randomIsotropic<4>(RG) * M(RG);
        const auto decomposition = decompose_by_direction(vector_to_decompose, chosen_direction);
        EXPECT_TRUE(abs(decomposition.n * (chosen_direction * 1.0)) < epsilon);
        EXPECT_TRUE(vector_to_decompose.CloseTo(decomposition.tau + decomposition.n, epsilon));
    }
}

TEST(LorentzVector, LorentzTransform1d)
{
    RANDOM RG;
    RandomUniform<> mr(0, 0.99), metrr(-5, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<1>(RG) * 1.0, metrr(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(desCartes(0.));
        EXPECT_TRUE(V1 == V0);
        const auto beta = randomIsotropic<1>(RG) * mr(RG);
        const auto V2 = V0.Transform(beta).Transform(-beta);
        EXPECT_TRUE(abs(V2.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V2.P().CloseTo(V0.P(), epsilon));
        EXPECT_ANY_THROW(V0.Transform(randomIsotropic<1>(RG) * (1.0 + mr(RG))));
        const auto L0 = V0.M(), L1 = V0.Transform(beta).M(),
                   L2 = V0.Transform(beta).Transform(randomIsotropic<1>(RG) * 0.2).M();
        EXPECT_TRUE(abs(L0 - L1) < epsilon);
        EXPECT_TRUE(abs(L2 - L1) < epsilon);
        EXPECT_TRUE(abs(L0 - L2) < epsilon);
        const auto V00 = lorentz_byPM(desCartes(0.), V0.M()).Transform(-V0.Beta());
        EXPECT_TRUE(abs(V00.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V00.P().CloseTo(V0.P(), epsilon));
    }
}
TEST(LorentzVector, LorentzTransform2d)
{
    RANDOM RG;
    RandomUniform<> mr(0, 0.99), metrr(-5, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<2>(RG) * 1.0, metrr(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(zero());
        EXPECT_TRUE(V1 == V0);
        const auto beta = randomIsotropic<2>(RG) * mr(RG);
        const auto V2 = V0.Transform(beta).Transform(-beta);
        EXPECT_TRUE(abs(V2.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V2.P().CloseTo(V0.P(), epsilon));
        EXPECT_ANY_THROW(V0.Transform(randomIsotropic<2>(RG) * (1.0 + mr(RG))));
        const auto L0 = V0.M(), L1 = V0.Transform(beta).M(),
                   L2 = V0.Transform(beta).Transform(randomIsotropic<2>(RG) * 0.2).M();
        EXPECT_TRUE(abs(L0 - L1) < epsilon);
        EXPECT_TRUE(abs(L2 - L1) < epsilon);
        EXPECT_TRUE(abs(L0 - L2) < epsilon);
        const auto V00 = lorentz_byPM(zero(), V0.M()).Transform(-V0.Beta());
        EXPECT_TRUE(abs(V00.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V00.P().CloseTo(V0.P(), epsilon));
    }
}
TEST(LorentzVector, LorentzTransform3d)
{
    RANDOM RG;
    RandomUniform<> mr(0, 0.99), metrr(-5, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<3>(RG) * 1.0, metrr(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(Zero());
        EXPECT_TRUE(V1 == V0);
        const auto beta = randomIsotropic<3>(RG) * mr(RG);
        const auto V2 = V0.Transform(beta).Transform(-beta);
        EXPECT_TRUE(abs(V2.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V2.P().CloseTo(V0.P(), epsilon));
        EXPECT_THROW(V0.Transform(randomIsotropic<3>(RG) * (1.0 + mr(RG))), Exception<LorentzVector<>>);
        const auto L0 = V0.M(), L1 = V0.Transform(beta).M(),
                   L2 = V0.Transform(beta).Transform(randomIsotropic<3>(RG) * 0.2).M();
        EXPECT_TRUE(abs(L0 - L1) < epsilon);
        EXPECT_TRUE(abs(L2 - L1) < epsilon);
        EXPECT_TRUE(abs(L0 - L2) < epsilon);
        const auto V00 = lorentz_byPM(Zero(), V0.M()).Transform(-V0.Beta());
        EXPECT_TRUE(abs(V00.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V00.P().CloseTo(V0.P(), epsilon));
    }
}
TEST(LorentzVector, LorentzTransform4d)
{
    RANDOM RG;
    RandomUniform<> mr(0, 0.99), metrr(-5, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<4>(RG) * 1.0, metrr(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(desCartes(0., 0., 0., 0.));
        EXPECT_TRUE(V1 == V0);
        const auto beta = randomIsotropic<4>(RG) * mr(RG);
        const auto V2 = V0.Transform(beta).Transform(-beta);
        EXPECT_TRUE(abs(V2.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V2.P().CloseTo(V0.P(), epsilon));
        EXPECT_ANY_THROW(V0.Transform(randomIsotropic<4>(RG) * (1.0 + mr(RG))));
        const auto L0 = V0.M(), L1 = V0.Transform(beta).M(),
                   L2 = V0.Transform(beta).Transform(randomIsotropic<4>(RG) * 0.2).M();
        EXPECT_TRUE(abs(L0 - L1) < epsilon);
        EXPECT_TRUE(abs(L2 - L1) < epsilon);
        EXPECT_TRUE(abs(L0 - L2) < epsilon);
        const auto V00 = lorentz_byPM(desCartes(0., 0., 0., 0.), V0.M()).Transform(-V0.Beta());
        EXPECT_TRUE(abs(V00.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V00.P().CloseTo(V0.P(), epsilon));
    }
}
TEST(LorentzVector, LorentzTransform1d_more)
{
    RANDOM RG;
    RandomUniform<> M(0, 5), P(0, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<1>(RG) * P(RG), M(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(V0.Beta());
        EXPECT_TRUE(V1.P().CloseTo(desCartes(0.), epsilon));
        EXPECT_TRUE(abs(V1.M() - V0.M()) < epsilon);
    }
}
TEST(LorentzVector, LorentzTransform2d_more)
{
    RANDOM RG;
    RandomUniform<> M(0, 5), P(0, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<2>(RG) * P(RG), M(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(V0.Beta());
        EXPECT_TRUE(V1.P().CloseTo(zero(), epsilon));
        EXPECT_TRUE(abs(V1.M() - V0.M()) < epsilon);
    }
}
TEST(LorentzVector, LorentzTransform3d_more)
{
    RANDOM RG;
    RandomUniform<> M(0, 5), P(0, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<3>(RG) * P(RG), M(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(V0.Beta());
        EXPECT_TRUE(V1.P().CloseTo(Zero(), epsilon));
        EXPECT_TRUE(abs(V1.M() - V0.M()) < epsilon);
    }
}
TEST(LorentzVector, LorentzTransform4d_more)
{
    RANDOM RG;
    RandomUniform<> M(0, 5), P(0, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<4>(RG) * P(RG), M(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(V0.Beta());
        EXPECT_TRUE(V1.P().CloseTo(desCartes(0., 0., 0., 0.), epsilon));
        EXPECT_TRUE(abs(V1.M() - V0.M()) < epsilon);
    }
}

TEST(LorentzVector, decays)
{
    RANDOM RG;
    RandomUniform<> IM(2, 3), M1(0, 1), M2(0, 1), THETA(0, PI()), PHI(0, 2.0 * PI());
    for (size_t i = 0; i < 50; i++) {
        const double im = IM(RG), m1 = M1(RG), m2 = M2(RG);
        const auto C2 = binaryDecay(im, m1, m2, randomIsotropic<2>(RG));
        EXPECT_TRUE(abs((C2.first + C2.second).M() - im) < epsilon);
        EXPECT_TRUE(abs((C2.first).M() - m1) < epsilon);
        EXPECT_TRUE(abs((C2.second).M() - m2) < epsilon);
        EXPECT_TRUE((C2.first.P() + C2.second.P()).M() < epsilon);
        const auto C3 = binaryDecay(im, m1, m2, randomIsotropic<3>(RG));
        EXPECT_TRUE(abs((C3.first + C3.second).M() - im) < epsilon);
        EXPECT_TRUE(abs(C3.first.M() - m1) < epsilon);
        EXPECT_TRUE(abs(C3.second.M() - m2) < epsilon);
        EXPECT_TRUE((C3.first.P() + C3.second.P()).M() < epsilon);
    }
}
TEST(LorentzVector, decays2)
{
    RANDOM RG;
    RandomUniform<> IM(2, 3), M1(0, 1), M2(0, 1), THETA(0, PI()), PHI(0, 2.0 * PI());
    for (size_t i = 0; i < 50; i++) {
        const double im = IM(RG), m1 = M1(RG), m2 = M2(RG);
        const auto C2 = binaryDecay(im, m1, m2, randomIsotropic<2>(RG));
        EXPECT_TRUE(abs((C2.first + C2.second).M() - im) < epsilon);
        EXPECT_TRUE(abs((C2.first).M() - m1) < epsilon);
        EXPECT_TRUE(abs((C2.second).M() - m2) < epsilon);
        EXPECT_TRUE((C2.first.P() + C2.second.P()).M() < epsilon);
        const auto C3 = binaryDecay(im, m1, m2, randomIsotropic<3>(RG));
        EXPECT_TRUE(abs((C3.first + C3.second).M() - im) < epsilon);
        EXPECT_TRUE(abs(C3.first.M() - m1) < epsilon);
        EXPECT_TRUE(abs(C3.second.M() - m2) < epsilon);
        EXPECT_TRUE((C3.first.P() + C3.second.P()).M() < epsilon);
    }
}

TEST(Vector, direction1d)
{
    RANDOM RG;
    RandomUniform<> V(-10.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto d = direction();
        const auto v = desCartes(V(RG));
        EXPECT_TRUE((v - d * v.x()).M() < epsilon);
    }
}
TEST(Vector, direction2d)
{
    RANDOM RG;
    RandomUniform<> PHI(-PI(), PI()), M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto phi = PHI(RG);
        const auto phi2 = direction(direction(phi) * M(RG)).phi();
        EXPECT_TRUE(abs(phi - phi2) < epsilon);
    }
}
TEST(Vector, direction3d)
{
    RANDOM RG;
    RandomUniform<> THETA(0.0, PI()), PHI(-PI(), PI()), M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto angles = direction(PHI(RG), THETA(RG));
        const auto angles2 = direction(angles * M(RG));
        EXPECT_TRUE(abs(angles.phi() - angles2.phi()) < epsilon);
        EXPECT_TRUE(abs(angles.th() - angles2.th()) < epsilon);
    }
}
TEST(Vector, direction4d)
{
    RANDOM RG;
    RandomUniform<> THETA(0.0, PI()), PHI(-PI(), PI()), M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto angles = direction(PHI(RG), THETA(RG), THETA(RG));
        const auto angles2 = direction(angles * M(RG));
        EXPECT_TRUE(abs(angles.phi() - angles2.phi()) < epsilon);
        EXPECT_TRUE(abs(angles.th<1>() - angles2.th<1>()) < epsilon);
        EXPECT_TRUE(abs(angles.th<2>() - angles2.th<2>()) < epsilon);
    }
}
TEST(Vector, direction5d)
{
    RANDOM RG;
    RandomUniform<> THETA(0.0, PI()), PHI(-PI(), PI()), M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto angles = direction(PHI(RG), THETA(RG), THETA(RG), THETA(RG));
        const auto angles2 = direction(angles * M(RG));
        EXPECT_TRUE(abs(angles.phi() - angles2.phi()) < epsilon);
        EXPECT_TRUE(abs(angles.th<1>() - angles2.th<1>()) < epsilon);
        EXPECT_TRUE(abs(angles.th<2>() - angles2.th<2>()) < epsilon);
        EXPECT_TRUE(abs(angles.th<3>() - angles2.th<3>()) < epsilon);
    }
}

TEST(Vector, rotations1d)
{
    RANDOM RG;
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto dir = randomIsotropic<1>(RG);
        const auto m = M(RG);
        const auto V1 = dir * m;
        const auto V2 = dir.Rotations() * (axis<1, 1>() * m);
        EXPECT_TRUE(V1.CloseTo(V2, epsilon));
        const auto V3 = dir.AntiRotations() * V2;
        EXPECT_TRUE(V3.CloseTo(axis<1, 1>()*m, epsilon));
    }
}
TEST(Vector, rotations2d)
{
    RANDOM RG;
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto dir = randomIsotropic<2>(RG);
        const auto m = M(RG);
        const auto V1 = dir * m;
        const auto V2 = dir.Rotations() * (x() * m);
        EXPECT_TRUE(V1.CloseTo(V2, epsilon));
        const auto V3 = dir.AntiRotations() * V2;
        EXPECT_TRUE(V3.CloseTo(x()*m, epsilon));
    }
}
TEST(Vector, rotations3d)
{
    RANDOM RG;
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto dir = randomIsotropic<3>(RG);
        const auto m = M(RG);
        const auto V1 = dir * m;
        const auto V2 = dir.Rotations() * (Z() * m);
        EXPECT_TRUE(V1.CloseTo(V2, epsilon));
        const auto V3 = dir.AntiRotations() * V2;
        EXPECT_TRUE(V3.CloseTo(Z()*m, epsilon));
    }
}
TEST(Vector, rotations4d)
{
    RANDOM RG;
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto dir = randomIsotropic<4>(RG);
        const auto m = M(RG);
        const auto V1 = dir * m;
        const auto V2 = dir.Rotations() * (axis<4, 4>() * m);
        EXPECT_TRUE(V1.CloseTo(V2, epsilon));
        const auto V3 = dir.AntiRotations() * V2;
        EXPECT_TRUE(V3.CloseTo(axis<4, 4>()*m, epsilon));
    }
}
TEST(Vector, rotations5d)
{
    RANDOM RG;
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto dir = randomIsotropic<5>(RG);
        const auto m = M(RG);
        const auto V1 = dir * m;
        const auto V2 = dir.Rotations() * (axis<5, 5>() * m);
        EXPECT_TRUE(V1.CloseTo(V2, epsilon));
        const auto V3 = dir.AntiRotations() * V2;
        EXPECT_TRUE(V3.CloseTo(axis<5, 5>()*m, epsilon));
    }
}

TEST(Vector, rotations3d_plane)
{
    RANDOM RG;
    RandomUniform<> M(-10.0, 10.0), TH(-PI<>(), PI<>());
    for (size_t i = 0; i < 50; i++) {
        const auto R = randomIsotropic<3>(RG).Rotations();
        const auto v1 = R * desCartes(M(RG), M(RG), 0.);
        const auto v2 = R * desCartes(M(RG), M(RG), 0.);
        const auto v3 = R * desCartes(M(RG), M(RG), 0.);
        EXPECT_TRUE(abs((v1 ^ v2)*v3) < epsilon);
    }
    for (size_t i = 0; i < 50; i++) {
        const auto R = randomIsotropic<3>(RG).Rotations();
        const auto v1 = R * desCartes(0., M(RG), M(RG));
        const auto v2 = R * desCartes(0., M(RG), M(RG));
        const auto v3 = R * desCartes(0., M(RG), M(RG));
        EXPECT_TRUE(abs((v1 ^ v2)*v3) < epsilon);
    }
    for (size_t i = 0; i < 10; i++) {
        const auto R = randomIsotropic<3>(RG).Rotations();
        const auto v1 = R * desCartes(M(RG), 0., M(RG));
        const auto v2 = R * desCartes(M(RG), 0., M(RG));
        const auto v3 = R * desCartes(M(RG), 0., M(RG));
        EXPECT_TRUE(abs((v1 ^ v2)*v3) < epsilon);
    }
    for (size_t i = 0; i < 50; i++) {
        const auto R1 = randomIsotropic<3>(RG).Rotations();
        const auto v1 = desCartes(M(RG), 0., M(RG));
        const auto R2 = Rotation(direction(v1), TH(RG));
        const auto R = R1 * R2;
        const auto v2 = desCartes(M(RG), 0., M(RG));
        const auto v3 = desCartes(M(RG), 0., M(RG));
        EXPECT_TRUE(abs(((R * v1) ^ (R * v2)) * (R * v3)) < epsilon);
    }
}
TEST(VectorTransformation, zero1)
{
    EXPECT_EQ(
        ZERO<1>(), line(0.)
    );
    EXPECT_EQ(
        ZERO<2>(), lines(
            desCartes(0., 0.),
            desCartes(0., 0.)
        )
    );
    EXPECT_EQ(
        ZERO<3>(), lines(
            desCartes(0., 0., 0.),
            desCartes(0., 0., 0.),
            desCartes(0., 0., 0.)
        )
    );
}
TEST(VectorTransformation, zero2)
{
    EXPECT_EQ(
        ZERO<1>(), column(0.)
    );
    EXPECT_EQ(
        ZERO<2>(), columns(
            desCartes(0., 0.),
            desCartes(0., 0.)
        )
    );
    EXPECT_EQ(
        ZERO<3>(), columns(
            desCartes(0., 0., 0.),
            desCartes(0., 0., 0.),
            desCartes(0., 0., 0.)
        )
    );
}
TEST(VectorTransformation, one1)
{
    EXPECT_EQ(
        ONE<1>(), line(1.)
    );
    EXPECT_EQ(
        ONE<2>(), lines(
            desCartes(1., 0.),
            desCartes(0., 1.)
        )
    );
    EXPECT_EQ(
        ONE<3>(), lines(
            desCartes(1., 0., 0.),
            desCartes(0., 1., 0.),
            desCartes(0., 0., 1.)
        )
    );
}
TEST(VectorTransformation, one2)
{
    EXPECT_EQ(
        ONE<1>(), column(1.)
    );
    EXPECT_EQ(
        ONE<2>(), columns(
            desCartes(1., 0.),
            desCartes(0., 1.)
        )
    );
    EXPECT_EQ(
        ONE<3>(), columns(
            desCartes(1., 0., 0.),
            desCartes(0., 1., 0.),
            desCartes(0., 0., 1.)
        )
    );
}
TEST(VectorTransformation, lines_columns)
{
    EXPECT_EQ(
        lines(
            desCartes(1, 2, 3),
            desCartes(4, 5, 6),
            desCartes(7, 8, 9)
        ),
        columns(
            desCartes(1, 4, 7),
            desCartes(2, 5, 8),
            desCartes(3, 6, 9)
        )
    );
    EXPECT_EQ(
        lines(
            desCartes(1, 2, 3),
            desCartes(4, 5, 6)
        ),
        columns(
            desCartes(1, 4),
            desCartes(2, 5),
            desCartes(3, 6)
        )
    );
}
TEST(VectorTransformation, elements)
{
    const auto M=lines(
	desCartes(1,2,3,4),
	desCartes(5,6,7,8),
	desCartes(9,0,1,2)
    );
    const auto&e11=M.element<1,1>();
    const auto&e12=M.element<1,2>();
    const auto&e13=M.element<1,3>();
    const auto&e14=M.element<1,4>();
    const auto&e21=M.element<2,1>();
    const auto&e22=M.element<2,2>();
    const auto&e23=M.element<2,3>();
    const auto&e24=M.element<2,4>();
    const auto&e31=M.element<3,1>();
    const auto&e32=M.element<3,2>();
    const auto&e33=M.element<3,3>();
    const auto&e34=M.element<3,4>();
    EXPECT_EQ(e11,1);
    EXPECT_EQ(e12,2);
    EXPECT_EQ(e13,3);
    EXPECT_EQ(e14,4);
    EXPECT_EQ(e21,5);
    EXPECT_EQ(e22,6);
    EXPECT_EQ(e23,7);
    EXPECT_EQ(e24,8);
    EXPECT_EQ(e31,9);
    EXPECT_EQ(e32,0);
    EXPECT_EQ(e33,1);
    EXPECT_EQ(e34,2);
}
TEST(VectorTransformation, mul1)
{
    RANDOM RG;
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto v = randomIsotropic<1>(RG) * M(RG);
        EXPECT_TRUE(desCartes(0.).CloseTo(ZERO<1>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<1>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<1>() * (ONE<1>()*v), epsilon));
        const auto R1 = randomIsotropic<1>(RG).Rotations();
        const auto R2 = randomIsotropic<1>(RG).Rotations();
        EXPECT_TRUE((R1 * (R2 * v)).CloseTo((R1 * R2)*v, epsilon));
    }
}
TEST(VectorTransformation, mul2)
{
    RANDOM RG;
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto v = randomIsotropic<2>(RG) * M(RG);
        EXPECT_TRUE(zero().CloseTo(ZERO<2>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<2>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<2>() * (ONE<2>()*v), epsilon));
        const auto R1 = randomIsotropic<2>(RG).Rotations();
        const auto R2 = randomIsotropic<2>(RG).Rotations();
        EXPECT_TRUE((R1 * (R2 * v)).CloseTo((R1 * R2)*v, epsilon));
    }
}
TEST(VectorTransformation, mul3)
{
    RANDOM RG;
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto v = randomIsotropic<3>(RG) * M(RG);
        EXPECT_TRUE(Zero().CloseTo(ZERO<3>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<3>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<3>() * (ONE<3>()*v), epsilon));
        const auto R1 = randomIsotropic<3>(RG).Rotations();
        const auto R2 = randomIsotropic<3>(RG).Rotations();
        EXPECT_TRUE((R1 * (R2 * v)).CloseTo((R1 * R2)*v, epsilon));
    }
}
