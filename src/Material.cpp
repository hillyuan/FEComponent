#include "material/MaterialBase.h"

/*Constitutive_type IsotropicElastic::ConstitutiveMatrix(std::vector<double>& depends)
{
	Constitutive_type matrix;
	double e, nu, ee;

    e  = this->E;
    nu = this->nu;

    ee = e / ( ( 1. + nu ) * ( 1. - 2. * nu ) );

    matrix(1, 1) =  1. - nu;
    matrix(1, 2) =  nu;
    matrix(1, 3) =  nu;
    matrix(2, 1) =  nu;
    matrix(2, 2) =  1. - nu;
    matrix(2, 3) =  nu;
    matrix(3, 1) =  nu;
    matrix(3, 2) =  nu;
    matrix(3, 3) =  1. - nu;

    matrix(4, 4) =  ( 1. - 2. * nu ) * 0.5;
    matrix(5, 5) =  ( 1. - 2. * nu ) * 0.5;
    matrix(6, 6) =  ( 1. - 2. * nu ) * 0.5;
}	*/

/*
void
HyperElasticMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode, GaussPoint *gp, TimeStep *tStep)
{
    double J2, c11, c22, c33, c12, c13, c23, A, B;
    FloatMatrix C(3, 3);
    FloatMatrix invC;

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    C.at(1, 1) = 1. + 2. * status->giveTempStrainVector().at(1);
    C.at(2, 2) = 1. + 2. * status->giveTempStrainVector().at(2);
    C.at(3, 3) = 1. + 2. * status->giveTempStrainVector().at(3);
    C.at(2, 3) = C.at(3, 2) = status->giveTempStrainVector().at(4);
    C.at(1, 3) = C.at(3, 1) = status->giveTempStrainVector().at(5);
    C.at(1, 2) = C.at(2, 1) = status->giveTempStrainVector().at(6);

    invC.beInverseOf(C);
    J2 = C.giveDeterminant();

    c11 = invC.at(1, 1);
    c22 = invC.at(2, 2);
    c33 = invC.at(3, 3);
    c12 = invC.at(1, 2);
    c13 = invC.at(1, 3);
    c23 = invC.at(2, 3);

    A = ( K - 2. / 3. * G ) * J2;
    B = -( K - 2. / 3. * G ) * ( J2 - 1. ) + 2. * G;

    answer.resize(6, 6);

    answer.at(1, 1) = ( A + B ) * c11 * c11;
    answer.at(2, 2) = ( A + B ) * c22 * c22;
    answer.at(3, 3) = ( A + B ) * c33 * c33;
    answer.at(4, 4) = A * c23 * c23 + B / 2. * ( c22 * c33 + c23 * c23 );
    answer.at(5, 5) = A * c13 * c13 + B / 2. * ( c11 * c33 + c13 * c13 );
    answer.at(6, 6) = A * c12 * c12 + B / 2. * ( c11 * c22 + c12 * c12 );
    answer.at(1, 2) = answer.at(2, 1) = A * c11 * c22 + B * c12 * c12;
    answer.at(1, 3) = answer.at(3, 1) = A * c11 * c33 + B * c13 * c13;
    answer.at(1, 4) = answer.at(4, 1) = A * c11 * c23 + B * c12 * c13;
    answer.at(1, 5) = answer.at(5, 1) = A * c11 * c13 + B * c11 * c13;
    answer.at(1, 6) = answer.at(6, 1) = A * c11 * c12 + B * c11 * c12;
    answer.at(2, 3) = answer.at(3, 2) = A * c22 * c33 + B * c23 * c23;
    answer.at(2, 4) = answer.at(4, 2) = A * c22 * c23 + B * c22 * c23;
    answer.at(2, 5) = answer.at(5, 2) = A * c22 * c13 + B * c12 * c23;
    answer.at(2, 6) = answer.at(6, 2) = A * c22 * c12 + B * c22 * c12;
    answer.at(3, 4) = answer.at(4, 3) = A * c33 * c23 + B * c33 * c23;
    answer.at(3, 5) = answer.at(5, 3) = A * c33 * c13 + B * c33 * c13;
    answer.at(3, 6) = answer.at(6, 3) = A * c33 * c12 + B * c13 * c23;
    answer.at(4, 5) = answer.at(5, 4) = A * c23 * c13 + B / 2. * ( c12 * c33 + c13 * c23 );
    answer.at(4, 6) = answer.at(6, 4) = A * c23 * c12 + B / 2. * ( c12 * c23 + c22 * c13 );
    answer.at(5, 6) = answer.at(6, 5) = A * c13 * c12 + B / 2. * ( c11 * c23 + c12 * c13 );
}*/