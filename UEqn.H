// Solve the Momentum equation

//MRF.correctBoundaryVelocity(U);
//volVectorField& resistance = KD * mu + rho * F * mag(U);
/**/

//Foam::pos0(Tc-Tm)
/*KdField = KD / 
(
(1+10*Eta)*Foam::neg(Eta-0.8)
+
(9+5*(KD/KdEnd)*(Eta-0.8))*Foam::pos0(Eta-0.8)
);*/
//KdField = KD;
KdField = 150 * (1 - a_g) * (1 - a_g) * KD / (a_g * a_g * a_g); //для формулы Эргуна
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(a_g, rho, U) 
    + fvm::div(phi, U)
    + a_g * rho * chi * (1 - a_g) * (
    fvc::ddt(U)
    + (U & fvc::grad(U)) 
//    fvm::ddt(U)
//	+ fvm::div(fvc::flux(U), U) - fvm::Sp(fvc::div(fvc::flux(U)), U)
    )
    + fvm::Sp(a_g*a_g*KdField * thermo.mu() + rho * F * mag(U), U)
//  + MRF.DDt(rho, U)
//  - fvm::laplacian(mu, U)
 ==
    fvOptions(a_g, rho, U)
);

//Формулировка относительно скорости фильтрации
/*tmp<fvVectorMatrix> tUEqn
(
    scalar(0)*fvm::ddt(a_g, rho, U) 
    + scalar(0)*fvm::div(phi, U)
    + 0*a_g * rho * chi * (1 - a_g) * (
    fvm::ddt(U)
    + fvm::div(fvc::flux(U), U) 
    - fvm::Sp(fvc::div(fvc::flux(U)), U)
    )
    + fvm::Sp(KD * thermo.mu() + rho * F * mag(U), U)
//  + MRF.DDt(rho, U)
//  - fvm::laplacian(mu, U)
 ==
    fvOptions(a_g, rho, U)
);*/
//tUEqn+=scalar(0);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -a_g *fvc::grad(p));
//    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
//    K = 0.5*magSqr(U);
}
