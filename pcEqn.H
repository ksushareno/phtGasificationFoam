if (!pimple.SIMPLErho())
{
    rho = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU(1.0/UEqn.A());
volScalarField rAtU(1.0/(1.0/rAU - UEqn.H1()));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::interpolate(a_g*rho)*fvc::flux(HbyA)
//		fvc::interpolate(rho)*fvc::flux(HbyA)
      /*+ MRF.zeroFilter
        (
            fvc::interpolate(rho*rAU)*fvc::ddtCorr(rho, U, phi, rhoUf)
        )*/
        + fvc::interpolate(rho*rAU)*fvc::ddtCorr(rho, U, phi, rhoUf)
    )
);

fvc::makeRelative(phiHbyA, rho, U);
//MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

volScalarField rhorAtU("rhorAtU", a_g*rho*rAtU);

// Update the pressure BCs to ensure flux consistency
//constrainPressure(p, rho, U, phiHbyA, rhorAtU, MRF);
constrainPressure(p, rho, U, phiHbyA, rhorAtU);

    phiHbyA += fvc::interpolate(a_g*rho*(rAtU - rAU)*a_g)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU)*a_g*fvc::grad(p);

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(a_g, rho) + psi*correction(fvm::ddt(p))
      + fvc::div(phiHbyA)
     ==
        scalar(0)*rho_cf0 * W_rate
      + fvOptions(psi, p, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAtU, p));

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - a_g*rAtU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
//K = 0.5*magSqr(U);

if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();
}

thermo.correctRho(psi*p - psip0, rhoMin, rhoMax) ;
rho = thermo.rho();

// Correct rhoUf if the mesh is moving
fvc::correctRhoUf(rhoUf, rho, U, phi);
/*
if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);

    if (mesh.moving())
    {
        dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
    }
}*/
