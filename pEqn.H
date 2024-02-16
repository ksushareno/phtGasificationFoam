if (!pimple.SIMPLErho())
{
    rho = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU("rAU", 1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(a_g*a_g*rho*rAU));
//surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
//volVectorField HbyA(rAU*UEqn.H());

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(a_g*rho)*fvc::flux(HbyA)
//	fvc::interpolate(rho)*fvc::flux(HbyA)
//  + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf))
//    + rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf)
);

fvc::makeRelative(phiHbyA, rho, U);
//MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
//constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);
constrainPressure(p, rho, a_g * U, phiHbyA, rhorAUf);

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(a_g, rho) + psi*correction(fvm::ddt(p))
	//    a_g*fvc::ddt(rho) + rho*fvc::ddt(a_g) + a_g * psi*correction(fvm::ddt(p))
      + fvc::div(phiHbyA) 
     ==
	    Mdot
      + fvOptions(psi, p, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(/*fvc::interpolate(a_g)*/rhorAUf, p));

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
		//	U = HbyA + fvc::reconstruct(pEqn.flux()/fvc::interpolate(a_g*rho));
		//	U.correctBoundaryConditions();
        //    fvOptions.correct(U);
		//volScalarField resP("resP", pEqn.A() * p - pEqn.H());
		//if (runTime.writeTime()) {resP.write();}
        }
    }

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - a_g*rAU*fvc::grad(p);
//U = fvc::reconstruct(phiHbyA/fvc::interpolate(a_g*rho)) - a_g*rAU*fvc::grad(p);
//U = fvc::reconstruct(phi/fvc::interpolate(a_g*rho));
U.correctBoundaryConditions();
fvOptions.correct(U);
//volVectorField U1("U1", fvc::reconstruct(phi/fvc::interpolate(a_g*rho)));
//if (runTime.writeTime()) {U1.write(); HbyA.write(); rAU.write(); phiHbyA.write();}
//K = 0.5*magSqr(U);

volVectorField resU("resU", UEqn.A() * U - UEqn.H() + a_g * fvc::grad(p));
if (runTime.writeTime()) {resU.write();}

if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();
}

thermo.correctRho(psi*p - psip0, rhoMin, rhoMax) ;

//Info<< "\nSolving density for continuity error checking" << endl;
//#include "rhoEqn.H"
//#include "compressibleContinuityErrs.H"

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
