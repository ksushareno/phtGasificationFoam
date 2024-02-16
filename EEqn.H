{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(a_g, rho, he) 
        + fvm::div(phi, he)
        - fvm::laplacian(a_g * thermo.alpha(), he)
	    -fvm::Sp(Mdot, he)
     ==
	   HTcoef * (1 - a_g) / a_cf0 * (thermoS.T() - thermo.T())
       + fvOptions(a_g, rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
}
