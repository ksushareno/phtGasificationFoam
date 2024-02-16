{
	volScalarField& Tc = thermoS.T();

    fvScalarMatrix TcEqn
    (
     //   solidCapacity*fvm::ddt(Tc)
	fvm::ddt(solidCapacity, Tc)
	    + fvm::Sp(c_cf*Mdot, Tc)
        - fvm::laplacian(thermo.kappa() * (1 - a_g), Tc)
        + fvm::Sp(HTcoef * (1 - a_g) / a_cf0, Tc)
        - HTcoef * (1 - a_g) / a_cf0 * thermo.T()
        ==
		- Qdot
        + fvOptions(solidCapacity, Tc)
    );

    TcEqn.relax();
    fvOptions.constrain(TcEqn);
    TcEqn.solve();
    fvOptions.correct(Tc);
}
