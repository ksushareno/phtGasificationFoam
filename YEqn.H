/**/
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi)")
    )
);

{
	volScalarField Yt(0.0*Y[0]);
	forAll(Y, i)
	{
        if(i != inertIndex)
        {
            volScalarField& Yi = Y[i];
            
            fvScalarMatrix YiEqn
            (
                fvm::ddt(a_g, rho,Yi)
            //    + fvm::div(phi,Yi)
                + mvConvection->fvmDiv(phi, Yi)
                ==
                ro[i-1] * Yc[i-1] * a_cf0 * W[i-1]
                + fvOptions(rho, Yi)
            );
            
            YiEqn.relax();
            fvOptions.constrain(YiEqn);
            YiEqn.solve(mesh.solver("Yi"));
            fvOptions.correct(Yi);
            
            Yi.max(0.0);
            Yt += Yi;
        }
	}
        
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
