/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    etaEqn

Description
    

\*---------------------------------------------------------------------------*/

{
	Eta = scalar(0);
	rho_cf = rho_cf0;
	Qdot = scalar(0)*Qdot;
	Mdot = scalar(0)*Mdot;
    forAll(eta, i)
	{
        volScalarField& etai = eta[i];
            
        fvScalarMatrix etaiEqn
        (
            fvm::ddt(etai)
            + fvm::Sp(k2[i] * exp(-E[i]/Runiv/thermoS.T()), etai)
            - k2[i] * exp(-E[i]/Runiv/thermoS.T())
        );
        etaiEqn.solve();
		Eta += etai*Yc[i];
		W[i] = (1 - etai) * k2[i] * exp(-E[i]/Runiv/thermoS.T());
		rho_cf -= etai * ro[i]*Yc[i]*a_cf0;
		Qdot += Q[i]*ro[i]*a_cf0*Yc[i]*W[i];
		Mdot += a_cf0*Yc[i]*ro[i]*W[i];
		
	}
    a_g = a_g0 + a_cf0* Eta;
    solidCapacity = rho_cf * c_cf  + rho_ci * c_ci;
	
	/**/
}

// ************************************************************************* //
