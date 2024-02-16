Info<< "Reading Fuel Properties\n" << endl;
IOdictionary fuelProperties
(
    IOobject
    (
        "fuelProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//Activation Energy

IOdictionary thermoDict
(
    IOobject
    (
        "thermophysicalProperties.solid",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

wordList speciesSolid(thermoDict.get<wordList>("species"));
PtrList<dimensionedScalar> E(speciesSolid.size());
PtrList<dimensionedScalar> k2(speciesSolid.size());
PtrList<dimensionedScalar> Q(speciesSolid.size());
PtrList<dimensionedScalar> ro(speciesSolid.size());
PtrList<dimensionedScalar> cp(speciesSolid.size());
PtrList<volScalarField> Yc(speciesSolid.size());
PtrList<volScalarField> W(speciesSolid.size());
PtrList<volScalarField> eta(speciesSolid.size());
forAll(speciesSolid, i){
	E.set
    (
        i,
        new dimensionedScalar
        (
            dimensionSet(1, 2, -2, 0, -1, 0, 0),
			thermoDict.subDict(speciesSolid[i]).subDict("reaction").get<scalar>("E")
        )
	);
	
	Q.set
    (
        i,
        new dimensionedScalar
        (
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
			thermoDict.subDict(speciesSolid[i]).subDict("reaction").get<scalar>("Q")
        )
	);
	
	k2.set
    (
        i,
        new dimensionedScalar
        (
            dimensionSet(0, 0, -1, 0, 0, 0, 0),
			thermoDict.subDict(speciesSolid[i]).subDict("reaction").get<scalar>("k")
        )
	);
	
	ro.set
    (
        i,
        new dimensionedScalar
        (
            dimensionSet(1, -3, 0, 0, 0, 0, 0),
			thermoDict.subDict(speciesSolid[i]).subDict("equationOfState").get<scalar>("rho")
        )
	);
	
	cp.set
    (
        i,
        new dimensionedScalar
        (
            dimensionSet(0, 2, -2, -1, 0, 0, 0),
			thermoDict.subDict(speciesSolid[i]).subDict("thermodynamics").get<scalar>("Cp")
        )
	);
	
	Yc.set
	(
        i,
		new volScalarField
		(
		    IOobject
			(
		        IOobject::groupName(speciesSolid[i], solidPhaseName),
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
		    ),
			mesh
		)
	);
	
	W.set
	(
        i,
		new volScalarField
		(
		    IOobject
			(
			    IOobject::groupName("W",IOobject::groupName(speciesSolid[i], solidPhaseName)),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar(dimless/dimTime, Zero)
		)
	);
	
	eta.set
	(
        i,
		new volScalarField
		(
		    IOobject
			(
			    IOobject::groupName("eta",IOobject::groupName(speciesSolid[i], solidPhaseName)),
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
			),
			mesh
		)
	);
				 
}

volScalarField Eta
(
    IOobject
    (
        "eta.solid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

//Overall heat-release intensity
volScalarField Qdot
(
    IOobject
    (
        "Qdot.solid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);
//Right-hand side in continuity eq

volScalarField Mdot
(
    IOobject
    (
        "Mdot.solid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
);

volScalarField rho_cf0("rhoEf0.solid",thermoS.rho() * a_cf0);

volScalarField rho_cf
(
    IOobject
    (
        "rhoEf.solid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
//    mesh,
    rho_cf0
);

dimensionedScalar T_norm
(
    "T_norm",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    273
);



