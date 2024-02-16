Info<< "Reading Porosity Properties\n" << endl;
IOdictionary porosityProperties
(
    IOobject
    (
        "porosityProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Porosity

dimensionedScalar a_g0
(
    "a0",
    dimless,
    porosityProperties
);

dimensionedScalar a_cf0
(
    "a_cf0",
//    dimless,
    (1 - a_g0) * 0.999999
);

dimensionedScalar a_cpEnd
(
    "a_cpEnd",
    dimless,
    0
);

// Porous resistance properties

dimensionedScalar chi
(
    "chi",
    dimless,
    porosityProperties
);

dimensionedScalar KD
(
    "Kd0",
    dimensionSet(0, -2, 0, 0, 0, 0, 0),
    porosityProperties
);

dimensionedScalar KdEnd
(
    "KdEnd",
    dimensionSet(0, -2, 0, 0, 0, 0, 0),
    0.1
);

dimensionedScalar F
(
    "F",
    dimensionSet(0, -1, 0, 0, 0, 0, 0),
    porosityProperties
);


dimensionedScalar Runiv
(
    "Runiv",
    dimensionSet(1, 2, -2, -1, -1, 0, 0),
    8.31446
);

// Heat transfer coefficient

dimensionedScalar HTcoef
(
    "Kh",
    dimensionSet(1, -1, -3, -1, 0, 0, 0),
    porosityProperties
);


volScalarField a_g
(
    IOobject
    (
        "a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    a_g0
);


volScalarField KdField
(
    IOobject
    (
        "rPerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    KD
);
