#include "createRDeltaT.H"

word fluidPhaseName = "fluid";
word solidPhaseName = "solid";

Info<< "Reading fluid thermophysical properties\n" << endl;

autoPtr<psiReactionThermo> pThermo
(
    psiReactionThermo::New(mesh,fluidPhaseName)
);
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.get<word>("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(Y, i)
{
    fields.add(Y[i]);
}

volScalarField rho
(
    IOobject
    (
        "rho.fluid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading solid thermophysical properties\n" << endl;

autoPtr<solidThermo> sThermo
(
    solidThermo::New(mesh,solidPhaseName)
);
solidThermo& thermoS = sThermo();
thermoS.validate(args.executable(), "h", "e");

#include "readPorosityProperties.H"
#include "readFuelProperties.H"



volScalarField c_cf = thermoS.Cp();

volScalarField solidCapacity
(
    IOobject
    (
        "solidCapacity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
//    mesh,
    rho_cf0 * c_cf
);



dimensionedScalar rho_ci
(
    "rho_ci",
    dimensionSet(1, -3, 0, 0, 0, 0, 0),
	1150
    //rho_cf0 * (1.0 - (a_g0 + a_cf0))
);

rho_ci = rho_ci * (1.0 - (a_g0 + a_cf0));

dimensionedScalar c_ci
(
    "c_ci",
    dimensionSet(0, 2, -2, -1, 0, 0, 0),
    1000
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);




#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());


const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

#include "createFvOptions.H"
