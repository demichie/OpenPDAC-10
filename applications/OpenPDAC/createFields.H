#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"

Info<< "Creating phaseSystem\n" << endl;

autoPtr<phaseSystem> fluidPtr
(
    phaseSystem::New(mesh)
);
phaseSystem& fluid = fluidPtr();
phaseSystem::phaseModelList& phases = fluid.phases();

volScalarField& p = phases[0].thermoRef().p();
volScalarField& T1 = phases[0].thermoRef().T();
volScalarField& T2 = phases[1].thermoRef().T();


volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::fvConstraints& fvConstraints = fluid.fvConstraints(mesh);


volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.U()
);

U.write();

if (fluid.found("pMin"))
{
    IOWarningInFunction(fluid)
        << "Pressure limits, pMin and pMax, are now read from "
        << pimple.dict().name() << endl;
}


volScalarField ph_rgh
    (
        IOobject
        (
            "ph_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
#include "phrghEqn.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


fluid.correctThermo();
rho = fluid.rho();

hRef.value() = 0.0;

Info << "hRef " << hRef.value() << endl;

#include "gh.H"

// this field has been written by phrghEqn.H, 
// if hydrostaticInitialisation is true
/*volScalarField ph_rgh
(
    IOobject
    (
        "ph_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);*/

p_rgh = ph_rgh;

Info<< "min p " << min(p).value() <<
  	               " max p " << max(p).value() << endl;
Info<< "min p_rgh " << min(p_rgh).value() <<
   	               " max p_rgh " << max(p_rgh).value() << endl;
Info<< "min rho " << min(rho).value() <<
   	               " max rho " << max(rho).value() << endl;

pressureReference pressureReference
(
    p,
    p_rgh,
    pimple.dict(),
    fluid.incompressible()
);

if (fluid.incompressible())
{
    p = p_rgh + fluid.rho()*gh;
}

if (p_rgh.needReference() && fluid.incompressible())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
}

p_rgh = p - fluid.rho()*gh;
mesh.schemes().setFluxRequired(p_rgh.name());

label carrierIdx;

carrierIdx = 0;

forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];
    if (!phase.incompressible())
    {
    	Info << phasei << " compressible" << endl;
    	carrierIdx = phasei;
    }
	else
    {
    	Info << phasei << " incompressible" << endl;
    }

    // Info << word(phase.thermo().subDict("thermoType").lookup("equationOfState")) << endl;
}
    
volScalarField muC
(
    IOobject
    (
        "muC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    phases[carrierIdx].thermo().mu() 
);    


Info<< "min muC " << min(muC).value() << " max muC " << max(muC).value() << endl;

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    muC * pow( 1.0 - ( 1.0 - phases[carrierIdx] ) / 0.62 , -1.55)
);

Info<< "min mu " << min(mu).value() << " max mu " << max(mu).value() << endl;

Info<< "Constructing clouds" << endl;
parcelCloudList clouds(rho, U, mu, g);

PtrList<volScalarField> rAUs;
PtrList<surfaceScalarField> rAUfs;
