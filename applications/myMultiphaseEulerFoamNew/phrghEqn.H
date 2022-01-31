Info<< "Hydrostatic Initialisation" << endl;
dimensionedScalar Uh_ref("Uh_ref",dimVelocity,scalar(1.0));
dimensionedScalar pRef("pRef",dimensionSet(1, -1, -2, 0, 0),scalar(101325.0));
dimensionedScalar hRef("hRef",dimLength,scalar(0.0));

if (pimple.dict().lookupOrDefault<bool>("hydrostaticInitialisation", false))
{
    const fvMesh& mesh = p.mesh();
		
    volVectorField centres = mesh.C();
    volScalarField x = centres.component(0);
    volScalarField y = centres.component(1);
    volScalarField z = centres.component(2);        

    dimensionedScalar xMin = min(mesh.Cf().component(0));
    dimensionedScalar xMax = max(mesh.Cf().component(0));
    dimensionedScalar yMin = min(mesh.Cf().component(1));
    dimensionedScalar yMax = max(mesh.Cf().component(1));
    dimensionedScalar zMin = min(mesh.Cf().component(2));
    dimensionedScalar zMax = max(mesh.Cf().component(2));
        
    if (!mesh.time().restart())
    {

        volScalarField& ph = regIOobject::store
       	(
            new volScalarField
            (
                IOobject
                (
                    "ph",
                    "0",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
               	),
               	mesh
            )
       	);

	label nFixValPatches = 0;
	label bdryPref = 0;

       	forAll(ph.boundaryField(), bdryID)     
   	{    		
	    Info << "magSfBdry" << mesh.Sf().boundaryField()[bdryID][0] << endl;
   	    //if ( ( mag(mesh.Sf().boundaryField()[bdryID][0]^g).value() == 0.0 ) 
   	    if ( ( mag(mesh.Sf().boundaryField()[bdryID][0]^g).value() <= 
   	           1.0e-8*mag(mesh.Sf().boundaryField()[bdryID][0])*mag(g).value() ) 
   				&& ( ph.boundaryField()[bdryID].type() == "fixedValue" ) )
   	    {
   	        nFixValPatches +=1;
   	        const fvPatchScalarField& ph_p = ph.boundaryField()[bdryID];
    			    
		if ( min(ph_p) == max(ph_p) )
		{
		    pRef.value() = min(ph_p);
		    bdryPref = bdryID;
		    Info << "pRef " << pRef.value() << endl;
		}					

   	    }
	}

	Info << "nFixValPatches " << nFixValPatches << endl;

	if ( nFixValPatches != 1 )
	{
	    FatalErrorInFunction
                << "Hydrostatic Initialisation requires exactly 1 patch with fixed value b.c."
                << endl
                << "The face has to be normal to the gravitational acceleration direction "
                << Foam::exit(FatalError);
	}    			
	else
	{
	    if ( (mesh.Sf().boundaryField()[bdryPref][0]&g).value() > 0 ) 
            {
                 Uh_ref.value()=1.0;
            }
        else
            {
                Uh_ref.value()=-1.0;
            }
        }
        	
        ph = pRef;
					
        // We first define a linear profile in the direction parallel to g
        // This helps in the solution of the ph eqn at the first iteration
        if ( g.component(0).value() != 0.0 )
        { 
            Info << "gravitational acceleration parallel to x-axis" << endl;
            for (label i=0; i<5; i++)
            {
                p = ph;
                fvConstraints.constrain(p);
                fluid.correctThermo();
                rho = fluid.rho();
                if ( (mesh.Sf().boundaryField()[bdryPref][0]&g).value() > 0 ) 
                {
                    hRef = xMin;
                }
                else
                {
                    hRef = xMax;
                }
                ph = pRef-(x-hRef)*(0.5*min(rho)+0.5*max(rho))*mag(g);
                Info<< "min ph " << min(ph).value() << 
                       " max ph " << max(ph).value() << endl;
            }
        }
        else
        {
            if ( g.component(1).value() != 0.0 )
            { 
                Info << "gravitational acceleration parallel to y-axis" << endl;
                for (label i=0; i<5; i++)
       	        {
                    p = ph;
                    fvConstraints.constrain(p);
                    fluid.correctThermo();
                    rho = fluid.rho();
                    if ( (mesh.Sf().boundaryField()[bdryPref][0]&g).value() > 0 )
                    { 
                        hRef = yMin;
                        Info << "yMin " << yMin << endl;
                    }
                    else
                    {
                        hRef = yMax;
                        Info << "yMax " << yMax << endl;
                    }
                        ph = pRef-(y-hRef)*(0.5*min(rho)+0.5*max(rho))*mag(g);
                        Info<< "min ph " << min(ph).value() << 
                               " max ph " << max(ph).value() << endl;
                    }
                }
                else
                {
                    Info << "Gravitational acceleration parallel to z-axis" << endl;            		
                    for (label i=0; i<5; i++)
                    {
                        p = ph;
                        fvConstraints.constrain(p);
                        fluid.correctThermo();
                        rho = fluid.rho();
                        if ( (mesh.Sf().boundaryField()[bdryPref][0]&g).value() > 0 )
                        { 
                            hRef = zMin;
                        }
                        else
                        { 
                            hRef = zMax;
                        }
                    ph = pRef-(z-hRef)*(0.5*min(rho)+0.5*max(rho))*mag(g);
                    Info<< "min ph " << min(ph).value() << 
                           "max ph " << max(ph).value() << endl;
                }
            }
        }

        Info << "hRef " << hRef.value() << endl;

        // Uh is the velocity field used for the upwinding in the divergence.
        // It is a vector with only one component !=0, in order to have
        // div(Uh*ph) equal to the derivative of ph in the direction parallel
        // to that of the gravitational acceleration 
        volVectorField& Uh = regIOobject::store
        (
            new volVectorField
            (
                IOobject
                (
                    "Uh",
                    "0",
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimVelocity
            )
        );
    
        Uh = Uh_ref * (-g) / mag(g);
        Info << "Upwind velocity for Hydrostatic Initialisation " << Uh[0] << endl;
   
        // phiUh is the volumetric flux of Uh, used in the fvm discretization
        // of the divergence term to have the desired upwinding direction 
        surfaceScalarField phiUh
        (
            IOobject
            (
       	        "phiUh",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::flux(Uh)
        );

        // this is the constant part of the known term of the equation	
        const dimensionedScalar& Uhg = -Uh_ref * mag(g);
	                   
        label nCorr
        (
            pimple.dict().lookupOrDefault<label>("nHydrostaticCorrectors", 5)
        );
            
        for (label i=0; i<nCorr; i++)
        { 
            p = ph;
       	    fluid.correctThermo();
       	    rho = fluid.rho();
                         
   	    fvScalarMatrix phEqn
            (
                fvm::div(phiUh, ph) - rho*Uhg
            );
                
   	    volScalarField err_ph("err_ph", fvc::div(phiUh, ph) - rho*Uhg);

            phEqn.solve();
   	    p = ph;
   	    fvConstraints.constrain(p);
   	    ph = p;

            Info<< "min ph " << min(ph).value() <<
                   " max ph " << max(ph).value() << endl;

            Info<< "min rho " << min(rho).value() <<
                   " max rho " << max(rho).value() << endl;
   	                  	               
        } 

        fluid.correctThermo();
        rho = fluid.rho();

                
        forAll (ph.boundaryField(), bdryID)     
        {   
   		     		
            if ( ph.boundaryField()[bdryID].type() == "fixedValue" )
            {
                nFixValPatches +=1;
                const fvPatchScalarField& ph_p = ph.boundaryField()[bdryID];
    			    
                if ( min(ph_p) == max(ph_p) )
                {
                    pRef.value() = min(ph_p);
                    bdryPref = bdryID;
                    Info << "pRef " << pRef.value() << endl;
                }			
            }
        }
        
        Info << "pBdry " << ph.boundaryField()[bdryPref][0] << endl;
        
        hRef.value() = 0.0;
        Info << "hRef " << hRef.value() << endl;
        
        #include "gh.H"
    
        // the new hydrostatic pressure profile is used to update the density field 		
        p = ph;
        fluid.correctThermo();
        rho = fluid.rho();
        	
        // we initialize the field ph_rgh with the computed pressure and density
        ph_rgh = ph - rho*gh;
       
        // Find the label of the patch that shall be changed
        label patchID = mesh.boundaryMesh().findPatchID("outlet");
        
        // we change the fixed value b.c. of ph_rgh at the top face, in order to be 
        // consistent with the values of ph, rho and gh
        forAll(ph_rgh.boundaryField()[patchID], faceI)
        {
            ph_rgh.boundaryFieldRef()[patchID][faceI] = ph.boundaryField()[bdryPref][faceI] 
                - rho.boundaryField()[bdryPref][faceI] * gh.boundaryField()[bdryPref][faceI];
        }

        surfaceScalarField rhof("rhof", fvc::interpolate(rho));
        surfaceScalarField phig
        (
            "phig",-rhof*ghf*fvc::snGrad(rho)*mesh.magSf()
        );

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(ph_rgh, rho, U, phig, rhof);
        p = ph_rgh + rho*gh;
            
        fluid.correctThermo();
        rho = fluid.rho();

        const fvPatchScalarField& ph_rgh_top = ph_rgh.boundaryField()[patchID];

        Info<< "min ph_rgh top " << min(ph_rgh_top) <<
   	       "max ph_rgh top " << max(ph_rgh_top) << endl;
  	    
        // This is the same loop of the original OpenFOAM initialization
        // but now at the first iteration both ph_rgh and rho come from the
        // pressure field computed above.
        // This should help the convergence with large domains.
                
        //for (label i=0; i<nCorr; i++)
        for (label i=0; i<0; i++)
        {
            surfaceScalarField rhof("rhof", fvc::interpolate(rho));
            surfaceScalarField phig
            (
                "phig",-rhof*ghf*fvc::snGrad(rho)*mesh.magSf()
            );

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(ph_rgh, rho, U, phig, rhof);

            fvScalarMatrix ph_rghEqn
            (
                fvm::laplacian(rhof, ph_rgh) == fvc::div(phig)
            );

            ph_rghEqn.solve();

            // with respect to the original hydrostatic initialization
            // here we do not add the contribution of pRef
            // in this way we are consistent with what is done in pEqn
            p = ph_rgh + rho*gh;
            
            fluid.correctThermo();
            rho = fluid.rho();


            Info<< "Hydrostatic pressure variation "
                << (max(ph_rgh) - min(ph_rgh)).value() << endl;
            Info<< "min p " << min(p).value() <<
                   "max p " << max(p).value() << endl;
            Info<< "min rho " << min(rho).value() <<
                   "max rho " << max(rho).value() << endl;

        }
        ph_rgh.write(); 
        p.write();
        rho.write(); 
	T1.write();
	T2.write();
        p = ph_rgh + rho*gh;
        
    }
}
else
{
    p.read();
    fluid.correctThermo();
    rho = fluid.rho();
    #include "gh.H"
    ph_rgh = p - rho*gh;
}    

    
