{
    surfaceScalarField phic = mag(phi_/mesh.magSf());
    phic = min(cAlpha_*phic, max(phic));

    volScalarField divU = fvc::div(phi_);

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime().deltaT();
        surfaceScalarField rhoPhiSum = 0.0*rhoPhi;

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqns.H"
            rhoPhiSum += (runTime().deltaT()/totalDeltaT)*rhoPhi;
        }

        rhoPhi = rhoPhiSum;
    }
    else
    {
        #include "alphaEqns.H"
    }
    
    #include "calculateK.H"
}
