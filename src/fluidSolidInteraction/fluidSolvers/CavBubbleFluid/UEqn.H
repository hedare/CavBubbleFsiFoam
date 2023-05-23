    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rhoField*turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rhoField, U_)
      + fvm::div(rhoPhi, U_)
      - fvm::laplacian(muEff, U_)
      - (fvc::grad(U_) & fvc::grad(muEff))
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    UEqn.relax();

    if (momentumPredictor)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                fvc::interpolate(rhoField)*(g & mesh.Sf())
              + (
                    fvc::interpolate(sigma*K_)*fvc::snGrad(alpha1)
                  - fvc::snGrad(prho)
                ) * mesh.magSf()
            )
        );
    }
