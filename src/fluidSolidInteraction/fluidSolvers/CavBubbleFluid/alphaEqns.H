{
    word alphaScheme("div(phi_,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir = phic*nHatf_;

    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {
        volScalarField::DimensionedInternalField Sp
        (
            IOobject
            (
                "Sp",
                runTime().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("Sp", dgdt.dimensions(), 0.0)
        );

        volScalarField::DimensionedInternalField Su
        (
            IOobject
            (
                "Su",
                runTime().timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            divU*min(alpha1, scalar(1))
        );

        forAll(dgdt, celli)
        {
            if (dgdt[celli] > 0.0 && alpha1[celli] > 0.0)
            {
                Sp[celli] -= dgdt[celli]*alpha1[celli];
                Su[celli] += dgdt[celli]*alpha1[celli];
            }
            else if (dgdt[celli] < 0.0 && alpha1[celli] < 1.0)
            {
                Sp[celli] += dgdt[celli]*(1.0 - alpha1[celli]);
            }
        }

        surfaceScalarField phiAlpha1 =
            fvc::flux
            (
                phi_,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            );

        MULES::explicitSolve(geometricOneField(), alpha1, phi_, phiAlpha1, Sp, Su, 1, 0);

        surfaceScalarField rho1f = fvc::interpolate(rho1);
        surfaceScalarField rho2f = fvc::interpolate(rho2);
        rhoPhi = phiAlpha1*(rho1f - rho2f) + phi_*rho2f;
        
        if (limitAlpha > 0.0)
        {
            if (limitAlpha < 1.0)
            {
            scalar limitAlpha2 = max(1 - (1-limitAlpha)*runTime().deltaT().value()/3e-9,limitAlpha);//weaker for smaller dt
                forAll(alpha1, celli)
                {
                    if (alpha1[celli] > limitAlpha2) {alpha1[celli] = 1.0;}
                    else if (alpha1[celli] < 1 - limitAlpha2) {alpha1[celli] = 0.0;}
                }
            }
            else {alpha1 = max(min(alpha1, 1.0), 0.0);}
        }
        alpha2 = scalar(1) - alpha1;
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Min(alpha2) = " << min(alpha2).value()
        << endl;
}
