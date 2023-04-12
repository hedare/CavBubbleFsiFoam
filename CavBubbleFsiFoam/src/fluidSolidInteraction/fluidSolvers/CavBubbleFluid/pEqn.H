{
    tmp<fvScalarMatrix> pEqnComp;

    if (transSonic)
    {
        pEqnComp =
            (fvm::ddt(prho) + fvm::div(phi_, prho) - fvm::Sp(fvc::div(phi_), prho));
    }
    else
    {
        pEqnComp =
            (fvm::ddt(prho) + fvc::div(phi_, prho) - fvc::Sp(fvc::div(phi_), prho));
    }


    U_ = rAU*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U_) & mesh.Sf()) + fvc::ddtPhiCorr(rAU, rhoField, U_, phi_)
    );

    phi_ = phiU +
        (
            fvc::interpolate(sigma*K_)*fvc::snGrad(alpha1)*mesh.magSf()
          + fvc::interpolate(rhoField)*(g & mesh.Sf())
        )*rAUf;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phi_)
          - fvm::laplacian(rAUf, prho)
        );

        solve
        (
            (
                max(alpha1, scalar(0))*(psi1/rho1)
              + max(alpha2, scalar(0))*(psi2/rho2)
            )
           *pEqnComp()
          + pEqnIncomp
        );

        if (nonOrth == nNonOrthCorr)
        {
            dgdt =
                (pos(alpha2)*(psi2/rho2) - pos(alpha1)*(psi1/rho1))
               *(pEqnComp & prho);
            phi_ += pEqnIncomp.flux();
        }
    }

    U_ += rAU*fvc::reconstruct((phi_ - phiU)/rAUf);
    U_.correctBoundaryConditions();

    prho.max(pMin);
    
    pmin = min(pmin,prho);
    pmax = max(pmax,prho);

    rho1 = rho10*pow(((prho+B1)/(p10+B1)),1/gamma1);
    rho2 = rho20*pow(((max(prho,pMinGas)+B2)/(p20+B2)),1/gamma2);
    rhoField = alpha1*rho1 + alpha2*rho2;
    psi1 = rho1/(prho+B1)/gamma1;
    psi2 = rho2/(max(prho,pMinGas)+B2)/gamma2;
    if (correctPsi>0)
    {
        psi1 = pow(alpha1, correctPsi)*psi1;
        psi2 = pow(alpha2, correctPsi)*psi2;
    }

    Info<< "max(U) " << max(mag(U_)).value() << endl;
    Info<< "min(p) " << min(prho).value() << endl;
}
