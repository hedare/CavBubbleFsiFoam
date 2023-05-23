dimensionedScalar Tinit
(
    transportProperties_.lookup("Tinit")
);

if (secondBubbleInfo == 0.0 && runTime().value() >= Tinit.value())
{
    secondBubbleInfo = runTime().value();
} else if (secondBubbleInfo < runTime().value() && secondBubbleInfo > 0.0)
{
    secondBubble = 0.0;
}

if (runTime().value() >= Tinit.value() && secondBubbleInfo <= runTime().value() && secondBubble > 0.0)
{

    Info << endl << "Initiating second bubble" << endl << endl;

    dimensionedScalar epsDtB
    (
        transportProperties_.lookup("epsDtB")
    );

    dimensionedVector BubbleCenter
    (
        transportProperties_.lookup("BubbleCenter2")
    );

    dimensionedScalar BubbleRadius
    (
        transportProperties_.lookup("BubbleRadius2")
    );

    dimensionedScalar epsionBubble
    (
        transportProperties_.lookup("epsionBubble")
    );

    dimensionedScalar pliquid
    (
        transportProperties_.lookup("pliquid")
    );

    dimensionedScalar pBubble
    (
        transportProperties_.lookup("pBubble")
    );
    
    volScalarField alphaBubble(alpha1);

    forAll(alpha1, cellI)
    {
    	vector x = mesh.C()[cellI];

        scalar kb = magSqr((x[0]-BubbleCenter.value()[0])/BubbleRadius.value())+magSqr((x[1]-BubbleCenter.value()[1])/BubbleRadius.value())+magSqr((x[2]-BubbleCenter.value()[2])/BubbleRadius.value());
        alphaBubble[cellI] = (1-Foam::tanh((kb-1)/epsionBubble.value()))/2;
    }

//smear the interface to avoid the numerical instability grows in the beginning
    alphaBubble = 1 - alphaBubble;
    volScalarField alphaBubble0(alphaBubble);

    fvScalarMatrix alphaBubbleSEqn
    (
        fvm::Sp(scalar(1),alphaBubble) - fvm::laplacian(epsDtB,alphaBubble) == alphaBubble0
    );

    alphaBubbleSEqn.solve();
    alphaBubble = 1 - alphaBubble;
    
    forAll(alpha1, cellI)
    {
    	vector x = mesh.C()[cellI];

        scalar kb = magSqr((x[0]-BubbleCenter.value()[0])/BubbleRadius.value())+magSqr((x[1]-BubbleCenter.value()[1])/BubbleRadius.value())+magSqr((x[2]-BubbleCenter.value()[2])/BubbleRadius.value());
        
        if (kb > 2)
        {
            alphaBubble[cellI] = 0.0;
        }
    }

    //alpha1 = max(min(1.0-alphaBubble,1.0),0.0);
    alpha2 = max(min(alpha2 + alphaBubble,1.0),0.0);
    alpha1 = 1.0 - alpha2;
    prho = prho + pBubble*alphaBubble;
    p_g = prho;
    rho1 = rho10*pow(((max(prho,pMin)+B1)/(p10+B1)),1/gamma1);
    rho2 = rho20*pow(((max(prho,pMinGas)+B2)/(p20+B2)),1/gamma2);
    rhoField = alpha1*rho1 + alpha2*rho2;
    p_ = prho/rho_;
    gradp_ = fvc::grad(p_);
    psi1 = rho1/(prho+B1)/gamma1;
    psi2 = rho2/(max(prho,pMinGas)+B2)/gamma2;
    if (correctPsi>0)
    {
        psi1 = pow(alpha1, correctPsi)*psi1;
        psi2 = pow(alpha2, correctPsi)*psi2;
    }
    if (correctBubbleMass)
    {
        bubbleV = Foam::gSum(mesh.V()*alpha2.field());
    }
}
