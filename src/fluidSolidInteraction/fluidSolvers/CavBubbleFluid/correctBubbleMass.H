 if (correctBubbleMass)
{
    alpha2rho = alpha2*rho2;
    bubblem = Foam::gSum(mesh.V()*alpha2rho.field());
    bubbleV = Foam::gSum(mesh.V()*alpha2.field());
    bubblet = runTime().value();
    rho2 = rho2*bubblem0/bubblem;
    if (bubblePhase == 0)
    {
    	if (bubblet > bubblet0)
    	{
		if (bubbleV < bubbleV0 && bubblet > 1e-6)
		{
		    bubblePhase = 1;
		    rho2 = rho2*pfact;
		    alpha2rho = alpha2*rho2;
		    bubblem0 = Foam::gSum(mesh.V()*alpha2rho.field());
		}
		bubbleV0 = bubbleV;
		bubblet0 = bubblet;
	}
    }/* else if (bubblePhase == 1)
    {
        if (bubbleV > bubbleV0)
        {
            bubblePhase = 2;
        }
        bubbleV0 = bubbleV;
    }*/
    p_g = (p20 + B2) * pow(rho2/rho20 , gamma2 ) - B2;
    prho = max(alpha1*prho + alpha2*p_g, pMin); // + rho*gh, pMin);
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
}
