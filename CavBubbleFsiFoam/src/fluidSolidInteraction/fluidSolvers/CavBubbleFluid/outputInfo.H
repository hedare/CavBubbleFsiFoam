alpha2p = alpha2*prho;
if (correctBubbleMass == false)
{
    alpha2rho = alpha2*rho2;
    bubblem = Foam::gSum(mesh.V()*alpha2rho.field());
    bubbleV = Foam::gSum(mesh.V()*alpha2.field());
}
infovstime << runTime().timeName()<<","
           << bubbleV << ","
           << Foam::gSum(mesh.V()*alpha2p.field()) << ","
           << bubblem << nl;
infovstime.flush();
