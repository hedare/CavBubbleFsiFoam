infovstime.open ("info.csv", std::ios::out|std::ios::trunc);
bubbleV0 = Foam::gSum(mesh.V()*alpha2.field());
bubblem0 = Foam::gSum(mesh.V()*alpha2rho.field());
infovstime << "%Time,bubbleVolume,bubblepressure,bubblemass" << nl;
infovstime << runTime().timeName()<<","
           << bubbleV0 << ","
           << Foam::gSum(mesh.V()*alpha2p.field())/bubbleV0 << ","
           << bubblem0 << nl;
infovstime.flush();
