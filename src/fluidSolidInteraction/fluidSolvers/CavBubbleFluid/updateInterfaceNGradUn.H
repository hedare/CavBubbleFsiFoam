{
    Info << "Update normal derivative of normal velocity" << endl;

    volVectorField phiU = fvc::reconstruct(phi_);

    vectorField n = mesh.boundary()[fsi.fluidPatchIndex()].nf();

    scalarField UnP =
        (n&phiU.boundaryField()[fsi.fluidPatchIndex()].patchInternalField());

    scalarField UnInt =
        phi_.boundaryField()[fsi.fluidPatchIndex()]
       /mesh.magSf().boundaryField()[fsi.fluidPatchIndex()];

    nGradUn_ =
        (UnInt - UnP)
       *mesh.deltaCoeffs().boundaryField()[fsi.fluidPatchIndex()];

    Info << "nGradUn, max: " << gMax(nGradUn_) << ", min: "
        << gMin(nGradUn_) << ", avg: " << gAverage(nGradUn_) << endl;
}
