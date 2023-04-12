/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "CavBubbleFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#include "fluidSolidInterface.H"
#include "fixedGradientFvPatchFields.H"

#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

#include "elasticSlipWallVelocityFvPatchVectorField.H"
#include "elasticWallVelocityFvPatchVectorField.H"
#include "elasticWallPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidSolvers
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(CavBubbleFluid, 0);
addToRunTimeSelectionTable(fluidSolver, CavBubbleFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CavBubbleFluid::CavBubbleFluid(const fvMesh& mesh)
:
    fluidSolver(this->typeName, mesh),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            runTime().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    B1(transportProperties_.subDict("phase1").lookup("B")),
    rho10(transportProperties_.subDict("phase1").lookup("rho0")),
    p10(transportProperties_.subDict("phase1").lookup("p0")),
    gamma1(transportProperties_.subDict("phase1").lookup("gamma")),
    mu1(transportProperties_.subDict("phase1").lookup("mu")),
    B2(transportProperties_.subDict("phase2").lookup("B")),
    rho20(transportProperties_.subDict("phase2").lookup("rho0")),
    p20(transportProperties_.subDict("phase2").lookup("p0")),
    gamma2(transportProperties_.subDict("phase2").lookup("gamma")),
    mu2(transportProperties_.subDict("phase2").lookup("mu")),
    pMin(transportProperties_.lookup("pMin")),
    pfact(transportProperties_.lookup("pfact")),
    sigma(transportProperties_.lookup("sigma")),
    nu_(mu1/rho10),
    rho_(rho10),
    alpha1
    (
        IOobject
        (
            "alpha1",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alpha2(1 - alpha1),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ),
    cAlpha_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PISO").lookup("cAlpha")
        )
    ),
    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1.time().timeName(),
            alpha1.mesh()
        ),
        alpha1.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),
    K_
    (
        IOobject
        (
            "K",
            alpha1.time().timeName(),
            alpha1.mesh()
        ),
        alpha1.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    ),
    p_
    (
        IOobject
        (
            "p",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    prho(p_*rho_),
    U_
    (
        IOobject
        (
            "U",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho1(rho10*pow((prho+B1)/(p10+B1),1/gamma1)),
    rho2(rho20*pow((prho+B2)/(p20+B2),1/gamma2)),
    rhoField
    (
        IOobject
        (
            "rho",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + alpha2*rho2
    ),
    gradp_(fvc::grad(p_)),
    gradU_(fvc::grad(U_)),
    phi_
    (
        IOobject
        (
            "phi",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_) & mesh.Sf()
    ),
    rhoPhi(fvc::interpolate(rhoField)*phi_),
    dgdt(pos(alpha2)*fvc::div(phi_)/max(alpha2, scalar(0.0001))),
    g
    (
        IOobject
        (
            "g",
            runTime().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    )
{
    #include "calculateK.H"
    mesh.schemesDict().setFluxRequired(p_.name());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const volVectorField& CavBubbleFluid::U() const
{
    return U_;
}


const volScalarField& CavBubbleFluid::p() const
{
    return p_;
}


//- Patch viscous force (N/m2)
tmp<vectorField> CavBubbleFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() = mu1.value()*U().boundaryField()[patchID].snGrad();

    vectorField n = mesh().boundary()[patchID].nf();
    tvF() -= n*(n&tvF());

    return tvF;
}

//- Patch pressure force (N/m2)
tmp<scalarField> CavBubbleFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = prho.boundaryField()[patchID];

    return tpF;
}

//- Patch viscous force (N/m2)
tmp<vectorField> CavBubbleFluid::faceZoneViscousForce
(
    const label zoneID,
    const label patchID
) const
{
    vectorField pVF = patchViscousForce(patchID);

    tmp<vectorField> tvF
    (
        new vectorField(mesh().faceZones()[zoneID].size(), vector::zero)
    );
    vectorField& vF = tvF();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(pVF, i)
    {
        vF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
            pVF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(vF, sumOp<vectorField>());


    return tvF;
}

//- Patch pressure force (N/m2)
tmp<scalarField> CavBubbleFluid::faceZonePressureForce
(
    const label zoneID,
    const label patchID
) const
{
    scalarField pPF = patchPressureForce(patchID);

    tmp<scalarField> tpF
    (
        new scalarField(mesh().faceZones()[zoneID].size(), 0)
    );
    scalarField& pF = tpF();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(pPF, i)
    {
        pF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
            pPF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(pF, sumOp<scalarField>());

    return tpF;
}

tmp<scalarField> CavBubbleFluid::faceZoneMuEff
(
    const label zoneID,
    const label patchID
) const
{
    tmp<scalarField> tMuEff
    (
        new scalarField
        (
            mesh().faceZones()[zoneID].size(),
            mu1.value()
        )
    );

    return tMuEff;
}

void CavBubbleFluid::evolve()
{
    twoPhaseMixture twoPhaseProperties(U_, phi_);
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U_, phi_, twoPhaseProperties)
    );
    
    Info << "Evolving fluid solver: " << this->type() << endl;

    const fvMesh& mesh = fluidSolver::mesh();
    
    // Read PISO controls
    int nCorr(readInt(fluidProperties().lookup("nCorrectors")));

    int nNonOrthCorr =
        readInt(fluidProperties().lookup("nNonOrthogonalCorrectors"));

    int nOuterCorr =
        readInt(fluidProperties().lookup("nOuterCorrectors"));

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p_, fluidProperties(), pRefCell, pRefValue);

    for (int oCorr = 0; oCorr < nOuterCorr; oCorr++)
    {
        #include "alphaEqnsSubCycle.H"
        
        if(mesh.moving())
        {
            // Make the fluxes relative
            phi_ -= fvc::meshPhi(U_);
        }

        rhoPhi = fvc::interpolate(rhoField)*phi_;
        solve(fvm::ddt(rhoField) + fvc::div(rhoPhi));

        // CourantNo
        {
          scalar CoNum = 0.0;
          scalar meanCoNum = 0.0;
          scalar velMag = 0.0;

          if (mesh.nInternalFaces())
          {
              surfaceScalarField SfUfbyDelta =
                  mesh.surfaceInterpolation::deltaCoeffs()*mag(phi_);

              CoNum = max(SfUfbyDelta/mesh.magSf())
                .value()*runTime().deltaT().value();

              meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf()))
                .value()*runTime().deltaT().value();

              velMag = max(mag(phi_)/mesh.magSf()).value();
          }

          Info<< "Courant Number mean: " << meanCoNum
              << " max: " << CoNum
              << " velocity magnitude: " << velMag << endl;
        }

        //#include "UEqn.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(U_)
          + fvm::div(phi_, U_)
          - fvm::laplacian(nu_, U_)
        );

        solve(UEqn == -gradp_);

        // --- PISO loop

        volScalarField rAU = 1.0/UEqn.A();
        surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

        for (int corr=0; corr<nCorr; corr++)
        {
            
            U_ = rAU*UEqn.H();
            phi_ = (fvc::interpolate(U_) & mesh.Sf());
//             + fvc::ddtPhiCorr(rUA, U_, phi_);

#           include "updateRobinFsiInterface.H"

            adjustPhi(phi_, U_, p_);

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                //#include "pEqn.H"
                
                fvScalarMatrix pEqn
                (
                    fvm::laplacian
                    (
                        rAUf, p_, "laplacian((1|A(U)),p)"
                    )
                 == fvc::div(phi_)
                   // fvm::laplacian(rAUf, p_) == fvc::div(phi_)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                gradp_ = fvc::grad(p_);
                prho = p_*rho_;

                if (nonOrth == nNonOrthCorr)
                {
                    phi_ -= pEqn.flux();
                }
            }

            // Continuity error
            {
                volScalarField contErr = fvc::div(phi_);

                scalar sumLocalContErr = runTime().deltaT().value()*
                    mag(contErr)().weightedAverage(mesh.V()).value();

                scalar globalContErr = runTime().deltaT().value()*
                    contErr.weightedAverage(mesh.V()).value();

                Info<< "time step continuity errors : sum local = "
                    << sumLocalContErr << ", global = "
                    << globalContErr << endl;
            }

            U_ -= rAU*gradp_;
            U_.correctBoundaryConditions();

            gradU_ = fvc::grad(U_);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolvers
} // End namespace Foam

// ************************************************************************* //
