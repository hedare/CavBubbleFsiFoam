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
    nCorr(readInt(fluidProperties().lookup("nCorrectors"))),
    nNonOrthCorr(readInt(fluidProperties().lookup("nNonOrthogonalCorrectors"))),
    nOuterCorr(readInt(fluidProperties().lookup("nOuterCorrectors"))),
    correctBubbleMass(mesh.solutionDict().subDict("PISO").lookupOrDefault<Switch>("correctBubbleMass", false)),
    cavModel(mesh.solutionDict().subDict("PISO").lookupOrDefault<Switch>("cavModel", false)),
    transSonic(mesh.solutionDict().subDict("PISO").lookupOrDefault<Switch>("transSonic", false)),
    momentumPredictor(mesh.solutionDict().subDict("PISO").lookupOrDefault<Switch>("momentumPredictor", false)),
    limitAlpha(mesh.solutionDict().subDict("PISO").lookupOrDefault<scalar>("limitAlpha", 0.0)),
    secondBubble(mesh.solutionDict().subDict("PISO").lookupOrDefault<scalar>("secondBubble", 0.0)),
    secondBubbleInfo(0.0),
    correctPsi(mesh.solutionDict().subDict("PISO").lookupOrDefault<scalar>("correctPsi", 0.0)),
    nAlphaCorr(readLabel(mesh.solutionDict().subDict("PISO").lookup("nAlphaCorr"))),
    nAlphaSubCycles(readLabel(mesh.solutionDict().subDict("PISO").lookup("nAlphaSubCycles"))),
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
    pMinGas(transportProperties_.lookup("pMinGas")),
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
    prho
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
    pmin
    (
        IOobject
        (
            "pmin",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        prho
    ),
    pmax
    (
        IOobject
        (
            "pmax",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        prho
    ),
    p_
    (
        IOobject
        (
            "p",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        prho/rho_
    ),
    p_g
    (
        IOobject
        (
            "p_g",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        prho
    ),
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
    rho1
    (
        IOobject
        (
            "rho1",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho10*pow((prho+B1)/(p10+B1),1/gamma1)
    ),
    rho2
    (
        IOobject
        (
            "rho2",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho20*pow((prho+B2)/(p20+B2),1/gamma2)
    ),
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
    psi1
    (
        IOobject
        (
            "psi1",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1/(prho+B1)/gamma1
    ),
    psi2
    (
        IOobject
        (
            "psi2",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho2/(prho+B2)/gamma2
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
    ),
    alpha2p(alpha2*prho),
    alpha2rho(alpha2*rho2),
    bubblem0(0),
    bubblem(0),
    bubbleV0(0),
    bubbleV(0),
    bubblet0(0),
    bubblet(0),
    bubblePhase(0)
{
    #include "calculateK.H"
    mesh.schemesDict().setFluxRequired(prho.name());
    #include "createOutputInfo.H"
    if (correctPsi>0)
    {
        psi1 = pow(alpha1, correctPsi)*psi1;
        psi2 = pow(alpha2, correctPsi)*psi2;
    }
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

    tvF() = (mu1.value()*alpha1.boundaryField()[patchID] + mu2.value()*alpha2.boundaryField()[patchID])*U().boundaryField()[patchID].snGrad();

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
            0
        )
    );
    tMuEff = mu1.value()*alpha1.boundaryField()[patchID] + mu2.value()*alpha2.boundaryField()[patchID];

    return tMuEff;
}

void CavBubbleFluid::evolve()
{
    Info << "Evolving fluid solver: " << this->type() << endl;

    twoPhaseMixture twoPhaseProperties(U_, phi_);
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U_, phi_, twoPhaseProperties)
    );
              
    const fvMesh& mesh = fluidSolver::mesh();

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

        #include "UEqn.H"
        
        volScalarField rAU = 1.0/UEqn.A();
        surfaceScalarField rAUf = fvc::interpolate(rAU);

        // --- PISO loop
        for (int corr=0; corr<nCorr; corr++)
        {

#           include "updateRobinFsiInterface.H"

            //adjustPhi(phi_, U_, p_);
            
            #include "pEqn.H"

            p_ = prho/rho_;
            gradp_ = fvc::grad(p_);

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

            gradU_ = fvc::grad(U_);
        }
        #include "correctBubbleMass.H"
    }
    rhoField = alpha1*rho1 + alpha2*rho2;

    turbulence->correct();

    #include "outputInfo.H"
    
    if (secondBubble > 0.0)
    {
        #include "initBubble.H"
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolvers
} // End namespace Foam

// ************************************************************************* //
