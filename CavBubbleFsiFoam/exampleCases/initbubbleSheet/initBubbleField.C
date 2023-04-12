/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    initbubble

Description
    Initalization utility for example case

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include <iostream>
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Info<< "Reading \n" << endl;

    dimensionedScalar epsDtB
    (
        transportProperties.lookup("epsDtB")
    );

    dimensionedVector BubbleCenter
    (
        transportProperties.lookup("BubbleCenter")
    );

    dimensionedScalar BubbleRadius
    (
        transportProperties.lookup("BubbleRadius")
    );

    dimensionedScalar epsionBubble
    (
        transportProperties.lookup("epsionBubble")
    );

    dimensionedScalar pliquid
    (
        transportProperties.lookup("pliquid")
    );

    dimensionedScalar pBubble
    (
        transportProperties.lookup("pBubble")
    );

    volScalarField alphaBubble(alpha1);

    volScalarField alphaDroplet(alpha1);

	
    forAll(alpha1, cellI)
    {
    	vector x = mesh.C()[cellI];

        scalar kb = magSqr((x[0]-BubbleCenter.value()[0])/BubbleRadius.value())+magSqr((x[1]-BubbleCenter.value()[1])/BubbleRadius.value())+magSqr((x[2]-BubbleCenter.value()[2])/BubbleRadius.value());
        //scalar kb = magSqr((x[0]-BubbleCenter.value()[0])/BubbleRadius.value())+magSqr((x[2]-BubbleCenter.value()[2])/BubbleRadius.value());
        
        alphaBubble[cellI] = (1-Foam::tanh((kb-1)/epsionBubble.value()))/2;
    }


//smear the interface to avoid the numerical instability grows in the beginning
    volScalarField alphaBubble0(alphaBubble);

    fvScalarMatrix alphaBubbleSEqn
    (
        fvm::Sp(scalar(1),alphaBubble) - fvm::laplacian(epsDtB,alphaBubble) == alphaBubble0
    );

    alphaBubbleSEqn.solve();
    
    forAll(alpha1, cellI)
    {
    	vector x = mesh.C()[cellI];

        scalar kb = magSqr((x[0]-BubbleCenter.value()[0])/BubbleRadius.value())+magSqr((x[1]-BubbleCenter.value()[1])/BubbleRadius.value())+magSqr((x[2]-BubbleCenter.value()[2])/BubbleRadius.value());
        //scalar kb = magSqr((x[0]-BubbleCenter.value()[0])/BubbleRadius.value())+magSqr((x[2]-BubbleCenter.value()[2])/BubbleRadius.value());
        
        if (kb > 2 || x[1] < -100e-6)
        {
            alphaBubble[cellI] = 0.0;
        }
    }

    p = pliquid*(1-alphaBubble) + pBubble*alphaBubble;
    
    forAll(alpha1, cellI)
    {
    	vector x = mesh.C()[cellI];
        if (x[1] < -0e-6)// || x[1] > 1e-3)
        {
            alphaBubble[cellI] = 1.0;
        }// else {
        //    alphaBubble[cellI] = 0.0;
        //}
    }
    alpha1 = max(min(1.0-alphaBubble,1.0),0.0);
    
    /*forAll(alpha1, cellI)
    {
    	vector x = mesh.C()[cellI];
        if (x[1] < 0)
        {
            alpha1[cellI] = 0;
        }
    }*/
         
    alpha1.write();

    p.write();
    
    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
