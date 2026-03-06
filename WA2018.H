/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RASModels::WA2018

Group
    grpRASTurbulence

Description
    Standard k-epsilon turbulence model for incompressible and compressible
    flows including rapid distortion theory (RDT) based compression term.

    Reference:
    \verbatim
        Standard model:
            Launder, B. E., & Spalding, D. B. (1972).
            Lectures in mathematical models of turbulence.

            Launder, B. E., & Spalding, D. B. (1974).
            The numerical computation of turbulent flows.
            Computer methods in applied mechanics and engineering,
            3(2), 269-289.

        For the RDT-based compression term:
            El Tahry, S. H. (1983).
            k-epsilon equation for compressible reciprocating engine flows.
            Journal of Energy, 7(4), 345-353.
    \endverbatim

    The default model coefficients are
    \verbatim
        WA2018Coeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            C3          -0.33;
            sigmak      1.0;
            sigmaEps    1.3;
        }
    \endverbatim

SourceFiles
    WA2018.C


\*---------------------------------------------------------------------------*/

#ifndef WA2018_H
#define WA2018_H

#include "RASModel.H"
#include "dimensionedScalarFwd.H"
#include "eddyViscosity.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class WA2018 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WA2018
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        WA2018(const WA2018&) = delete;

        //- No copy assignment
        void operator=(const WA2018&) = delete;


protected:

    // Protected data

        // Model coefficients
            dimensionedScalar C1kOmega_;
            dimensionedScalar C1kEps_;
            dimensionedScalar sigmakEps_;
            dimensionedScalar sigmakOmega_;
            dimensionedScalar C2kOmega_;
            dimensionedScalar C2kEps_;
            dimensionedScalar Cw_;
            dimensionedScalar Cmu_;
            dimensionedScalar Cm_;
            
        // Fields

            volScalarField k_;
            volScalarField omega_;
            volScalarField Rwa_; // Rwa_ => R Wray-Aggarwal, corresponds to R in the model



    // Protected Member Functions

        virtual void correctNut();

        //! For some reason, including this func in the fvMatrix RwaEqn casused an unknown symbol error at runtime
        // virtual tmp<fvScalarMatrix> RwaSource() const;

        // Coressponds to f_{1} param
        // arg1IP => arg1 input
        tmp<volScalarField> F1WA();

        tmp<volScalarField> arg1();

        tmp<volScalarField> C1
        (
            const volScalarField& f1
        );

        tmp<volScalarField> sigmaR
        (
            const volScalarField& f1
        );

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("WA2018");


    // Constructors

        //- Construct from components
        WA2018
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~WA2018() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        // - Return R value 
        // I use R_wa , since apparently R() is an openfoam reserved function that calcs reynolds stress
        // This function is not used in field calculation anymore (1-03-26)
        tmp<volScalarField> R_wa() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "R",
                    (k_ / omega_) 
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "WA2018.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
