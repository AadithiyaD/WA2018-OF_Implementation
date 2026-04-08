/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "GeometricScalarField.H"
#include "WA2018.H"
#include "fvOptions.H"
#include "bound.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void WA2018<BasicTurbulenceModel>::correctNut()
{
    // Calc and init chi and fmu
    volScalarField chi(Rwa_/this->nu());

    // fmu is the damping function accounting for wall blocking
    volScalarField fmu(pow3(chi)/( pow3(chi) + pow3(Cw_) ) );
    
    // Calculate damped turbulent viscosity
    this->nut_ = fmu*Rwa_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WA2018<BasicTurbulenceModel>::sigmaR
(
    const volScalarField& f1
)
{
    // Directly return sigmaR
    return ((f1*(sigmakOmega_ - sigmakEps_)) + sigmakEps_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WA2018<BasicTurbulenceModel>::C1
(
    const volScalarField& f1
)
{   
    // Directly return C1
    return (f1*(C1kOmega_ - C1kEps_) + C1kEps_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WA2018<BasicTurbulenceModel>::F1WA
(
    const volScalarField& magS,
    const volScalarField& magW
)
{
    // calculate eta
    volScalarField eta(magS * max(scalar(1), magW/magS));
    
    // Calculate arg1 
    volScalarField omega = WA_omega(magS);
    volScalarField k = WA_k(omega);
    volScalarField arg1( ((this->nu() + Rwa_)/2) * ((sqr(eta))/max(Cmu_*k*omega,dimensionedScalar("SMALL", 
                                                                    dimensionSet(0, 2, -3, 0, 0), 
                                                                    SMALL)
                                                                     ) ) );

    // return tanh(arg1^4)
    return tanh(pow4(arg1));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WA2018<BasicTurbulenceModel>::WA_k
(
    const volScalarField& omega
)
{
    return this->nut_*omega;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WA2018<BasicTurbulenceModel>::WA_omega
(
    const volScalarField& magS
)
{    
    return magS/sqrt(Cmu_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WA2018<BasicTurbulenceModel>::WA2018
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    C1kOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1kOmega",
            this->coeffDict_,
            0.0829
        )
    ),
    C1kEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1kEps",
            this->coeffDict_,
            0.1284
        )
    ),
    sigmakEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmakEps",
            this->coeffDict_,
            1.0
        )
    ),
    sigmakOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmakOmega",
            this->coeffDict_,
            0.72
        )
    ),
    C2kOmega_(C1kOmega_ / sqr(0.41) + sigmakOmega_), // Directly init 'derived' constants
    C2kEps_(C1kEps_ / sqr(0.41) + sigmakEps_),
    Cw_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw",
            this->coeffDict_,
            8.54
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    Cm_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cm",
            this->coeffDict_,
            8
        )
    ),

    Rwa_
    (
        IOobject
        (
            IOobject::groupName("Rwa", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(Rwa_, dimensionedScalar("small", Rwa_.dimensions(), SMALL));

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WA2018<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        C1kOmega_.readIfPresent(this->coeffDict());
        C1kEps_.readIfPresent(this->coeffDict());
        sigmakEps_.readIfPresent(this->coeffDict());
        sigmakOmega_.readIfPresent(this->coeffDict());
        C2kOmega_.readIfPresent(this->coeffDict());
        C2kEps_.readIfPresent(this->coeffDict());
        Cw_.readIfPresent(this->coeffDict());
        Cmu_.readIfPresent(this->coeffDict());
        Cm_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void WA2018<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    // Create a tmp gradU object
    tmp<volTensorField> tgradU = fvc::grad(U);

    //* Calc the variable magS, which represents S from the model
    //* Take max() to avoid div by zero errors
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField magS(sqrt(S2)); 
    S2 = max(S2, dimensionedScalar("0", S2.dimensions(), SMALL));
    magS = max(magS, dimensionedScalar("0", magS.dimensions(), SMALL));

    //* Similarly, calculate W, which is the antisymmetric part of the velocity gradient
    volScalarField W2(2*magSqr(skew(tgradU())));
    volScalarField magW(sqrt(W2));

    // Clear tgradu since not needed after this point
    tgradU.clear();

    //TODO- Can all of these be refs or tmp objects to be more memory efficient ?
    // Calc f1, sigmaR,and C1
    volScalarField f1 = F1WA(magS, magW);
    volScalarField C1Var = C1(f1);
    volScalarField sigmaRVar = sigmaR(f1);

    // R equation
    //* Note for last term (i.e term before fvOptions) in the eqn - grad(R or S) gives vector,
    //* so you do magSqr and convert to scalar to make it compatible with Cm and S2
    tmp<fvScalarMatrix> RwaEqn
    (
          fvm::ddt(Rwa_)
        + fvm::div(alphaRhoPhi, Rwa_)
        - fvm::laplacian((sigmaRVar*Rwa_) + this->nu(), Rwa_)
        ==
          C1Var*fvm::Sp(magS, Rwa_)
        + fvm::Sp((f1* C2kOmega_/magS)*(fvc::grad(Rwa_) & fvc::grad(magS)),Rwa_)
        - (1.0-f1)*min(C2kEps_*(Rwa_*Rwa_)*((magSqr(fvc::grad(magS))/S2)),
                                         Cm_*magSqr(fvc::grad(Rwa_)))
        + fvOptions(alpha, rho, Rwa_)
    );

    RwaEqn.ref().relax();
    fvOptions.constrain(RwaEqn.ref());

    solve(RwaEqn);
    fvOptions.correct(Rwa_);
    bound(Rwa_, dimensionedScalar("small", Rwa_.dimensions(), SMALL));
    Rwa_.correctBoundaryConditions();

    // Update turbulent viscosity
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
