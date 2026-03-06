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
    volScalarField fmu(pow3(chi)/( pow3(chi) + pow3(Cm_) ) );

    // Calculate turbulent viscosity
    this->nut_ = this->rho_*fmu*Rwa_;
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
    return ((f1*(C1kOmega_ - C1kEps_)) + C1kEps_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WA2018<BasicTurbulenceModel>::C1
(
    const volScalarField& f1
)
{   
    // Directly return C1
    return (f1*(C1kOmega_ + C1kEps_) + C1kEps_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WA2018<BasicTurbulenceModel>::F1WA()
{
    //! Would it be better if i put all of these as refs instead of direct init? Would that be faster / more memory efficient?
    // First calculate eta
    // Calculate grad U and get S
    tmp<volTensorField> tgradU = fvc::grad(this->U_);

    // We now get the variable magS, which represents S from the model
    volScalarField S2(2*magSqr(devSymm(tgradU())));

    //* using max (to avoid div by 0 errors) led to some issues at runtime
    //* Therefore, use bound instead of max
    volScalarField magS(sqrt(S2)); 
    bound(magS, dimensionedScalar("0", magS.dimensions(), SMALL));

    // Similarly, calculate W, which is the antisymmetric part of the velocity grad
    volScalarField W2(2*magSqr(skew(tgradU())));
    volScalarField magW(sqrt(W2));

    // Now we can calculate eta
    volScalarField eta(magS * max(scalar(1), magW/magS));
    
    // Calculate arg1 
    volScalarField arg1( ((this->nu() + Rwa_)/2) * ((sqr(eta))/(Cmu_*k_*omega_) ) );

    // return tanh(arg1^4)
    return tanh(pow4(arg1));
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

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
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
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    // kMin and omegaMin are defined in RASModel.H . So if I wanted to do similar for R, would
    // I have to modify that file?
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
    // const volScalarField& nut = this->nut_; //! WTF if nut is not used in the eqn then what is its point?

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    // Create a tmp gradU object
    tmp<volTensorField> tgradU = fvc::grad(U);

    //! Can all of these be refs or tmp objects to be more memory efficient ?
    // Calc f1, sigmaR,and C1
    volScalarField f1Var = F1WA();
    volScalarField C1Var = C1(f1Var);
    volScalarField sigmaRVar = sigmaR(f1Var);

    // Calc S (Remember, magS => S)
    volScalarField S2(2*magSqr(devSymm(tgradU()))); 
    volScalarField magS(max(sqrt(S2), dimensionedScalar("small", dimensionSet(0, 0, -1, 0, 0), SMALL)));
    // volScalarField magS(sqrt(S2)); 
    // bound(magS, dimensionedScalar("0", magS.dimensions(), SMALL));

    // Clear tgradu since not needed after this point
    //! Can tgradU be safely cleared here?
    // tgradU.clear();

    //! Does R need to be updated at the wall? 
    // Update R at the wall
    // Rwa_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    // Rwa_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    // DIAGNOSTIC: Print field statistics BEFORE solving
    Info << nl << "=== WA2018 correct() START ===" << nl;
    Info << "Rwa   min/max: " << gMin(Rwa_) << " / " << gMax(Rwa_) << endl;
    Info << "k     min/max: " << gMin(k_) << " / " << gMax(k_) << endl;
    Info << "omega min/max: " << gMin(omega_) << " / " << gMax(omega_) << endl;
    Info << "nut   min/max: " << gMin(this->nut_) << " / " << gMax(this->nut_) << endl;

    // R equation
    // Note for last term (i.e term before fvOptions) in the eqn -  
    // grad(R or S) gives vector, so you do magSqr to make it compatible with Cm and sqr(magS)
    tmp<fvScalarMatrix> RwaEqn
    (
          fvm::ddt(Rwa_)
        + fvm::div(alphaRhoPhi, Rwa_)
        - fvm::laplacian((sigmaRVar*Rwa_) + this->nu(), Rwa_)
        ==
        (C1Var*Rwa_*magS)
        + (f1Var* C2kOmega_*fvm::Sp((fvc::grad(Rwa_) & fvc::grad(magS))/magS,Rwa_) )
        - (1-f1Var)* min(C2kEps_*sqr(Rwa_)*((magSqr(fvc::grad(magS))/sqr(magS))), Cm_*magSqr(fvc::grad(Rwa_)) ) 
        + fvOptions(alpha, rho, Rwa_)
    );

    RwaEqn.ref().relax();
    fvOptions.constrain(RwaEqn.ref());

    //! R is not specified at the boundary, so is this necesary?
    // RwaEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    
    solve(RwaEqn);
    fvOptions.correct(Rwa_);
    bound(Rwa_, dimensionedScalar("small", Rwa_.dimensions(), SMALL));
    Rwa_.correctBoundaryConditions();

    // Update turbulent viscosity
    correctNut();

    // DIAGNOSTIC: Print AFTER solving
    Info << nl << "=== AFTER SOLVING ===" << nl;
    Info << "Rwa   min/max: " << gMin(Rwa_) << " / " << gMax(Rwa_) << endl;
    Info << "k     min/max: " << gMin(k_) << " / " << gMax(k_) << endl;
    Info << "omega min/max: " << gMin(omega_) << " / " << gMax(omega_) << endl;
    Info << "nut   min/max: " << gMin(this->nut_) << " / " << gMax(this->nut_) << endl;
    Info << "=== WA2018 correct() END ===" << nl << endl;

    //! Should f1, C1, sigmaR, S also be cleared at the end?
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
