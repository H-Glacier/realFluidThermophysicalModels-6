/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "PengRobinson.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);
    return p/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::H(scalar p, scalar T) const
{
    scalar B = b_*p/(RR*T); 
    if(B <=0 ) { B = 1e-16; } 

    const scalar aAlpha = coef1_ - coef2_*sqrt(T) + coef3_*T;
    const scalar daAlpha = -coef2_/(2*sqrt(T))+coef3_;

    const scalar Z = this->Z(p, T);
    if (B == -Z){
    return (RR*T*(Z - 1))/this->W();
    }
    else 
    {
        return (RR*T*(Z - 1) +
               (T*daAlpha-aAlpha)/(2*sqrt(2.0)*b_) *
               log((Z + 2.414*B)/(Z-0.414*B)))/this->W();
    }
}


template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::Cp(scalar p, scalar T) const
{
    const scalar aAlpha = coef1_ - coef2_*sqrt(T) + coef3_*T;
    const scalar daAlpha = -coef2_/(2*sqrt(T))+coef3_;
    const scalar ddaAlpha = coef2_/(4*T*sqrt(T));
    const scalar A = aAlpha*p/sqr(RR*T);
    scalar B = b_*p/(RR*T); 
    if(B <=0 ) { B = 1e-16; }

    const scalar Z = this->Z(p, T);

    // const scalar M = (sqr(Z) + B*Z)/(Z - B);
    const scalar M = (sqr(Z) + 2*B*Z-sqr(B))/(Z - B); //bug finding by Nicola
    const scalar N = daAlpha*B/(b_*RR);
 
    return
    (
        (T/(2*sqrt(2.0)*b_))*ddaAlpha*log((Z + 2.414*B)/(Z-0.414*B))
      + RR*sqr(M - N)/(sqr(M) - 2*A*(Z + B))
      - RR
    )/this->W();
   
}


// This function has not been validated yet
template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::S
(
    scalar p,
    scalar T
) const
{
    const scalar daAlpha = -coef2_/(2*sqrt(T))+coef3_;
    scalar B = b_*p/(RR*T); 
    if(B <=0 ) { B = 1e-16; }

    const scalar Z = this->Z(p, T);

    return
        (
        RR*log(Z - B)
        -daAlpha/(2*sqrt(2.0)*b_)*log((Z + 2.414*B)/(Z-0.414*B))
        )/this->W();
}


template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);

    return 1.0/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    if
    (
        hasMixtureState_
     && mixingB_
     && mixingCoef1_
     && mixingCoef2_
     && mixingCoef3_
     && mixtureX_.size()
    )
    {
        const scalar relTol = 1e-9;

        if
        (
            hasCachedZ_
         && mag(p - cachedP_) <= relTol*(1.0 + mag(p))
         && mag(T - cachedT_) <= relTol*(1.0 + mag(T))
        )
        {
            return cachedZ_;
        }

        const scalar root = solveWithMixtureState(p, T);

        cachedP_ = p;
        cachedT_ = T;
        cachedZ_ = root;
        hasCachedZ_ = true;

        if (lastZPtr_)
        {
            *lastZPtr_ = root;
            if (lastZValidPtr_)
            {
                *lastZValidPtr_ = true;
            }
        }

        return root;
    }

    const scalar aAlpha = coef1_ - coef2_*sqrt(T) + coef3_*T;

    const scalar A = aAlpha*p/sqr(RR*T);
    const scalar B = b_*p/(RR*T);

    const scalar a2 = - (1.0 - B);
    const scalar a1 = A - 2*B - 3*sqr(B);
    const scalar a0 = -(A*B - sqr(B) - pow(B, 3.0));

    const scalar Q = (3*a1 - a2*a2)/9.0;
    const scalar Rl = (9*a2*a1 - 27*a0 - 2*a2*a2*a2)/54.0;

    const scalar Q3 = Q*Q*Q;
    const scalar D = Q3 + Rl*Rl;

    scalar root = -1;

    if (D <= 0)
    {
        const scalar th = ::acos(Rl/sqrt(-Q3));
        const scalar qm = 2*sqrt(-Q);
        const scalar r1 = qm*cos(th/3.0) - a2/3.0;
        const scalar r2 =
            qm*cos((th + 2*constant::mathematical::pi)/3.0) - a2/3.0;
        const scalar r3 =
            qm*cos((th + 4*constant::mathematical::pi)/3.0) - a2/3.0;

        root = max(r1, max(r2, r3));
    }
    else
    {
        // One root is real
        const scalar D05 = sqrt(D);
        scalar S = 0;
    
        if (Rl + D05 < 0)
        {
            S = -pow(mag(Rl + D05), 1.0/3.0);
        }
        else
        {
            S = pow(Rl + D05, 1.0/3.0);
        }
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0/3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0/3.0);
        }

        root = S + Tl - a2/3.0;
    }
    
    return root;
}


template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::solveWithMixtureState
(
    scalar p,
    scalar T
) const
{
    const scalar TSafe = max(T, SMALL);
    const scalar sqrtT = sqrt(TSafe);
    const scalar aAlpha = coef1_ - coef2_*sqrtT + coef3_*TSafe;

    const scalar A = aAlpha*p/sqr(RR*TSafe);
    const scalar B = b_*p/(RR*TSafe);

    List<scalar> roots;
    cubicRoots(A, B, roots);

    if (!roots.size())
    {
        const scalar a2 = -(1.0 - B);
        const scalar a1 = A - 2*B - 3*sqr(B);
        const scalar a0 = -(A*B - sqr(B) - pow(B, 3.0));

        const scalar Q = (3*a1 - a2*a2)/9.0;
        const scalar Rl = (9*a2*a1 - 27*a0 - 2*a2*a2*a2)/54.0;

        const scalar Q3 = Q*Q*Q;
        const scalar D = Q3 + Rl*Rl;

        if (D <= 0)
        {
            const scalar ratio = max(min(Rl/sqrt(-Q3), scalar(1)), scalar(-1));
            const scalar th = ::acos(ratio);
            const scalar qm = 2*sqrt(-Q);
            const scalar r1 = qm*cos(th/3.0) - a2/3.0;
            const scalar r2 =
                qm*cos((th + 2*constant::mathematical::pi)/3.0) - a2/3.0;
            const scalar r3 =
                qm*cos((th + 4*constant::mathematical::pi)/3.0) - a2/3.0;

            return max(r1, max(r2, r3));
        }
        else
        {
            const scalar D05 = sqrt(D);
            scalar S = 0;

            if (Rl + D05 < 0)
            {
                S = -pow(mag(Rl + D05), 1.0/3.0);
            }
            else
            {
                S = pow(Rl + D05, 1.0/3.0);
            }
            scalar Tl = 0;
            if (D05 > Rl)
            {
                Tl = -pow(mag(Rl - D05), 1.0/3.0);
            }
            else
            {
                Tl = pow(Rl - D05, 1.0/3.0);
            }

            return S + Tl - a2/3.0;
        }
    }

    if (roots.size() == 1)
    {
        return roots[0];
    }

    scalar selected = roots[roots.size()-1];

    if (lastZPtr_ && lastZValidPtr_ && *lastZValidPtr_)
    {
        scalar minDiff = great;
        forAll(roots, i)
        {
            const scalar diff = mag(roots[i] - *lastZPtr_);
            if (diff < minDiff)
            {
                minDiff = diff;
                selected = roots[i];
            }
        }

        if (minDiff < 1e-4)
        {
            return selected;
        }
    }

    scalar bestG = great;
    bool found = false;

    forAll(roots, i)
    {
        const scalar g = residualGibbsForRoot(roots[i], A, B, p, TSafe);

        if (g < bestG)
        {
            bestG = g;
            selected = roots[i];
            found = true;
        }
    }

    if (!found)
    {
        // fallback: pick root closest to previous if available, else largest
        if (lastZPtr_ && lastZValidPtr_ && *lastZValidPtr_)
        {
            scalar minDiff = great;
            forAll(roots, i)
            {
                const scalar diff = mag(roots[i] - *lastZPtr_);
                if (diff < minDiff)
                {
                    minDiff = diff;
                    selected = roots[i];
                }
            }
        }
        else
        {
            selected = roots[roots.size()-1];
        }
    }

    return selected;
}


template<class Specie>
inline void Foam::PengRobinson<Specie>::cubicRoots
(
    scalar A,
    scalar B,
    List<scalar>& roots
) const
{
    roots.setSize(0);

    const scalar a2 = -(1.0 - B);
    const scalar a1 = A - 2*B - 3*sqr(B);
    const scalar a0 = -(A*B - sqr(B) - pow(B, 3.0));

    const scalar Q = (3*a1 - a2*a2)/9.0;
    const scalar Rl = (9*a2*a1 - 27*a0 - 2*a2*a2*a2)/54.0;

    const scalar Q3 = Q*Q*Q;
    const scalar D = Q3 + Rl*Rl;

    if (D <= rootVSmall && Q < -rootVSmall)
    {
        const scalar qm = 2*sqrt(-Q);
        scalar ratio = Rl/sqrt(-Q3);
        ratio = max(min(ratio, scalar(1)), scalar(-1));

        const scalar th = ::acos(ratio);

        const scalar r1 = qm*cos(th/3.0) - a2/3.0;
        const scalar r2 =
            qm*cos((th + 2*constant::mathematical::pi)/3.0) - a2/3.0;
        const scalar r3 =
            qm*cos((th + 4*constant::mathematical::pi)/3.0) - a2/3.0;

        scalar rMin = min(r1, min(r2, r3));
        scalar rMax = max(r1, max(r2, r3));
        scalar rMid = r1 + r2 + r3 - rMin - rMax;

        roots.setSize(3);
        roots[0] = rMin;
        roots[1] = rMid;
        roots[2] = rMax;
    }
    else
    {
        const scalar D05 = sqrt(max(D, scalar(0)));
        scalar S = 0;

        if (Rl + D05 < 0)
        {
            S = -pow(mag(Rl + D05), 1.0/3.0);
        }
        else
        {
            S = pow(Rl + D05, 1.0/3.0);
        }

        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0/3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0/3.0);
        }

        roots.setSize(1);
        roots[0] = S + Tl - a2/3.0;
    }
}


template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::residualGibbsForRoot
(
    scalar Z,
    scalar A,
    scalar B,
    scalar p,
    scalar T
) const
{
    const scalar sqrt2 = sqrt(2.0);
    const scalar Bsafe = max(B, rootVSmall);
    const scalar ZminusB = Z - Bsafe;
    const scalar arg1 = Z + (1.0 + sqrt2)*Bsafe;
    const scalar arg2 = Z + (1.0 - sqrt2)*Bsafe;

    if (ZminusB <= SMALL || arg1 <= SMALL || arg2 <= SMALL)
    {
        return great;
    }

    if (!mixingB_ || !mixingCoef1_ || !mixingCoef2_ || !mixingCoef3_)
    {
        return great;
    }

    const label nSpec = mixtureX_.size();

    if (!nSpec)
    {
        return great;
    }

    if (aMixWork_.size() != nSpec)
    {
        aMixWork_.setSize(nSpec, 0.0);
    }

    scalar bMix = 0;
    for (label i=0; i<nSpec; ++i)
    {
        bMix += mixtureX_[i]*(*mixingB_)[i];
    }
    bMix = max(bMix, rootVSmall);

    const scalar RT = RR*T;
    const scalar aMix = max(A*sqr(RT)/max(p, rootVSmall), rootVSmall);
    const scalar sqrtT = sqrt(T);

    for (label i=0; i<nSpec; ++i)
    {
        scalar sumAij = 0.0;
        for (label j=0; j<nSpec; ++j)
        {
            const scalar aij =
                (*mixingCoef1_)[i][j]
              - (*mixingCoef2_)[i][j]*sqrtT
              + (*mixingCoef3_)[i][j]*T;

            sumAij += mixtureX_[j]*aij;
        }
        aMixWork_[i] = sumAij;
    }

    const scalar logTerm = log(arg1/arg2);
    const scalar logZminusB = log(ZminusB);

    scalar g = 0.0;

    for (label i=0; i<nSpec; ++i)
    {
        const scalar xi = max(mixtureX_[i], rootVSmall);
        const scalar bi = (*mixingB_)[i];
        const scalar sumAij = aMixWork_[i];

        const scalar term1 = bi/bMix*(Z - 1.0) - logZminusB;
        const scalar term2 =
            A/(Bsafe*sqrt2)*(2.0*sumAij/aMix - bi/bMix)*logTerm;

        const scalar lnPhi = term1 - term2;

        g += xi*lnPhi;
    }

    return g;
}


template<class Specie>
inline void Foam::PengRobinson<Specie>::setMixingTables
(
    const List<scalar>& bList,
    const List<List<scalar>>& coef1,
    const List<List<scalar>>& coef2,
    const List<List<scalar>>& coef3
) const
{
    mixingB_ = &bList;
    mixingCoef1_ = &coef1;
    mixingCoef2_ = &coef2;
    mixingCoef3_ = &coef3;
}


template<class Specie>
inline void Foam::PengRobinson<Specie>::setMixtureState
(
    const UList<scalar>& X,
    scalar* lastZPtr,
    bool* lastZValidPtr
) const
{
    mixtureX_.setSize(X.size());
    scalar sumX = 0.0;
    forAll(X, i)
    {
        mixtureX_[i] = max(X[i], rootVSmall);
        sumX += mixtureX_[i];
    }

    if (sumX > rootVSmall)
    {
        mixtureX_ /= sumX;
    }

    if (aMixWork_.size() != mixtureX_.size())
    {
        aMixWork_.setSize(mixtureX_.size(), 0.0);
    }

    lastZPtr_ = lastZPtr;
    lastZValidPtr_ = lastZValidPtr;

    hasMixtureState_ = true;
    hasCachedZ_ = false;
}


template<class Specie>
inline void Foam::PengRobinson<Specie>::clearMixtureState() const
{
    mixtureX_.clear();
    aMixWork_.clear();
    lastZPtr_ = nullptr;
    lastZValidPtr_ = nullptr;
    hasMixtureState_ = false;
    hasCachedZ_ = false;
}


template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    const scalar aAlpha = coef1_ - coef2_*sqrt(T) + coef3_*T;
    const scalar daAlpha = -coef2_/(2*sqrt(T))+coef3_;

    const scalar A = aAlpha*p/sqr(RR*T);
    scalar B = b_*p/(RR*T); 
    if(B <=0 ) { B = 1e-16; } 

    const scalar Z = this->Z(p, T);

    //const scalar M = (sqr(Z) + B*Z)/(Z - B);
    const scalar M = (sqr(Z) + 2*B*Z-sqr(B))/(Z - B); //bug finding by Nicola
    const scalar N = daAlpha*B/(b_*RR);

    return this->R()*sqr(M - N)/(sqr(M) - 2*A*(Z + B));
}


// Functions returns member data for checking 
template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::bM() const
{
    return b_; 
}

template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::coef1() const
{
    return coef1_; 
}

template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::coef2() const
{
    return coef2_; 
}

template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::coef3() const
{
    return coef3_; 
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::PengRobinson<Specie>::operator+=
(
    const PengRobinson<Specie>& srk
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(srk);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        b_ = b_; 
    }
}


template<class Specie>
inline void Foam::PengRobinson<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::PengRobinson<Specie> Foam::operator+
(
    const PengRobinson<Specie>& srk1,
    const PengRobinson<Specie>& srk2
)
{
    Specie sp
    (
        static_cast<const Specie&>(srk1)
      + static_cast<const Specie&>(srk2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return PengRobinson<Specie>
        (
            sp,
            srk1.b_,
            srk1.coef1_,
            srk1.coef2_,
            srk1.coef3_
        );
    }
    else
    {
        const scalar b = srk1.b_;
        return PengRobinson<Specie>
        (
            sp,
            b,
            srk1.coef1_,
            srk1.coef2_,
            srk1.coef3_
        );
    }
}


template<class Specie>
Foam::PengRobinson<Specie> Foam::operator*
(
    const scalar s,
    const PengRobinson<Specie>& srk
)
{
    return PengRobinson<Specie>
    (
        s*static_cast<const Specie&>(srk),
        srk.b_,
        srk.coef1_,
        srk.coef2_,
        srk.coef3_
    );
}


template<class Specie>
Foam::PengRobinson<Specie> Foam::operator==
(
    const PengRobinson<Specie>& srk1,
    const PengRobinson<Specie>& srk2
)
{
    Specie sp
    (
        static_cast<const Specie&>(srk1)
     == static_cast<const Specie&>(srk2)
    );

    const scalar b = srk1.b_;
    return PengRobinson<Specie>
    (
        sp,
        b,
        srk1.coef1_,
        srk1.coef2_,
        srk1.coef3_
    );
}


// ************************************************************************* //
