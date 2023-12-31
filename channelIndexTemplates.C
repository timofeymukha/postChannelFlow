/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation, 2020 Timofey Mukha
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

#include "channelIndex.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::Field<T> Foam::channelIndex::regionSum(const Field<T>& cellField) const
{
    Field<T> regionField(cellRegion_().nRegions(), Zero);

    forAll(cellRegion_(), celli)
    {
        regionField[cellRegion_()[celli]] += cellField[celli];
    }

    // Global sum
    Pstream::listCombineGather(regionField, plusEqOp<T>());
    Pstream::listCombineScatter(regionField);

    return regionField;
}


template<class T>
Foam::Field<T> Foam::channelIndex::collapse
(
    const Field<T>& cellField
) const
{
    // Average and order
    const Field<T> summedField(regionSum(cellField));

    Field<T> regionField
    (
        summedField
      / regionCount_,
        sortMap_
    );

    return regionField;
}

template<class T>
Foam::Pair<T> Foam::channelIndex::collapseBoundary
(
    const polyBoundaryMesh& bMesh,
    const typename GeometricField<T, fvPatchField, volMesh>::Boundary & boundaryField    
) const
{

    // A pair of values, corresponding to bottom and top patches
    Pair<T> result(pTraits<T>::zero, pTraits<T>::zero);
    
    Pair<labelList> patchIndices(bottomPatchIndices_, topPatchIndices_);

    for (label i=0; i<2; ++i)
    {
        scalar totalArea = 0;
        T areaWeightedValues = pTraits<T>::zero;

        for (label pI : patchIndices[i])
        {
            scalarField magSf = mag(bMesh[pI].faceAreas());
            areaWeightedValues += gSum(magSf*boundaryField[pI]);
            totalArea += gSum(magSf);
        }
        
        result[i] = areaWeightedValues/totalArea;
    }

    return result;
}

// ************************************************************************* //
