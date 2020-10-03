/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |                    2020 Timofey Mukha
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
    Foam::channelIndex

Description
    Does averaging of fields over layers of cells. Assumes layered mesh.

SourceFiles
    channelIndex.C

\*---------------------------------------------------------------------------*/

#ifndef channelIndex_H
#define channelIndex_H

#include "regionSplit.H"
#include "direction.H"
#include "scalarField.H"
#include "polyMesh.H"
#include "GeometricField.H"
#include "volMesh.H"
#include "fvPatchField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class channelIndex Declaration
\*---------------------------------------------------------------------------*/

class channelIndex
{

    // Private data

        static const Enum<vector::components> vectorComponentsNames_;

        //- Is mesh symmetric
        //const bool symmetric_;

        //- Direction to sort
        const direction dir_;

        //- Per cell the global region
        autoPtr<regionSplit> cellRegion_;

        //- Per global region the number of cells (scalarField so we can use
        //  field algebra)
        scalarField regionCount_;

        //- From sorted region back to unsorted global region
        labelList sortMap_;

        //- Sorted component of cell centres
        scalarField y_;

        //- Start patch index
        labelList bottomPatchIndices_;

        //- Opposite patch index
        labelList topPatchIndices_;



    // Private Member Functions

        void walkOppositeFaces
        (
            const polyMesh& mesh,
            const labelList& startFaces,
            boolList& blockedFace
        );

        void calcLayeredRegions
        (
            const polyMesh& mesh,
            const boolList& blockedFace
        );

        //- No copy construct
        channelIndex(const channelIndex&) = delete;

        //- No copy assignment
        void operator=(const channelIndex&) = delete;


public:

    // Constructors

        //- Construct from dictionary
        channelIndex(const polyMesh&, const dictionary&);

        //- Construct from supplied starting faces
        channelIndex
        (
            const polyMesh& mesh,
            const labelList& startFaces,
            //const bool symmetric,
            const direction dir
        );


    // Member Functions

        //- Find the indices of the seed patches
        void findBottomPatchIndices
        (
            const polyMesh& mesh,
            const wordList& patchNames
        );

        void findBottomPatchIndices
        (
            const polyMesh& mesh,
            const labelList& startFaces
        );

        //- Find the indices of the patches opposite the seed patches
        void findTopPatchIndices
        (
            const polyMesh& mesh,
            const boolList& blockedFace
        );

        //- Check that the bottom and top patches have the same number
        //  of faces in total
        void checkPatchSizes
        (
            const polyBoundaryMesh& bMesh
        );

        //- Sum field per region
        template<class T>
        Field<T> regionSum(const Field<T>& cellField) const;

        //- Collapse a field to a line
        template<class T>
        Field<T> collapse
        (
            const Field<T>& vsf
        //    const bool asymmetric=false
        ) const;

        //- Collapse a field to a value on the bottom and top patches
        template<class T>
        Pair<T> collapseBoundary
        (
            const polyBoundaryMesh&,
            const typename GeometricField<T, fvPatchField, volMesh>::Boundary&
        );

        // Access

            //- Return the field of Y locations from the cell centres
            const scalarField& y() const
            {
                return y_;
            }

            labelList bottomPatchIndices() const
            {
                return bottomPatchIndices_;
            }

            labelList topPatchIndices() const
            {
                return topPatchIndices_;
            }

            direction dir() const
            {
                return dir_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "channelIndexTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
