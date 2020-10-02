/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "channelIndex.H"
#include "boolList.H"
#include "syncTools.H"
#include "OFstream.H"
#include "meshTools.H"
#include "Time.H"
#include "SortableList.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::vector::components
>
Foam::channelIndex::vectorComponentsNames_
({
    { vector::components::X, "x" },
    { vector::components::Y, "y" },
    { vector::components::Z, "z" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Determines face blocking
void Foam::channelIndex::walkOppositeFaces
(
    const polyMesh& mesh,
    const labelList& startFaces,
    boolList& blockedFace
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const label nBnd = mesh.nBoundaryFaces();

    DynamicList<label> frontFaces(startFaces);
    forAll(frontFaces, i)
    {
        label facei = frontFaces[i];
        blockedFace[facei] = true;
    }

    while (returnReduce(frontFaces.size(), sumOp<label>()) > 0)
    {
        // Transfer across.
        boolList isFrontBndFace(nBnd, false);
        forAll(frontFaces, i)
        {
            label facei = frontFaces[i];

            if (!mesh.isInternalFace(facei))
            {
                isFrontBndFace[facei-mesh.nInternalFaces()] = true;
            }
        }
        syncTools::swapBoundaryFaceList(mesh, isFrontBndFace);

        // Add
        forAll(isFrontBndFace, i)
        {
            label facei = mesh.nInternalFaces()+i;
            if (isFrontBndFace[i] && !blockedFace[facei])
            {
                blockedFace[facei] = true;
                frontFaces.append(facei);
            }
        }

        // Transfer across cells
        DynamicList<label> newFrontFaces(frontFaces.size());

        forAll(frontFaces, i)
        {
            label facei = frontFaces[i];

            {
                const cell& ownCell = cells[mesh.faceOwner()[facei]];

                label oppositeFacei = ownCell.opposingFaceLabel(facei, faces);

                if (oppositeFacei == -1)
                {
                    FatalErrorInFunction
                        << "Face:" << facei << " owner cell:" << ownCell
                        << " is not a hex?" << abort(FatalError);
                }
                else
                {
                    if (!blockedFace[oppositeFacei])
                    {
                        blockedFace[oppositeFacei] = true;
                        newFrontFaces.append(oppositeFacei);
                    }
                }
            }

            if (mesh.isInternalFace(facei))
            {
                const cell& neiCell = mesh.cells()[mesh.faceNeighbour()[facei]];

                label oppositeFacei = neiCell.opposingFaceLabel(facei, faces);

                if (oppositeFacei == -1)
                {
                    FatalErrorInFunction
                        << "Face:" << facei << " neighbour cell:" << neiCell
                        << " is not a hex?" << abort(FatalError);
                }
                else
                {
                    if (!blockedFace[oppositeFacei])
                    {
                        blockedFace[oppositeFacei] = true;
                        newFrontFaces.append(oppositeFacei);
                    }
                }
            }
        }

        frontFaces.transfer(newFrontFaces);
    }
}


// Calculate regions.
void Foam::channelIndex::calcLayeredRegions
(
    const polyMesh& mesh,
    const boolList& blockedFace
)
{


    if (false)
    {
        OFstream str(mesh.time().path()/"blockedFaces.obj");
        label vertI = 0;
        forAll(blockedFace, facei)
        {
            if (blockedFace[facei])
            {
                const face& f = mesh.faces()[facei];
                forAll(f, fp)
                {
                    meshTools::writeOBJ(str, mesh.points()[f[fp]]);
                }
                str<< 'f';
                forAll(f, fp)
                {
                    str << ' ' << vertI+fp+1;
                }
                str << nl;
                vertI += f.size();
            }
        }
    }


    // Do analysis for connected regions
    cellRegion_.reset(new regionSplit(mesh, blockedFace));

    Info<< "Detected " << cellRegion_().nRegions() << " layers." << nl << endl;

    // Sum number of entries per region
    regionCount_ = regionSum(scalarField(mesh.nCells(), 1.0));

    // Average cell centres to determine ordering.
    pointField regionCc
    (
        regionSum(mesh.cellCentres())
      / regionCount_
    );

    SortableList<scalar> sortComponent(regionCc.component(dir_));

    sortMap_ = sortComponent.indices();

    y_ = sortComponent;

    //if (symmetric_)
    //{
        //y_.setSize(cellRegion_().nRegions()/2);
    //}
}


void Foam::channelIndex::findBottomPatchIndices
(
    const polyMesh& mesh,
    const wordList& patchNames
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    for (word i : patchNames)
    {
        const label patchI = bMesh.findPatchID(i);

        if (patchI == -1)
        {
            FatalErrorInFunction
                << "Illegal patch " << i
                << ". Valid patches are " << bMesh.name()
                << exit(FatalError);
        }

        bottomPatchIndices_.append(patchI);
    }
}

void Foam::channelIndex::findBottomPatchIndices
(
    const polyMesh& mesh,
    const labelList& startFaces
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    for (label i : startFaces)
    {
        const label patchI = bMesh.whichPatch(i);

        if (!bottomPatchIndices_.found(patchI))
        {
            bottomPatchIndices_.append(patchI);
        }
        
    }
}


void Foam::channelIndex::findTopPatchIndices
(
    const polyMesh& mesh,
    const boolList& blockedFace
)
{
    const polyBoundaryMesh & bMesh = mesh.boundaryMesh();

    for (label i=0; i<bMesh.nFaces(); i++)
    {
        label faceI = mesh.nInternalFaces() + i;

        if (blockedFace[faceI])
        {
            label patchI = bMesh.whichPatch(faceI);
            if ((!bottomPatchIndices_.found(patchI)) && 
                (!topPatchIndices_.found(patchI)))
            {
                topPatchIndices_.append(patchI); 
            }

        }
    }

    if (topPatchIndices_.size() == 0)
    {
        FatalErrorInFunction
            << "Could not find the top patch(es)."
            << exit(FatalError);
    }

    Info<< "The top patches are: ";
    wordList patchNames = bMesh.names();
    for (label i : topPatchIndices_)
    {
        Info<< patchNames[i] << " ";
    }
}

void Foam::channelIndex::checkPatchSizes
(
    const polyBoundaryMesh& bMesh
)
{
    label nBottom = 0;
    label nTop = 0;

    for (label i : bottomPatchIndices_)
    {
        nBottom += bMesh[i].size();
    }
    for (label i : topPatchIndices_)
    {
        nTop += bMesh[i].size();
    }
    
    if (nBottom != nTop)
    {
        FatalErrorInFunction
            << "The total number of faces on the bottom and top patches are"
            << "unequal. Bottom faces: "<< nBottom << " Top faces: "
            << nTop << "."
            << exit(FatalError);
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::channelIndex::channelIndex
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    //symmetric_(dict.get<bool>("symmetric")),
    dir_(vectorComponentsNames_.get("component", dict))
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    const wordList patchNames(dict.get<wordList>("patches"));

    // Get the seed patch indices from the patch names
    findBottomPatchIndices(mesh, patchNames);

    // Sum the number of faces on the seed patches
    label nFaces = 0;

    forAll(patchNames, i)
    {
        nFaces += bMesh[bottomPatchIndices_[i]].size();
    }

    labelList startFaces(nFaces);
    nFaces = 0;

    forAll(patchNames, i)
    {
        const polyPatch& pp = bMesh[patchNames[i]];

        forAll(pp, j)
        {
            startFaces[nFaces++] = pp.start()+j;
        }
    }

    boolList blockedFace(mesh.nFaces(), false);
    walkOppositeFaces
    (
        mesh,
        startFaces,
        blockedFace
    );


    // Find 
    findTopPatchIndices(mesh, blockedFace);

    checkPatchSizes(bMesh);

    // Calculate regions.
    calcLayeredRegions(mesh, blockedFace);
}


Foam::channelIndex::channelIndex
(
    const polyMesh& mesh,
    const labelList& startFaces,
    //const bool symmetric,
    const direction dir
)
:
    //symmetric_(symmetric),
    dir_(dir)
{
    boolList blockedFace(mesh.nFaces(), false);
    walkOppositeFaces
    (
        mesh,
        startFaces,
        blockedFace
    );
    

    findBottomPatchIndices(mesh, startFaces);
    findTopPatchIndices(mesh, blockedFace);

    checkPatchSizes(mesh.boundaryMesh());
    
    // Calculate regions.
    calcLayeredRegions(mesh, blockedFace);
}


// ************************************************************************* //
