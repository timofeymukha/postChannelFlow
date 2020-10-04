/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation, 2020 Timofey Mukha	
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

Application
    postChannel

Group
    grpPostProcessingUtilities

Description
    Post-processes data from channel flow calculations.

    Assuming that the mesh is periodic in the x and z directions, collapse
    fields to a line and print them to postProcesing/collapsedFields.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "channelIndex.H"
#include "makeGraph.H"

#include "OSspecific.H"
#include "IOobjectList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
Foam::wordList comptNames()
{
    return wordList();
}


template<>
Foam::wordList comptNames<vector>()
{
    return {"_X", "_Y", "_Z"};
}


template<>
Foam::wordList comptNames<symmTensor>()
{
    return {"_XX", "_XY", "_XZ", "_YY", "_YZ", "_ZZ"};
}


template<>
Foam::wordList comptNames<tensor>()
{
    return {"_XX", "_XY", "_XZ", "_YX", "_YY", "_YZ", "_ZX", "_ZY", "_ZZ"};
}

template<>
Foam::wordList comptNames<sphericalTensor>()
{
    return {"_II"};
}


template<class T>
void writeToFile
(
    const scalarField& y,
    const Field<T>& values,
    const word name,
    const fileName path,
    const word format,
    const wordList comptNames)
{
    for (label i=0; i<comptNames.size(); ++i)
    {
        makeGraph(y, values.component(i), name+comptNames[i], path, format);
    }
}

template<>
void writeToFile
(
    const scalarField& y,
    const Field<scalar>& values,
    const word name,
    const fileName path,
    const word format,
    const wordList comptNames)
{
    makeGraph(y, values, name, path, format);
}

using HashType = Foam::HashTable
<
    Foam::IOobject*,
    Foam::word,
    Foam::string::hash
>::const_iterator;

template<class T>
void collapse
(
    HashType & fieldIter,
    const fvMesh & mesh,
    const channelIndex & channelIndexing,
    const word format
)
{
    using FieldType=GeometricField<T, fvPatchField, volMesh>;

    if (fieldIter()->typeHeaderOk<FieldType>(true, true, false))
    {
        const word fieldName = fieldIter()->name();
        Info<<"    " << fieldName << endl;

        fileName path
        (
            fieldIter()->rootPath()/fieldIter()->caseName()/
            "postProcessing"/"collapsedFields"/fieldIter()->instance()
        );

        mkDir(path);
        FieldType field
        (
            *fieldIter(),
            mesh
        );

        Pair<T> boundaryValues =
            channelIndexing.collapseBoundary<T>
            (
                mesh.boundaryMesh(),
                field.boundaryField()
            );

        Field<T> internalValues =
            channelIndexing.collapse(field);

        Field<T> allValues(internalValues.size() + 2);
        allValues[0] = boundaryValues[0];
        allValues[allValues.size()-1] = boundaryValues[1];

        for (label i=0; i<internalValues.size(); ++i)
        {
            allValues[i+1] = internalValues[i];
        }

        tmp<scalarField> y = channelIndexing.y(mesh.boundaryMesh());

        const wordList cNames = comptNames<T>();
        writeToFile<T>(y, allValues, fieldName, path, format, cNames);
    }
}




int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Post-process data from channel flow calculations"
    );

    argList::noParallel();
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createNamedMesh.H"
#   include "readTransportProperties.H"

    const word& gFormat = runTime.graphFormat();

    // Setup channel indexing for averaging over channel down to a line

    IOdictionary channelDict
    (
        IOobject
        (
            "postChannelDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    channelIndex channelInd(mesh, channelDict);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Collapsing fields for time " << runTime.timeName() << endl;

        IOobjectList fieldList
        (
            mesh,
            runTime.timeName()
        );

        forAllConstIters(fieldList, fieldIter)
        {
            const word fieldName = fieldIter()->name();

            collapse<scalar>(fieldIter, mesh, channelInd, gFormat);
            collapse<vector>(fieldIter, mesh, channelInd, gFormat);
            collapse<sphericalTensor>(fieldIter, mesh, channelInd, gFormat);
            collapse<symmTensor>(fieldIter, mesh, channelInd, gFormat);
            collapse<tensor>(fieldIter, mesh, channelInd, gFormat);

        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
