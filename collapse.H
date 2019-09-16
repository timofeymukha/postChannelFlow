    //Get the coordinates

    const scalarField& y = channelIndexing.y();
    scalarField allY(y.size()+2);

    allY[0] = yBot;
    forAll(y, yI)
    {
        allY[yI+1] = y[yI];
    }
    allY[allY.size()-1] = yTop;

    IOobjectList fieldList
    (
        mesh,
        runTime.timeName()
    );

    forAllConstIters(fieldList, fieldIter)
    {
         const word fieldName = fieldIter()->name();
         Info << fieldName << nl;


         fileName path(fieldIter()->rootPath()/fieldIter()->caseName()/"postProcessing"/"collapsedFields"/fieldIter()->instance());

        if (fieldIter()->typeHeaderOk<volScalarField>(true, true, false))
        {
            mkDir(path);
            volScalarField field
            (
                *fieldIter(),
                mesh
            );

            scalarField values = channelIndexing.collapse(field);

            scalar valBot = patchAverage<volScalarField, scalar>(mesh, field,
                                                         areaBot, patchBot);
            scalar valTop = patchAverage<volScalarField, scalar>(mesh, field,
                                                         areaTop, patchTop);
            scalarField allValues(values.size()+2);

            allValues[0] = valBot;
            forAll(values, valI)
            {
                allValues[valI+1] = values[valI];
            }
            allValues[allValues.size()-1] = valTop;

            makeGraph(allY, allValues, fieldName, path, gFormat);
        }
        else if (fieldIter()->typeHeaderOk<volVectorField>(true, true, false))
        {
            mkDir(path);
            volVectorField field
            (
                *fieldIter(),
                mesh
            );

            scalarField valuesX =
                channelIndexing.collapse(field.component(vector::X)());
            scalarField valuesY =
                channelIndexing.collapse(field.component(vector::Y)());
            scalarField valuesZ =
                channelIndexing.collapse(field.component(vector::Z)());

            vector valBot = patchAverage<volVectorField, vector>
            (
                mesh,
                field,
                areaBot,
                patchBot
            );
            vector valTop = patchAverage<volVectorField, vector>
            (
                mesh,
                field,
                areaTop,
                patchTop
            );
            scalarField allValuesX(valuesX.size()+2);
            scalarField allValuesY(valuesY.size()+2);
            scalarField allValuesZ(valuesZ.size()+2);

            allValuesX[0] = valBot.x();
            allValuesY[0] = valBot.y();
            allValuesZ[0] = valBot.z();

            forAll(valuesX, valI)
            {
                allValuesX[valI+1] = valuesX[valI];
                allValuesY[valI+1] = valuesY[valI];
                allValuesZ[valI+1] = valuesZ[valI];
            }
            allValuesX[allValuesX.size()-1] = valTop.x();
            allValuesY[allValuesY.size()-1] = valTop.y();
            allValuesZ[allValuesZ.size()-1] = valTop.z();

            makeGraph(allY, allValuesX, fieldName+"_X", path, gFormat);
            makeGraph(allY, allValuesY, fieldName+"_Y", path, gFormat);
            makeGraph(allY, allValuesZ, fieldName+"_Z", path, gFormat);
        }
        else if (fieldIter()->typeHeaderOk<volSymmTensorField>(true, true, false))
        {
            mkDir(path);
            volSymmTensorField field
            (
                *fieldIter(),
                mesh
            );

            scalarField valuesXX =
                channelIndexing.collapse(field.component(symmTensor::XX)());
            scalarField valuesXY =
                channelIndexing.collapse(field.component(symmTensor::XY)());
            scalarField valuesXZ =
                channelIndexing.collapse(field.component(symmTensor::XZ)());
            scalarField valuesYY =
                channelIndexing.collapse(field.component(symmTensor::YY)());
            scalarField valuesYZ =
                channelIndexing.collapse(field.component(symmTensor::YZ)());
            scalarField valuesZZ =
                channelIndexing.collapse(field.component(symmTensor::ZZ)());

            symmTensor valBot = patchAverage<volSymmTensorField, symmTensor>
            (
                mesh,
                field,
                areaBot,
                patchBot
            );
            symmTensor valTop = patchAverage<volSymmTensorField, symmTensor>
            (
                mesh,
                field,
                areaTop,
                patchTop
            );

            scalarField allValuesXX(valuesXX.size()+2);
            scalarField allValuesXY(valuesXY.size()+2);
            scalarField allValuesXZ(valuesXZ.size()+2);
            scalarField allValuesYY(valuesXX.size()+2);
            scalarField allValuesYZ(valuesXY.size()+2);
            scalarField allValuesZZ(valuesXZ.size()+2);

            allValuesXX[0] = valBot.component(symmTensor::XX);
            allValuesXY[0] = valBot.component(symmTensor::XY);
            allValuesXZ[0] = valBot.component(symmTensor::XZ);
            allValuesYY[0] = valBot.component(symmTensor::YY);
            allValuesYZ[0] = valBot.component(symmTensor::YZ);
            allValuesZZ[0] = valBot.component(symmTensor::ZZ);

            forAll(valuesXX, valI)
            {
                allValuesXX[valI+1] = valuesXX[valI];
                allValuesXY[valI+1] = valuesXY[valI];
                allValuesXZ[valI+1] = valuesXZ[valI];
                allValuesYY[valI+1] = valuesYY[valI];
                allValuesYZ[valI+1] = valuesYZ[valI];
                allValuesZZ[valI+1] = valuesZZ[valI];
            }

            allValuesXX[allValuesXX.size()-1] = valTop.component(symmTensor::XX);
            allValuesXY[allValuesXY.size()-1] = valTop.component(symmTensor::XY);
            allValuesXZ[allValuesXZ.size()-1] = valTop.component(symmTensor::XZ);
            allValuesYY[allValuesYY.size()-1] = valTop.component(symmTensor::YY);
            allValuesYZ[allValuesYZ.size()-1] = valTop.component(symmTensor::YZ);
            allValuesZZ[allValuesZZ.size()-1] = valTop.component(symmTensor::ZZ);

            makeGraph(allY, allValuesXX, fieldName+"_XX", path, gFormat);
            makeGraph(allY, allValuesYY, fieldName+"_YY", path, gFormat);
            makeGraph(allY, allValuesZZ, fieldName+"_ZZ", path, gFormat);
            makeGraph(allY, allValuesXY, fieldName+"_XY", path, gFormat);
            makeGraph(allY, allValuesXZ, fieldName+"_XZ", path, gFormat);
            makeGraph(allY, allValuesYZ, fieldName+"_YZ", path, gFormat);
        }
        else if (fieldIter()->typeHeaderOk<volTensorField>(true, true, false))
        {
            mkDir(path);

            volTensorField field
            (
                *fieldIter(),
                mesh
            );

            scalarField valuesXX =
                channelIndexing.collapse(field.component(tensor::XX)());
            scalarField valuesXY =
                channelIndexing.collapse(field.component(tensor::XY)());
            scalarField valuesXZ =
                channelIndexing.collapse(field.component(tensor::XZ)());
            scalarField valuesYX =
                channelIndexing.collapse(field.component(tensor::YX)());
            scalarField valuesYY =
                channelIndexing.collapse(field.component(tensor::YY)());
            scalarField valuesYZ =
                channelIndexing.collapse(field.component(tensor::YZ)());
            scalarField valuesZX =
                channelIndexing.collapse(field.component(tensor::ZX)());
            scalarField valuesZY =
                channelIndexing.collapse(field.component(tensor::ZY)());
            scalarField valuesZZ =
                channelIndexing.collapse(field.component(tensor::ZZ)());

            tensor valBot = patchAverage<volTensorField, tensor>
            (
                mesh,
                field,
                areaBot,
                patchBot
            );
            tensor valTop = patchAverage<volTensorField, tensor>
            (
                mesh,
                field,
                areaTop,
                patchTop
            );

            scalarField allValuesXX(valuesXX.size()+2);
            scalarField allValuesXY(valuesXY.size()+2);
            scalarField allValuesXZ(valuesXZ.size()+2);
            scalarField allValuesYX(valuesYX.size()+2);
            scalarField allValuesYY(valuesYY.size()+2);
            scalarField allValuesYZ(valuesYZ.size()+2);
            scalarField allValuesZX(valuesZX.size()+2);
            scalarField allValuesZY(valuesZY.size()+2);
            scalarField allValuesZZ(valuesZZ.size()+2);

            allValuesXX[0] = valBot.component(tensor::XX);
            allValuesXY[0] = valBot.component(tensor::XY);
            allValuesXZ[0] = valBot.component(tensor::XZ);
            allValuesYX[0] = valBot.component(tensor::YX);
            allValuesYY[0] = valBot.component(tensor::YY);
            allValuesYZ[0] = valBot.component(tensor::YZ);
            allValuesZX[0] = valBot.component(tensor::ZX);
            allValuesZY[0] = valBot.component(tensor::ZY);
            allValuesZZ[0] = valBot.component(tensor::ZZ);

            forAll(valuesXX, valI)
            {
                allValuesXX[valI+1] = valuesXX[valI];
                allValuesXY[valI+1] = valuesXY[valI];
                allValuesXZ[valI+1] = valuesXZ[valI];
                allValuesYX[valI+1] = valuesYX[valI];
                allValuesYY[valI+1] = valuesYY[valI];
                allValuesYZ[valI+1] = valuesYZ[valI];
                allValuesZX[valI+1] = valuesZX[valI];
                allValuesZY[valI+1] = valuesZY[valI];
                allValuesZZ[valI+1] = valuesZZ[valI];
            }

            allValuesXX[allValuesXX.size()-1] = valTop.component(tensor::XX);
            allValuesXY[allValuesXY.size()-1] = valTop.component(tensor::XY);
            allValuesXZ[allValuesXZ.size()-1] = valTop.component(tensor::XZ);
            allValuesYY[allValuesYX.size()-1] = valTop.component(tensor::YX);
            allValuesYY[allValuesYY.size()-1] = valTop.component(tensor::YY);
            allValuesYZ[allValuesYZ.size()-1] = valTop.component(tensor::YZ);
            allValuesZZ[allValuesZX.size()-1] = valTop.component(tensor::ZX);
            allValuesZZ[allValuesZY.size()-1] = valTop.component(tensor::ZY);
            allValuesZZ[allValuesZZ.size()-1] = valTop.component(tensor::ZZ);

            makeGraph(allY, allValuesXX, fieldName+"_XX", path, gFormat);
            makeGraph(allY, allValuesXY, fieldName+"_XY", path, gFormat);
            makeGraph(allY, allValuesXZ, fieldName+"_XZ", path, gFormat);
            makeGraph(allY, allValuesYX, fieldName+"_YX", path, gFormat);
            makeGraph(allY, allValuesYY, fieldName+"_YY", path, gFormat);
            makeGraph(allY, allValuesYZ, fieldName+"_YZ", path, gFormat);
            makeGraph(allY, allValuesZX, fieldName+"_ZX", path, gFormat);
            makeGraph(allY, allValuesZY, fieldName+"_ZY", path, gFormat);
            makeGraph(allY, allValuesZZ, fieldName+"_ZZ", path, gFormat);
        }
        else
        {
            Info<< "    Field " << fieldIter()->name() 
                << " not found or its type is unsupported "
                << nl;
        }
    }
