	//Get the corrdinates
    const scalarField& y = channelIndexing.y();

	forAll(fieldNames, i)
	{
		//Create the IOobject for the field	
    	IOobject fieldHeader
    	(
        	fieldNames[i],
        	runTime.timeName(),
        	mesh,
        	IOobject::MUST_READ
    	);

		//Check if the field was found, if not, warn and do nothing
		if (!fieldHeader.headerOk())
		{
			Info<< "    No " << fieldNames[i] << " field" << endl;
			continue;
		}

		//set the path
		fileName path(fieldHeader.rootPath()/fieldHeader.caseName()/"postProcessing"/"collapsedFields"/fieldHeader.instance());
    mkDir(path);

		//Check the type of the field
		word type = fieldHeader.headerClassName();

		if ( type == "volScalarField" )
		{
			volScalarField field
			(
				fieldHeader,
				mesh
			);

			scalarField values = channelIndexing.collapse(field);
			makeGraph(y, values, fieldNames[i], path, gFormat);
		}
		else if ( type == "volVectorField" )
		{
			volVectorField field
			(
				fieldHeader,
				mesh
			);

			scalarField valuesX = channelIndexing.collapse(field.component(vector::X)());
			scalarField valuesY = channelIndexing.collapse(field.component(vector::Y)());
			scalarField valuesZ = channelIndexing.collapse(field.component(vector::Z)());
			makeGraph(y, valuesX, fieldNames[i]+"_X", path, gFormat);
			makeGraph(y, valuesY, fieldNames[i]+"_Y", path, gFormat);
			makeGraph(y, valuesZ, fieldNames[i]+"_Z", path, gFormat);
		}
		else if ( type == "volSymmTensorField" )
		{
			volSymmTensorField field
			(
				fieldHeader,
				mesh
			);

			scalarField valuesXX = channelIndexing.collapse(field.component(symmTensor::XX)());
			scalarField valuesXY = channelIndexing.collapse(field.component(symmTensor::XY)());
			scalarField valuesXZ = channelIndexing.collapse(field.component(symmTensor::XZ)());
			scalarField valuesYY = channelIndexing.collapse(field.component(symmTensor::YY)());
			scalarField valuesYZ = channelIndexing.collapse(field.component(symmTensor::YZ)());
			scalarField valuesZZ = channelIndexing.collapse(field.component(symmTensor::ZZ)());

			makeGraph(y, valuesXX, fieldNames[i]+"_XX", path, gFormat);
			makeGraph(y, valuesYY, fieldNames[i]+"_YY", path, gFormat);
			makeGraph(y, valuesZZ, fieldNames[i]+"_ZZ", path, gFormat);
			makeGraph(y, valuesXY, fieldNames[i]+"_XY", path, gFormat);
			makeGraph(y, valuesXZ, fieldNames[i]+"_XZ", path, gFormat);
			makeGraph(y, valuesYZ, fieldNames[i]+"_YZ", path, gFormat);
		}
		else
		{
			Info << "    Unsupported field type, ignoring field " << fieldNames[i];
		}

	}
