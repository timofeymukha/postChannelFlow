# README #

**I am currently in the process of trying to merge this utility into the main code of the .com fork of OF.
Therefore, things are rapidly changing right now, I work directly on the master branch.
I recommend checking out an older commit for a stable version of the utility that works with older versions
of OF**

This is an enhancment of the exisitng OpenFOAM utility called postChannel, which allows one to average the solution of a channel flow simulation over the stream- and spanwise directions.
The original utility is improved by allowing the user to specify the fields that should be averaged.
To run the utility just write postChannelFlow in the command-line.
The utility expectes a file called postChannelDict in the constant directory.
A sample dictionary can be found in the repository.


**This offering is not approved or endorsed by OpenCFD Limited, producer
and distributor of the OpenFOAM software and owner of the OPENFOAM®  and
OpenCFD®  trade marks.**