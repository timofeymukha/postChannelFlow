# README #

This is an enhancment of the exisitng OpenFOAM utility called `postChannel`, which allows one to average the solution of a channel flow simulation over the stream- and spanwise directions.
To run the utility just write postChannelFlow in the command-line.
The utility expectes a file called postChannelDict in the constant directory.
A sample dictionary can be found in the repository.

Some differences with the `postChannel`
- Averages all the fields you have in the time directory.
- Averages data on the wall patches.
- Does *not* average across the channel centerline, so you get the full profile across the channel.

Works with OpenFOAM v2212, with the ambition to always support the latest version from OpenCFD.

![Cirrus CI - Task and Script Build Status](https://img.shields.io/cirrus/github/timofeymukha/postChannelFlow)

**This offering is not approved or endorsed by OpenCFD Limited, producer
and distributor of the OpenFOAM software and owner of the OPENFOAM®  and
OpenCFD®  trade marks.**
