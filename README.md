# README #

This is an enhancment of the exisitng OpenFOAM utility called postChannel, which allows one to average the solution of a channel flow simulation over the stream- and spanwise directions.
The original utility is improved by allowing the user to specify the fields that should be averaged.
To run the utility just write postChannelFlow in the command-line.
The utility expectes a file called postChannelDict in the constant directory.
A sample dictionary can be found in the repository.