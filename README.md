# HarvesterCPP
Harvester C++

Implementation of Harvester in C++. For undirected graphs for now.

Version of Harvester in Python can be found [here](https://github.com/nd7141/influence-maximization/blob/e9dd1b354c7ce90cbe6fbcb4a866d43bb178a9ad/IC/ArbitraryP/Harvester.py). 

Compile:
--------
g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"Harvester.d" -MT"Harvester.d" -o "Harvester.o" "../Harvester.cpp"

Execute:
--------
./Harvester [filename of the weighted graph] [number of vertices] [folder with PWs] [number of seeds] [filename of the seeds]

#####Example:
./Harvester dblp317080_mv1.txt 317080 PW/Random/dblp317080_mv1/ 50 seeds_dblp317080_mv1.txt

#####Format of _the weighted graph_:
Each line has:
node1 node2 probability
#####Example:
######0 1 0.04
######0 2 0.02

#####Format of _PW (Possible World)_:
Each line has:
node1 node2
#####Example:
######1 17
######2 112

The file for the seeds will be overwritten is exists. 
