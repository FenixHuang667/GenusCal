# Gcal

Gcal is a program to compute the topological genus of an input RNA structure. 

The program is written in C. Please cite the software as specified at the bottom of the paper. 


### Prerequisites

GNU Automake
C Standard Library


### Installing

Use the command: 

```
./configure
make
```
The executive file is Gcal in ./src. 


## Running 

### Running HamSampler 

There are several options for HamSampler. Use the following command to see all options. 
```
./src/Gcal -h  
```

To specify an input file, use the command
```
./src/Gcal -i input.in 
```
The default input file is input.in. 

To specify an output file, use the command
```
./src/Gcal -o output.out 
```
The default output file is output.out. 


### Input file style

We accept the input structure in dot/bracket form. Bracket "(" pairs to ")", "[" to "]", "{" to "}". Unpaired vertex is denoted by ".". Currently we support at most three types of bracket. A valid input can be found in input.in. 


### Output file
We output the genus after each structure. An example of output file can be found in output.out.    


## Contact

If you have any question or bug report about HamSampler, please feel free to conatct the author Fenix Huang by fenixprotoss@gmail.com.  




