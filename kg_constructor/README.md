# Knowledge Graph Constructor
This folder contains source code for constructing the inconsistency-free knowledge graph.

## Directories
* <code>[./data](./data)</code>: All input data go here.
* <code>[./output](./output)</code>: All output files go here.
* <code>[./src](./src)</code>: Contains all source code used for knowledge graph construction.

## Clean the output directory.
Current output directory consists of results used in the paper. If you wish to run the code and obtain new results, please remove all files and directories under it.
```
rm -r ./output/*
```

## Identify the current branch.
The current repo has two directories master and refactor. Refactor is the updated one and we will run this. 
Find out your current brach by using the following. 

```
git branch --show-current 
```
You can also use the following.
```
git branch
```
The branch highlighted in green is your current branch.


## How to Switch different branches 
Switch to the desired branch using the following. 
git checkout <your desired branch>
```
git checkout refactor
```


rm -r ./output/*
## How to Run
Using the toy graph makes sure everything is setup correctly.
```
./construct_toy_kg.sh
```

Ecoli knowledge graph can be constructed using the script below.
```
./construct_ecoli_kg.sh
```
