# CellRegMap Speed-up Testing

This repository is for testing the speed-up of the CellRegMap package. Note: This is not the repository for the CellRegMap package itself.

## Purpose
* Use cell profiling to identify bottlenecks in the original CellRegMap package
* Improve bottlenecks
  * within python (better structures, better python)
  * different languages (e.g. C) 
  * different external libraries & functions
  * dependency 

## CellRegMap Installation
* For both original & speed-up version: follow the installation instructions from [Here](https://github.com/YPZ404/CellRegMap_Optimized)

## Python Version
Using python3 between 3.7 and 3.9 (3.9 is preferred)

## Usage Instructions
* Clone the repository
```
git clone https://github.com/annacuomo/CellRegMap_speedup_testing.git
```
* The [Profiling](./profiling) folder contains the simple cProfile scripts
* The [test_data](./test_data) folder contains the simulated data for testing and validation of any changes made to the original CellRegMap package
* The [basic1](./basic1) folder is used for testing basic usage of CellRegMap with random generated data
* The [more_data](./more_data) folder is used for testing larger scale simulated data
