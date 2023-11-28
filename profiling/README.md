


### Single processng using module

/usr/bin/time -v python3 -m cProfile -s cumtime profiling_single.py

### Multiprocessing and complex code

/usr/bin/time -v python3 profiling_mp.py

### For MacOS users, use the following command

brew install gnu-time

gtime --verbose python3 -m cProfile -s cumtime profiling_single.py