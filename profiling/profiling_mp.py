import os, sys
import time

import multiprocessing as mp

import cProfile, pstats, io

# Number of processes to spawn
procs = 4

# function to be multiprocessed
def worker_process(N):
    # turn profiling on
    pr = cProfile.Profile()
    pr.enable()

    # thing being measured
    time.sleep(5)

    # turn off profiling and dump info into log files in current dir
    pr.disable()
    s = io.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    with open("worker_{}.log".format(N), 'w') as f:
        print(s.getvalue(), file=f)
            

def main():

    mp.set_start_method('spawn')
    processes = []


    single = mp.Process(target=worker_process, args=("single",), name='single_worker')
    single.start()
    single.join()

    for i in range(procs):
        worker = mp.Process(target=worker_process, args=(i,), daemon=True, name='worker_{}'.format(i))
        worker.start()
        processes.append(worker)
        
    for p in processes:
        p.join()


if __name__ == '__main__':
    main()