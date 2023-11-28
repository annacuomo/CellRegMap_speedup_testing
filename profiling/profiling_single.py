import os, sys
import time

def worker(N):
    time.sleep(N)


def main():
    N = 5
    worker(N)



if __name__ == '__main__':
    main()