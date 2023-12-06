import cProfile
import pstats
import io
from numpy.random import RandomState

from cellregmap import compute_maf

def profile_compute_maf():
    # Sample data for profiling
    random = RandomState(10)
    n = 30                               # number of samples (cells)
    G = 1.0 * (random.rand(n, 1) < 0.2)  # SNP vector

    pr = cProfile.Profile()
    pr.enable()

    compute_maf(G)

    pr.disable()
    s = io.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue())

def main():
    profile_compute_maf()

if __name__ == "__main__":
    main()
