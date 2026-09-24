"""Run the three ASE-neighbor pilot benchmarks with mandatory LAMMPS checks.

This entry point shares the implementation and options of benchmark_lammps.py.
Historical cluster/native-FIRE timings remain in validation/zno and Git history.
"""

from benchmark_lammps import main


if __name__ == "__main__":
    main()
