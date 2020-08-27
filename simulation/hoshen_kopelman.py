import simulation.percolation as pr
import numpy as np
import time


def main_hoshen_kopelman(lattice_size: int, trials: int, p_start: float,
                         p_end: float, dp: float):

    file_name = f"clustering_L_{lattice_size}_T_{trials}.dat"
    probability, p_end, dp = p_start*100, p_end*100, dp*100

    with open(file_name, "w") as file_to_save:
        while probability <= p_end:
            new_forest = pr.Forest(lattice_size, probability, trials)
            new_forest.simulation_hoshen_kopelman()
            pattern = f"{round(probability/100,3)} \t {int(new_forest.max_average_cluster)} \t {int(new_forest.average_cluster)} \n"
            file_to_save.write(pattern)
            print(probability)
            probability += dp
