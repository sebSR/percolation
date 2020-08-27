import simulation.percolation as pr
import numpy as np



def main_burning(lattice_size: int, trials: int, p_start: float,
                 p_end: float, dp: float):

    file_name = f"burning_L_{lattice_size}_T_{trials}.dat"
    probability, p_end, dp = p_start*100, p_end*100, dp*100

    with open(file_name, "w") as file_to_save:
        while probability <= p_end:
            new_forest = pr.Forest(lattice_size, probability, trials)
            new_forest.simulation_burning()
            pattern = f"{round(probability/100,3)} \t {round(new_forest.percolation_p,2)} \t {new_forest.min_average_path} \n"
            file_to_save.write(pattern)
            print(probability)
            probability += dp
