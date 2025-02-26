""""
Run A3FE calculations for all ligands. Run interactively in ipython
using %run run_all_a3fe.py.
"""

import os
import a3fe as a3
import logging
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)
fh = logging.FileHandler("run_a3fe.log")
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
fh.setFormatter(formatter)
logger.addHandler(fh)


def get_cfg(fast: bool = True) -> a3.SystemPreparationConfig:
    cfg = a3.SystemPreparationConfig()
    cfg.forcefields["ligand"] = "openff_unconstrained-2.2.1"
    cfg.forcefields["protein"] = "ff14SB"
    lambda_values = {
        a3.LegType.BOUND: {
            a3.StageType.RESTRAIN: [0.0, 1.0],
            a3.StageType.DISCHARGE: [0.0, 0.291, 0.54, 0.776, 1.0],
            a3.StageType.VANISH: [
                0.0,
                0.026,
                0.054,
                0.083,
                0.111,
                0.14,
                0.173,
                0.208,
                0.247,
                0.286,
                0.329,
                0.373,
                0.417,
                0.467,
                0.514,
                0.564,
                0.623,
                0.696,
                0.833,
                1.0,
            ],
        },
        a3.LegType.FREE: {
            a3.StageType.DISCHARGE: [0.0, 0.222, 0.447, 0.713, 1.0],
            a3.StageType.VANISH: [
                0.0,
                0.026,
                0.055,
                0.09,
                0.126,
                0.164,
                0.202,
                0.239,
                0.276,
                0.314,
                0.354,
                0.396,
                0.437,
                0.478,
                0.518,
                0.559,
                0.606,
                0.668,
                0.762,
                1.0,
            ],
        },
    }
    cfg.lambda_values = lambda_values
    cfg.slurm = True
    if fast:  # Drop the equilibration times
        cfg.runtime_npt_unrestrained = 50 # ps
        cfg.runtime_npt = 50 # ps
        cfg.ensemble_equilibration_time = 100 # ps

    return cfg

# Sort calculation paths first by the ligand number, then by the conformer number
calc_paths = [d for d in os.listdir() if "lig" in d]
calc_paths.sort(key=lambda x: (int(x.split('_')[1]), int(x.split('_')[2])))

# Set the calcuations up with a "fast" config
cfg = get_cfg(fast=True)

# Try and run every calculation, but tolerate failures
failed_calcs = []
for i, calc_path in enumerate(calc_paths):
    calc_name = calc_path.split("/")[-1]
    percent_complete = i / len(calc_paths) * 100
    logger.info(f"Running calculation {calc_name}, {percent_complete:.2f}% complete")
    try:
        logger.info(50 * "=")
        # Skip calculations which are already complete (e.g. overall_stats.dat in the output directory)
        if os.path.exists(os.path.join(calc_path, "output", "overall_stats.dat")):
            logger.info(f"Calculation {calc_name} already complete")
            continue
        logger.info(f"Setting up calculation {calc_name}")
        calc = a3.Calculation(base_dir=calc_path)
        calc.setup()
        logger.info(f"Running calculation {calc_name}")
        calc.run(adaptive=False, runtime=0.1)
        calc.wait()
        # Don't keep the heavy trajectory and restart files
        for leg in calc.legs:
            for stage in leg.stages:
                stage.lighten()
        logger.info(f"Analysing calculation {calc_name}")
        calc.set_equilibration_time(0.02)
        calc.analyse(slurm=True)
        calc.save()
        # Clean up
        calc._close_logging_handlers()
        del calc
        plt.close("all")
        logger.info(f"Finished calculation {calc_name}")

    except Exception as e:
        logger.error(f"Failed to run calculation {calc_name} due to {e}")
        failed_calcs.append((calc_path, e))
