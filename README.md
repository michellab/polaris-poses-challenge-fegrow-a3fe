# polaris-poses-challenge-fegrow-a3fe

Code and inputs to reproduce the ABFE component of the "newcastle-edinburgh-fegrow-a3fe" submission to the [Polaris Pose Prediction Challenge](https://polarishub.io/competitions/asap-discovery/antiviral-ligand-poses-2025), from Finlay Clark, Asma Feriel Khoualdi, Josh Horton, Julien Michel and Daniel Cole. The total cost of the ABFE runs was ~ 7000 GPU hours (for ~ 700 ABFE calculations of 5 replicate runs each), mainly on RTX3080 GPUs.

## Method Details

 - Poses generated with constrained geometry optimisation with ML/MM (ANI), using [FEGrow](https://github.com/cole-group/FEgrow)
 - Poses scored using fast ABFE runs (0.1 ns / lambda window, 5 replicate runs with independent equilibration) using [A3FE](https://github.com/michellab/a3fe)

For more details, see

 - `run_a3fe.py`: The script used to set-up and run all ABFE calculations. Details e.g. the force fields used.
 - The SOMD config files e.g. `mers-040225/input/template_config.cfg` and `mers-040225/input/template_config_charged.cfg` - these are the config files used to run the ABFE calculations for neutral and charged ligands, respectively, using [SOMD](https://siremol.org/tutorials/somd/Binding_free_energy/Production.html).
 - The notebooks used to select the input for the ABFE calculations, e.g. `mers-040225/setup.ipynb`. These generate ABFE input files for distinct poses from the FEGrow input, as judged by RMSD > 0.5 Å from any other poses for the same ligand.

## To Reproduce

Note that this repository only contains the ABFE component of the submission, and begins from the poses generated with FEGrow.

 - Ensure that you have GROMACS and SLURM installed
 - Create and activate the environment from `environment.yml` e.g. `mamba env create -f environment.yml`, `mamba activate polaris-challenge-a3fe`. Alternatively, follow the install instructions for [A3FE](https://github.com/michellab/a3fe), then install polaris.
 - FEGrow output poses to be used for ABFE input are in:
    - `mers-040225`: MERS poses generated by FEGrow
    - `sars-290125-plato`: SARS poses generated by FEGrow
    - `sars-extra-poses-210225`: Additional SARS poses generated by FEGrow (for ligands which were not successfully generated in the initial FEGrow run)
- For each of these directories:
    - Modify `run_somd.sh` in the `input` directories (e.g. `mers-040225/input`) to so that the SLURM settings are correct for your cluster.
    - Create the input directories for the ABFE runs by running `setup.ipynb`. This only generates ABFE input files for distinct poses, as judged by RMSD > 0.5 Å from any other poses for the same ligand.
    - Run all of the ABFE runs by running `python ../run_a3fe.py` (e.g. in a tmux session so that it persists after you log out).
- Extract the highest-scored poses and submit with `submission/submit.ipynb`. Note that for a few ligands, the highest-scored docked pose was directly submitted with no attempt to run ABFE (inputs from `sars-gnina` and `mers-gnina`).
