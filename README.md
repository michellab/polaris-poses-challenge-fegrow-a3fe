# Polaris Pose Prediction Challenge - ABFE Component of "newcastle-edinburgh-fegrow-a3fe-equilibrated-no-training" Submission

Code and inputs to reproduce the absolute binding free energg (ABFE) component of the "newcastle-edinburgh-fegrow-a3fe-equilibrated-no-training" submission to the [Polaris Pose Prediction Challenge](https://polarishub.io/competitions/asap-discovery/antiviral-ligand-poses-2025), from Finlay Clark, Asma Feriel Khoualdi, Josh Horton, Julien Michel and Daniel Cole. This involved growing the compounds into the pocket using [FEGrow](https://github.com/cole-group/FEgrow), followed by scoring of the poses with ABFE using [A3FE](https://github.com/michellab/a3fe). The FEGrow component is described in detail [here](https://github.com/cole-group/polaris-fegrow/tree/main). The total cost of the ABFE runs was ~ 7000 GPU hours (for ~ 700 ABFE calculations of 5 replicate runs each), mainly on RTX3080 GPUs.

## Method Details

 1. Poses were generated with constrained geometry optimisation with ML/MM (ANI), using [FEGrow](https://github.com/cole-group/FEgrow). This is described in detail [here](https://github.com/cole-group/polaris-fegrow/tree/main). **Very little training data were used** - the core structures were taken from very few complexes, which was likely a poor choice. A better approach would likely have been to search the training data by maximum common substructure, and use the most similar ligand as the core to start the FEGrow runs from.
 2. Poses were scored using fast ABFE runs (0.1 ns / lambda window, 5 replicate runs with independent equilibration) using [A3FE](https://github.com/michellab/a3fe). GROMACS 2024 (.2 and .5) (through A3FE) was used for system preparation. Systems were parameterised with the AMBER FF14SB, OpenFF Sage 2.2.1,and TIP3P force fields, solvated with 0.15 M NaCl and energy minimized for 1000 steps, followed by equilibration in the NVT ensemble (5 ps with all nonsolvent atoms restrained and heating from 0 to 298 K, followed by 50 ps with restraints on all backbone atoms for the complexes only, then 50 ps with no restraints). NPT equilibration was then performed at 1 atm and 298 K (50 ps with restraints on nonsolvent heavy atoms, followed by 50 ps with no restraints). Finally, independent 100 ps NPT equilibration runs were carried out for each of 5 independent replicate runs for all systems to provide varied starting conformations. All restraints used a force constant of 10 kcal mol–1 Å–2. Boresch restraints (including force constants) were fit to the independent 100 ps NPT simulations. Production simulations were run for 0.1 ns/ lambda window (see [here](https://github.com/michellab/polaris-poses-challenge-fegrow-a3fe/blob/main/run_a3fe.py) for the lambda schedules). Electrostatics were treated with RF with a 12 A cutoff for neutral ligands. Charged ligands were treated with a co-alchemical ion approach with PME electrostatics with a 10 A cutoff. HMR with a factor of 3 was applied and a 4 fs timestep was used. OpenMM’s LangevinMiddleIntegrator was used (friction coefficient of 1 ps–1, coupled to a heat bath at 298 K), along with a Monte Carlo barostat (1 atm). A Zacharias-type soft-core potential was used. The final ABFE was calculated using the last 80 % of the data using MBAR, and the mean ABFE from the 5 replicate runs was used as the score. The pose with the lowest predicted ABFE was selected. Its coordinates were extracted from after the first NPT equilibration stage with no restraints, but before the independent NPT equilibration. The protocol was generally as given [here](https://pubs.acs.org/doi/full/10.1021/acs.jctc.4c00806) for the non-adaptive runs, other than the reduced run times.

For more details, see

 - `run_a3fe.py`: The script used to set-up and run all ABFE calculations. Provides details of e.g. the lambda schedules.
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
