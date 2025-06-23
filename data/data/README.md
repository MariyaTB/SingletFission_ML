Experimental Calibration
========================
The molecules used for calibration to experimental data are provided in the
`expt_calib.csv` file, containing:

- DOI
- SMILES representation
- experimental S1 energy (eV)
- experimental T1 energy (eV)
- computed S1 energy (eV, BLYP35/def2-TZVP geometry, M06-2X/def2-TZVP)
- computed T1 energy (eV, BLYP35/def2-TZVP geometry, M06-2X/def2-TZVP)

Groups A1-A3
============
Molecules belonging to groups A1-A3 are provided in the `group_A_opt.csv` file
with the results of our calculations, and some data from the CSD, containing:

- CSD identifier
- DOI
- SMILES representation
- computed S1 energy (eV, X-ray geometry, M06-2X/def2-SVP)
- computed S1 oscillator strength (X-ray geometry, M06-2X/def2-SVP)
- computed T1 energy (eV, X-ray geometry, M06-2X/def2-SVP)
- computed S1 energy (eV, BLYP35/def2-TZVP geometry, M06-2X/def2-TZVP)
- computed S1 oscillator strength (BLYP35/def2-TZVP geometry, M06-2X/def2-TZVP)
- computed T1 energy (eV, BLYP35/def2-TZVP geometry, M06-2X/def2-TZVP)
- computed T2 energy (eV, BLYP35/def2-TZVP geometry, M06-2X/def2-TZVP)

The 262 optimised geometries are gathered in the `opt_structs.xyz` file.
Geometries have been optimised at BLYP35/def2-TZVP level with the
[Gaussian16](http://gaussian.com/gaussian16/) software, starting from experimental
X-ray geometries provided within the CSD. To split the geometries in single
files, install [OpenBabel](https://openbabel.org) and run the shell command

`obabel opt_structs.xyz -OCCDC_.xyz --separate -m`

This will create 262 files named `CCDC_N.xyz` containing the geometries. The
CSD identifier is provided as a comment in each `.xyz` entry.
The 262 molecules have been classified in families, and `.smi` files for easy
visualisation are provided in the `families` folder.

Group B
=======
Molecules belonging to group B are provided in the `group_B.csv` file with
the results of our calculations, and some data from the CSD, containing:

- CSD identifier
- DOI
- SMILES representation
- computed S1 energy (eV, X-ray geometry, M06-2X/def2-SVP)
- computed S1 oscillator strength (X-ray geometry, M06-2X/def2-SVP)
- computed T1 energy (eV, X-ray geometry, M06-2X/def2-SVP)

Molecules are listed according to increasing distance from the decision
boundary, so ideally they should be studied in the order provided, checking
the oscillator strength of S1.
