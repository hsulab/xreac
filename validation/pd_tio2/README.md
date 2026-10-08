# Pd4/rutile TiO2 input

**No relaxation performed: the complete H/O/Ti/Pd potential is unavailable.**
The relevant publication is Addou et al., ACS Nano 8, 6321-6333 (2014),
DOI:10.1021/nn501817w; SCM lists its potential as `HOTiPd.ff`. The publisher
supplement contains adsorption energies, not the parameter file. No replacement
parameters were inferred. Parsing, export and force checks have not run.

`initial.xyz` contains Ti48O96Pd4 on unreconstructed rutile (110), not the
paper's reconstructed (011)-(2x1) surface. Three stoichiometric trilayers
are used, with 48 bottom atoms fixed. The cell is 12.994 x 11.836 Angstrom,
xy periodic, with 15 Angstrom vacuum per face. Pd4 edges are 2.75 Angstrom
and its base starts 2.3 Angstrom above surface oxygen. Bulk input uses
a=4.594, c=2.959 Angstrom and u(O)=0.305. No hydroxyls, vacancies or
adsorption-site search are included; this input alone validates no potential.

```sh
python examples/pd_tio2_build.py --output validation/runs/pd-tio2/initial.xyz
```

Only the input and this summary are retained. The next prerequisite is the
unchanged complete `HOTiPd.ff` from an authenticated source.
