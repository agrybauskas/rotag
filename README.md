# rotag

*rotag* &ndash; a collection of tools that generate and analyse side-chain rotamer
libraries from protein structure data.

## Dependencies

Dependencies have to be met before starting the installation. Some of the
installation scripts for specific Linux distributions are provided.

```bash
bash dependencies/Ubuntu-22.04/install.sh
```

## Installation
In order to build *rotag*, the following commands have to be executed in the
top *rotag* directory.

```bash
make
```

After the installation, the environment has to be set. The commands can be placed
in `~/.bashrc` file.

```bash
ROTAG_SRC=~/src/rotag
export PERL5LIB=${ROTAG_SRC}/lib:${PERL5LIB}
export PATH=${ROTAG_SRC}/scripts:${PATH}
```

It is recommended to perform tests before using *rotag* tools.

```
make test
```

## Usage

To select specific side-chain atoms, ```rotag_select``` script is used.

```bash
rotag_select -t 'resid 20 && chain A' -s 'target around 5' -k tests/inputs/4dhw.cif > 4dhw-rotag-select.cif
```

The rotations around the dihedral angles of selected residues are performed with ```rotag_scan```.

```bash
rotag_scan 4dhw-rotag-select.cif > 4dhw-rotag-scan.cif
```

There is a way to control the dihedral angle steps with ```-a``` option:

```bash
rotag_scan -a 'THR: -180.0..36.0..180.0' 4dhw-rotag-select.cif > 4dhw-rotag-scan.cif
```

If too many dihedral angle choices are suggested, there are some filtering capabilities.

```bash
rotag_library -M 0.1 4dhw-rotag-scan.cif > 4dhw-rotag-library.cif
```

The outputs are stored in separate CIF data blocks in order not to interfere with geometry data. The addition of atoms has to be explicitly stated.

```bash
rotag_add -S 4dhw-rotag-library.cif > 4dhw-rotag-add.cif
```

One more important factor when calculating side-chain rotations is hydrogens. Sometimes they are missing in the structural data. It could be added with the ```-H``` option in ```rotag_add``` script.

```bash
rotag_add -H -c -s 4dhw-rotag-select.cif > 4dhw-rotag-select-H.cif
```

However, when doing more precise and sensitive calculations, we suggest using external tools for this application.

The whole workflow can be streamlined with ```STDIN``` and ```STDOUT``` using pipes.

```bash
rotag_select -t 'resid 20 && chain A' -s 'target around 5' -k tests/inputs/4dhw.cif \
    | rotag_scan \
    | rotag_library -M 0.1 \
    | rotag_add -S \
    > 4dhw-rotag-add.cif
```

For more detailed information about scripts, each script has ```--help``` argument.

## Citing

If you use *rotag* in your research, please cite:
* Algirdas Grybauskas, Saulius Gražulis, Building protein structure-specific rotamer libraries, Bioinformatics, Volume 39, Issue 7, July 2023, btad429. https://doi.org/10.1093/bioinformatics/btad429
