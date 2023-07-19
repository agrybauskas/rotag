## Dependencies

Dependencies have to be met before starting the installation. Some of the
installation scripts for specific Linux distributions are provided.

```bash
bash dependencies/Ubuntu-20.04/install.sh
```

## Installation
In order to build rotag, the following commands have to be executed in the
top rotag directory.

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

It is recommended to perform tests before using rotag tools.

```
make test
```