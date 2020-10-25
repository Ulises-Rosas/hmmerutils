# hmmerutils

Requirements:

* python3
* git (to download this repo)

### Installation

```shell
git clone https://github.com/Ulises-Rosas/hmmerutils.git && cd hmmerutils
python setup.py install
```

### Usage

```
hmmerparser -h
```

```
usage: hmmerparser [-h] -t  [...] -f  [...] -s  [-m] [-n]

                           Parse HMMER output
    Example:

        $ hmmerparse -t [hit tables] -f [scaffolds]  -s [species mame]

optional arguments:
  -h, --help            show this help message and exit
  -t  [ ...], --tables  [ ...]
                        Hit tables from HMMER
  -f  [ ...], --scaffolds  [ ...]
                        Scaffolds or genomes sequences previously used with
                        HMMER
  -s , --species        Species name or label to name sequences
  -m , --min            [Optional] Minimum length for each hit [Default: 100]
  -n , --threads        [Optional] number of cpus [Default: 1]
```
