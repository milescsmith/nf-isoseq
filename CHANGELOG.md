# nf-core/isoseq: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## 1.0.0 - [2021/06/22]

### `Changed

* All containers updated to newest versions
* Data now saved based on barcode name
* Added several optional parameters (such as log levels, whether to polish, to chunk, etc...)
* Removed all container/conda definitions to the config file
* Using a conda env for isoseq3 because the current biocontainers container has a bug when used by singularity

## 0.2.1 - [2020/07/13]

### `Changed`

* Update container versions.  Switch to using milescsmith/sqanti3-1.4.0

## 0.2.0 - [2020/06/29]

### `Added`

* Added SQANTI3 processing and reporting

### `Changed`

* Switched many of the conda environment definitions to the docker image
  provided by Bioconda


## 0.1.0 - [2020/06/01]

### `Added`

* Initial working release of nf-core/isoseq.


### `Fixed`

### `Dependencies`

### `Deprecated`
