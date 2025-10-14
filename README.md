# Snippy-NG



[![PyPI - Version](https://img.shields.io/pypi/v/snippy-nextgen.svg)](https://pypi.org/project/snippy-nextgen)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/snippy-nextgen.svg)](https://pypi.org/project/snippy-nextgen)
[![Test Coverage](https://raw.githubusercontent.com/centre-pathogen-genomics/snippy-ng/refs/heads/_xml_coverage_reports/data/tests/badge.svg)](https://app.codecov.io/github/centre-pathogen-genomics/snippy-ng)
![Downloads](https://img.shields.io/github/downloads/centre-pathogen-genomics/snippy-ng/total)
[![Benchmark](https://byob.yarr.is/centre-pathogen-genomics/snippy-ng/benchmark)](https://cpg.org.au/snippy-blog/posts/2024-10-10-snappy-snippy)
-----

🚨 Snippy-NG is under construction and should not replace Snippy 🚨
----

Check out our progress in the [Snippy-NG Development Blog](https://snippy.cpg.org.au/)!

**Table of Contents**

- [Installation](#installation)
- [License](#license)

## Installation

`snippy` is available on [PyPI](https://pypi.org/project/snippy-nextgen/) and can be installed using `pip` (without external dependencies):
```console
pip install snippy-nextgen
```

The complete snippy-ng environment (including all dependencies) can be installed using the `snippy-ng` installer script. This script will install the latest version of `snippy-ng` and all its dependencies.

```console
curl -sSL https://cpg.org.au/snippy-ng/install.sh | bash -s -- --force
```

## Development

To set up a development environment, clone the repository and install `pixi` and `hatch`. Pixi is used to manage external dependencies, and Hatch is used to manage the Python package development.

```console
git clone git@github.com:centre-pathogen-genomics/snippy-ng.git && cd snippy-ng
# install pixi if not already installed
curl -fsSL https://pixi.sh/install.sh | bash
# install hatch if not already installed
pixi global add hatch
```

Activate the pixi environment and launch a hatch shell. THis will install all dependencies and set up the development environment.

```console
pixi shell
hatch shell
```

```console
snippy-ng --help
```

## License

`snippy-ng` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.
