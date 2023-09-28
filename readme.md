
# Astra: Scalable HMM-based Sequence Search and Retrieval

<img src="img/astra_logo.png" width="50%">
	
## Overview

Astra is a Python package designed to facilitate bioinformatic workflows involving Hidden Markov Models (HMMs). It serves as a wrapper around (PyHMMER)[https://pyhmmer.readthedocs.io/en/stable/index.html] and (HMMER)[http://hmmer.org/], as well as a . It automates the process of downloading, installing, and searching custom or pre-installed HMM databases. Astra aims to streamline bioinformatic analyses, allowing for greater flexibility and ease of use.

## Features

- **Initialize**: Download and install HMM databases directly from various sources with a simple command.
- **Search**: Perform advanced HMM searches on sequence data with customizable options.
- **Flexibility**: Combine both custom and pre-installed HMM databases in a single search run.
- **Configurable**: Customize your workflow with a variety of command-line options.

## Installation

To install the package, use pip (Not yet live! Don't try this!):

```bash
pip install astra
```
Dependencies:

- Python packages (automatically installed via pip):
	- pyhmmer
	- tqdm
	- urllib
	- pandas
	- biopython
- External dependencies:
	- (pullseq)[https://github.com/bcthomas/pullseq]

## Usage

### Initialization

This is not necessary if you have locally installed HMMs. You can specify those without ever running 'Astra initialize'.

If you would like to view the available HMM databases for install:
```bash
Astra initialize --show_available
```


To initialize and download one of these databases:

```bash
Astra initialize --hmms database_name
```

Or to install them all (takes quite a bit of time!):

```bash
Astra initialize --hmms all_prot
```


### Search

To perform an HMM search:

```bash
Astra search --prot_in your_fasta_file --installed_hmms database_name  --cut_ga --outdir example_output
```

#### Combined Search

To perform an HMM search using custom HMM files:

```bash
Astra search --prot_in your_fasta_file --hmm_in custom_db --installed_hmms pre_installed_db --cut_ga --outdir example_output
```

## Contributing & License

Contributions are welcome! Especially if you have database suggestions. Feel free to raise an issue if you'd like to add a database to the installable list. Feature requests will be considered and implemented if I have the time and ability.

---

MIT License

Copyright (c) 2023 Jacob West-Roberts

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.