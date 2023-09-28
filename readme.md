
# Astra: Scalable HMM-based Sequence Search and Retrieval

<img src="img/astra_logo.png" width="50%">

## Overview

Astra is a Python package designed to facilitate bioinformatic workflows, especially those involving Hidden Markov Models (HMMs). It automates the process of downloading, installing, and searching custom or pre-installed HMM databases. Astra aims to streamline bioinformatic analyses, allowing for greater flexibility and ease of use.

## Features

- **Initialize**: Download and install HMM databases directly from various sources with a simple command.
- **Search**: Perform advanced HMM searches on sequence data with customizable options.
- **Flexibility**: Combine both custom and pre-installed HMM databases in a single search run.
- **Configurable**: Customize your workflow with a variety of command-line options.

## Installation

To install the package, use pip (Not yet live! Don't try this! If you do I'll make fun of you because you didn't read!:

```bash
pip install astra
```

## Usage

### Initialization

To initialize and download HMM databases:

```bash
Astra initialize --hmms database_name
```

### Search

To perform an HMM search:

```bash
Astra search --fasta your_fasta_file --hmms database_name
```

### Combined Search

To perform an HMM search using both custom and pre-installed HMM databases:

```bash
Astra search --fasta your_fasta_file --hmms custom_db --installed_hmms pre_installed_db
```

## Contributing

Contributions are welcome! Please read the [contributing guidelines](CONTRIBUTING.md) for more information.
