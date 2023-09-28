# Initialize.py Usage Guide

## Introduction

The `initialize.py` script is designed to manage the installation and querying of Hidden Markov Model (HMM) databases. It provides utilities for configuring the installation directory, downloading databases, and listing available and installed databases.

## Functions

### `load_config()`

- **Purpose**: Loads or initializes the configuration JSON file (`astra_config.json`).
- **Behavior**: 
    - If the configuration file doesn't exist, it creates one with default settings.
    - If a custom path is provided by the user, validates and saves it.
- **Output**: Returns a dictionary containing the configuration settings.

### `load_json()`

- **Purpose**: Loads the `hmm_databases.json` file.
- **Behavior**: If the file doesn't exist in the specified directory, it copies the file from the package directory.
- **Output**: Returns a dictionary containing the database information.

### `install_databases(parsed_json, db_name, db_path)`

- **Purpose**: Installs the specified HMM database.
- **Parameters**: 
    - `parsed_json`: The content of `hmm_databases.json` as a dictionary.
    - `db_name`: The name of the database to install.
    - `db_path`: The directory to install the database into.
- **Behavior**: 
    - Downloads and extracts the database files.
    - Updates the `installed` and `installation_dir` fields in `hmm_databases.json`.

### `show_available_databases(parsed_json)`

- **Purpose**: Lists the available databases that are not yet installed.
- **Parameters**: 
    - `parsed_json`: The content of `hmm_databases.json` as a dictionary.
- **Output**: Prints the list of available databases to the console.

### `show_installed_databases(parsed_json)`

- **Purpose**: Lists the databases that are installed.
- **Parameters**: 
    - `parsed_json`: The content of `hmm_databases.json` as a dictionary.
- **Output**: Prints the list of installed databases to the console.

## Usage Flow

1. Run `load_config()` to set up or load the configuration file.
2. Run `load_json()` to prepare the databases JSON file.
3. Use `install_databases()` to install the desired databases.
4. Utilize `show_available_databases()` and `show_installed_databases()` to view database status.

## Customizing Database Installation

If you want to customize the database installation directory, modify the `db_path` field in `astra_config.json`.

## Note

This guide assumes that you are familiar with Hidden Markov Models, database management, and have a working Python environment.

