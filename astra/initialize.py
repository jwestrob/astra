import json
import os, sys, requests
import shutil
import tarfile
import gzip
import urllib.request
import argparse
import subprocess
import textwrap
from tqdm import tqdm 

class TqdmUpTo(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)

def load_config():
    script_dir = os.path.dirname(os.path.realpath(__file__))
    default_hmm_db_json = os.path.join(script_dir, 'hmm_databases.json')

    config_path = os.path.join(os.path.expanduser('~'), '.astra', 'astra_config.json')
    if not os.path.exists(config_path):
        print("Welcome to Astra! It looks like this is your first time running the package.")
        print("A configuration file will be created at ~/.astra/astra_config.json.")
        print("By default, HMM databases will be installed in ~/.astra/databases/.")

        user_path = input("Would you like to specify an alternative installation directory? (Press Enter to use the default): ")

        if user_path:
            if not os.path.exists(user_path):
                try:
                    os.makedirs(user_path, exist_ok=True)
                except Exception as e:
                    print(f"Error: Could not create directory at {user_path}. Please make sure the path is valid.")
                    print(str(e))
                    exit(1)
            db_path = user_path
        else:
            db_path = os.path.join(os.path.expanduser('~'), '.astra', 'databases')
        
        # Create the directory if it doesn't exist
        os.makedirs(db_path, exist_ok=True)
        
        # Copy the default hmm_databases.json to the new directory
        shutil.copy(default_hmm_db_json, db_path)

        os.makedirs(os.path.join(os.path.expanduser('~'), '.astra'), exist_ok=True)
        config = {"db_path": db_path}
        
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=4)
        
        with open(config_path, 'r') as f:
            config = json.load(f)
    else:
        with open(config_path, 'r') as f:
            config = json.load(f)
        
        # Ensure that hmm_databases.json exists in the specified directory
        target_hmm_db_json = os.path.join(config['db_path'], 'hmm_databases.json')
        if not os.path.exists(target_hmm_db_json):
            shutil.copy(default_hmm_db_json, config['db_path'])

    return config





def load_json():
    #Loads controller json file storing database location; initializes database directory if this is the first run
    config = load_config()

    #Database path
    db_path = config['db_path']
    json_path = os.path.join(db_path, 'hmm_databases.json')

    with open(json_path, 'r') as f:
        hmm_databases = json.load(f)

    return hmm_databases



def show_available_databases(parsed_json):
    print("Available databases:")
    
    # Group by molecule type
    protein_dbs = [db for db in parsed_json['db_urls'] if not db['installed'] and db['molecule_type'] == 'protein']
    nucleotide_dbs = [db for db in parsed_json['db_urls'] if not db['installed'] and db['molecule_type'] == 'nucleotide']
    
    if protein_dbs:
        print("  Protein Databases:")
        for db in protein_dbs:
            print(f"    - {db['name']}")
            if 'notes' in db:
                wrapped_notes = textwrap.fill(db['notes'], initial_indent='      Notes: ', subsequent_indent='            ')
                print(wrapped_notes)
            if 'citation' in db:
                wrapped_citation = textwrap.fill(db['citation'], initial_indent='      Citation: ', subsequent_indent='               ')
                print(wrapped_citation)
            print('\n')
    
    if nucleotide_dbs:
        print("  Nucleotide Databases:")
        for db in nucleotide_dbs:
            print(f"    - {db['name']}")
            if 'notes' in db:
                wrapped_notes = textwrap.fill(db['notes'], initial_indent='      Notes: ', subsequent_indent='            ')
                print(wrapped_notes)
            if 'citation' in db:
                wrapped_citation = textwrap.fill(db['citation'], initial_indent='      Citation: ', subsequent_indent='               ')
                print(wrapped_citation)
            print('\n')
    print("HMM databases requiring licenses: SUPERFAMILY (https://supfam.mrc-lmb.cam.ac.uk/SUPERFAMILY/models.html), PRISM (https://prism.adapsyn.com/)")

    return


def download_progress_hook(count, block_size, total_size):
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write("\r%2d%%" % percent)
    sys.stdout.flush()

def install_databases(parsed_json=None, db_name, db_path=None):

    #Did you call this as a function from an external script?
    #Want to model that function call as 'intialize.install_databases(db_name)'
    #So must leave out parsed_json/db_path
    if parsed_json == None:
        parsed_json = load_json()
    if db_path = None:
        config = load_config()
        db_path = config['db_path']        


    for db in parsed_json['db_urls']:
        if db['installed'] == True or db['installation_dir'] != '':
            continue
        if db['name'] == db_name:
            target_folder = os.path.abspath(os.path.join(db_path, db_name))  # Use absolute path
            
            # Check if the folder exists and if it's empty
            if os.path.exists(target_folder):
                if len(os.listdir(target_folder)) != 0:
                    print(f"Folder for {db_name} exists and is not empty. Skipping download.")
                    if "installation_dir" not in db:
                        db["installation_dir"] = target_folder
                        # Update the JSON file to reflect the new installation_dir
                        json_path = os.path.join(db_path, 'hmm_databases.json')  
                        with open(json_path, 'w') as f:
                            json.dump(parsed_json, f, indent=4)
                    elif db['installation_dir'] == '':
                        db["installation_dir"] = target_folder
                        # Update the JSON file to reflect the new installation_dir
                        json_path = os.path.join(db_path, 'hmm_databases.json')  
                        with open(json_path, 'w') as f:
                            json.dump(parsed_json, f, indent=4)
                    return
            else:
                os.makedirs(target_folder, exist_ok=True)
            
            print(f"Downloading {db_name} to {target_folder}...")
            
            # Download the database
            url = db['url']
            if "github.com" in url:
                # Special case for GitHub URLs
                url = url.rstrip("/")
                branch_name = "main"  # or get it dynamically
                repo_api_url = url.replace(f"github.com", f"api.github.com/repos").replace(f"/tree/{branch_name}", "/contents")
                response = requests.get(repo_api_url)
                if response.status_code == 200:
                    files = response.json()
                    for file in files:
                        file_name = file['name']
                        if file_name.lower().endswith('.hmm'):  # Only download .hmm or .HMM files
                            file_url = file['download_url']
                            download_path = os.path.join(target_folder, file_name)
                            with TqdmUpTo(unit='B', unit_scale=True, miniters=1, desc=file_name) as t:  
                                urllib.request.urlretrieve(file_url, download_path, reporthook=t.update_to)
                else:
                    print(f"Failed to fetch GitHub directory: {response.status_code}")


            # Check if the URL points to a directory (ends with '/')
            elif url.endswith('/'):
                subprocess.run(["wget", "-r", "-nH", "--cut-dirs=1", "-P", target_folder, url])
            else:
                file_name = url.split('/')[-1]
                download_path = os.path.join(target_folder, file_name)
                
                # Download the file with progress
                with TqdmUpTo(unit='B', unit_scale=True, miniters=1, desc=file_name) as t:  
                    urllib.request.urlretrieve(url, download_path, reporthook=t.update_to)
                
                # Extract the file if it's a tar archive
                if tarfile.is_tarfile(download_path):
                    print(f"Extracting {file_name}...")
                    with tarfile.open(download_path, 'r') as tar_ref:
                        tar_ref.extractall(target_folder)
                    os.remove(download_path)  # Remove the original tar file

                # Extract the file if it's a gz archive
                elif file_name.endswith('.gz'):
                    print(f"Decompressing {file_name}...")
                    with gzip.open(download_path, 'rb') as f_in:
                        with open(download_path[:-3], 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(download_path)  # Remove the original gz file
                            
                # Check for other file types and add .hmm extension if necessary
                elif not file_name.endswith(('.hmm', '.gz', '.tgz')):
                    new_file_name = file_name + '.hmm'
                    os.rename(download_path, os.path.join(target_folder, new_file_name))
            
            print(f"{db_name} successfully downloaded and extracted.")
            db['installed'] = True
            db["installation_dir"] = target_folder  # Add installation directory
            
            # Update the JSON file to reflect the installed status
            json_path = os.path.join(db_path, 'hmm_databases.json')  # Use db_path here
            with open(json_path, 'w') as f:
                json.dump(parsed_json, f, indent=4)
                
            return
    print(f"Database {db_name} not found.")



def show_installed_databases(parsed_json):
    print("Installed databases:")
    
    installed_protein_dbs = [db for db in parsed_json['db_urls'] if db['installed'] and db['molecule_type'] == 'protein']
    installed_nucleotide_dbs = [db for db in parsed_json['db_urls'] if db['installed'] and db['molecule_type'] == 'nucleotide']
    
    if installed_protein_dbs:
        print("  Protein Databases:")
        for db in installed_protein_dbs:
            print(f"    - {db['name']} (Installed in: {db['installation_dir']})")
    
    if installed_nucleotide_dbs:
        print("  Nucleotide Databases:")
        for db in installed_nucleotide_dbs:
            print(f"    - {db['name']} (Installed in: {db['installation_dir']})")



def main(args):
    # Load JSON containing information about available HMM databases
    parsed_json = load_json()
    
    # Load configuration to get the database path
    config = load_config()
    db_path = config['db_path']
    
    # Extract HMM database names from command line arguments
    hmms = args.hmms
    
    # Show available databases if the flag is set
    if args.show_available:
        show_available_databases(parsed_json)
        sys.exit()
        
    # Show installed databases if the flag is set
    if args.show_installed:
        show_installed_databases(parsed_json)
        sys.exit()
        
    # If no HMM database names are provided, show available databases and return
    if hmms is None or len(hmms) == 0:
        show_available_databases(parsed_json)
        return

    # Split the HMM names if multiple databases are provided in a comma-delimited string
    if ',' in hmms:
        hmms = hmms.split(',')
    else:
        hmms = [hmms]

    # Install the databases
    # Iterate through user-provided database names or special keywords for batch installation
    for db_name in hmms:
        # Case for installing all protein databases
        if db_name == 'all_prot':
            for db in parsed_json['db_urls']:
                if db['molecule_type'] == 'protein' and not db['installed']:
                    install_databases(parsed_json, db['name'], db_path)
                    
        # Case for installing all nucleotide databases
        elif db_name == 'all_nuc':
            for db in parsed_json['db_urls']:
                if db['molecule_type'] == 'nucleotide' and not db['installed']:
                    install_databases(parsed_json, db['name'], db_path)
                    
        # Case for installing a specific database
        else:
            install_databases(parsed_json, db_name, db_path)