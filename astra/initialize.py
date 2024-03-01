import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tarfile
import textwrap
import pyhmmer.plan7
from tqdm import tqdm
import urllib.request
from platformdirs import user_config_dir
import pandas as pd
import requests
from tqdm import tqdm

class TqdmUpTo(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)

def load_config():
    app_name = "Astra"
    app_author = "YourOrg"  # Replace with the actual name of your organization or app author

    # The default directory for the HMM databases.json (inside the package directory)
    package_dir = os.path.dirname(os.path.abspath(__file__))
    default_db_json_path = os.path.join(package_dir, 'hmm_databases.json')

    # Attempt to load the existing hmm_databases.json
    if os.path.exists(default_db_json_path):
        with open(default_db_json_path, 'r') as f:
            hmm_databases = json.load(f)
    else:
        print("hmm_databases.json not found!! Please raise an issue on github.")
        sys.exit()

    # If 'db_path' is empty, prompt the user for the directory to store HMM databases
    if not hmm_databases.get('db_path'):
        # Use platformdirs to get the standard configuration directory
        default_db_path = user_config_dir(app_name, app_author)
        print(f"The default directory for HMM databases is: {default_db_path}")
        user_input = input("Would you like to use the default directory for the HMM databases? [Y/n] ").strip().lower()
        if user_input == 'n':
            hmm_databases['db_path'] = input("Please enter the full path to the desired HMM database directory: ")
        else:
            hmm_databases['db_path'] = default_db_path
        
        # Ensure the HMM database directory exists
        os.makedirs(hmm_databases['db_path'], exist_ok=True)

        # Update the hmm_databases.json file with the new db_path
        with open(default_db_json_path, 'w') as f:
            json.dump(hmm_databases, f)

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

def install_KOFAM():
    
    #Separate function to install KOFAM because we need to manually add bitscore cutoffs to HMM models
    #which drastically reduces size of the output (288M from a 36M metaproteome vs. Pfam's 24M output!!)
    db_name = 'KOFAM'
    parsed_json = load_config()
    config = load_config()
    db_path = config['db_path']    


    for db in parsed_json['db_urls']:
        #Look. I have to organize this as a for loop because of the JSON structure.
        #I know we're only accessing one element. Leave me alone
        #Feel free to tell me alternatives if you read this far and know what to do! 
        #I'm too lazy to ask GPT. This works just fine

        if (db['installed'] == True and db['installation_dir'] != '') and db['name'] == db_name:
            #Database is already installed; just say so.
            print("Database {} already installed.".format(db['name']))
            continue
        if db['name'] == db_name:
            target_folder = os.path.abspath(os.path.join(db_path, db_name))  # Use absolute path
            
            os.makedirs(target_folder, exist_ok=True)
            
            print(f"Downloading {db_name} to {target_folder}...")
            
            # Download the database
            url = db['url']

            file_name = url.split('/')[-1]
            download_path = os.path.join(target_folder, file_name)
            
            if not os.path.exists(download_path):
                # Download the file with progress bar
                with TqdmUpTo(unit='B', unit_scale=True, miniters=1, desc=file_name) as t:  
                    urllib.request.urlretrieve(url, download_path, reporthook=t.update_to)
            else:
                print("KOFAM profiles already detected in download directory. Decompressing...")
                
            # Extract the file if it's a tar archive
            if tarfile.is_tarfile(download_path):
                print(f"Extracting {file_name}...")
                with tarfile.open(download_path, 'r') as tar_ref:
                    tar_ref.extractall(target_folder)
                os.remove(download_path)  # Remove the original tar file
                
                # Check if a single directory was extracted
                extracted_files = os.listdir(target_folder)
                if len(extracted_files) == 1 and os.path.isdir(os.path.join(target_folder, extracted_files[0])):
                    single_dir = os.path.join(target_folder, extracted_files[0])
                    for file_to_move in os.listdir(single_dir):
                        shutil.move(os.path.join(single_dir, file_to_move), target_folder)
                    os.rmdir(single_dir)  # Remove the now-empty directory
            

            # BUG IS HERE COME FIND ME
            # Mark installation as complete in hmm_databases.json
            db['installed'] = True
            db["installation_dir"] = target_folder  # Add installation directory
            # Write changes to the JSON file 
            json_path = os.path.join(db_path, 'hmm_databases.json')  # Use db_path here
            print(json_path)
            with open(json_path, 'w') as f:
                json.dump(parsed_json, f, indent=4)
            print(f"{db_name} successfully downloaded and extracted.")

    ####
    #Now we've installed the HMMs, let's grab profiles.gz.
    ####

    print("Downloading KOFAM thresholds list...")
    ko_list_url = 'https://www.genome.jp/ftp/db/kofam/ko_list.gz'
    file_name = ko_list_url.split('/')[-1]
    ko_list_path = os.path.join(db_path, file_name)
    # Download the file with progress bar
    with TqdmUpTo(unit='B', unit_scale=True, miniters=1, desc=file_name) as t:  
        urllib.request.urlretrieve(ko_list_url, ko_list_path, reporthook=t.update_to)

    print("Decompressing ko_list...")
    with gzip.open(ko_list_path, 'rb') as f_in:
        with open(ko_list_path[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    print("Adding kofam bitscore thresholds to HMM files...")
    #Now parse ko_list as a pandas DF
    ko_list = pd.read_csv(os.path.join(db_path, 'ko_list'), sep='\t')

    # Not all models actually have thresholds; insist that they do
    # All this extraneous code is for the purposes of avoiding a pandas warning. Don't judge me
    # Create a mask for rows where the threshold is not '-'
    mask = ko_list['threshold'] != '-'

    # Use the mask to filter rows and convert the 'threshold' column to float
    ko_list.loc[mask, 'threshold'] = ko_list.loc[mask, 'threshold'].astype(float)


    
    kofam_dir = os.path.join(db_path, 'KOFAM')
    #Iterate only on rows of ko_list as all other HMMs will lack thresholds
    for index, row in ko_list.iterrows():
        threshold = row.threshold
        knum = row.knum
        #Get path to HMM file
        hmm_file = os.path.join(kofam_dir,'{}.hmm'.format(knum))
        #Add the threshold and overwrite the file!!
        add_threshold(hmm_file, threshold)
    return

def add_threshold(hmm_file_path, threshold):
    # Some thresholds in ko_list are specified for HMMs not provided by the package...
    
    if not os.path.exists(hmm_file_path) or threshold == 0.0 or threshold == '-':
        # '-' values shouldn't still exist in this set, but if the threshold is 0
        # or the HMM is specified in KO_list but not provided in the HMM set,
        # let's just ignore it and not add bad thresholds
        return

    with pyhmmer.plan7.HMMFile(hmm_file_path) as hmm_file:
        hmm = hmm_file.read()

    hmm.cutoffs.gathering = threshold, threshold
    hmm.cutoffs.trusted = threshold, threshold
    hmm.cutoffs.noise = threshold, threshold

    with open(hmm_file_path, "wb") as dst:
        hmm.write(dst)
    
    return

def install_databases(db_name, parsed_json=None, db_path=None):

    # Are you trying to install KOFAM? Let's have separate logic for that.
    if db_name == 'KOFAM':
        return install_KOFAM()

    # Did you call this as a function from an external script?
    # Want to model that function call as 'initialize.install_databases(db_name)'
    # So must leave out parsed_json/db_path
    if parsed_json is None:
        parsed_json = load_json()
    if db_path is None:
        config = load_config()
        db_path = config['db_path']

    # Flag to check if we need to update the JSON file
    update_required = False

    for index, db in enumerate(parsed_json['db_urls']):
        if db['name'] == db_name:
            # Check if database is already installed
            if db['installed'] and db['installation_dir']:
                print(f"Database {db_name} already installed.")
                continue

            target_folder = os.path.abspath(os.path.join(db_path, db_name))

            if os.path.exists(target_folder) and os.listdir(target_folder):
                print(f"Folder for {db_name} exists and is not empty. Skipping download.")
                if not db.get('installation_dir'):
                    db['installation_dir'] = target_folder
                    update_required = False
                continue

            # Create target directory if it does not exist
            os.makedirs(target_folder, exist_ok=True)
            
            print(f"Downloading {db_name} to {target_folder}...")
            
            # Download the database
            url = db['url']

            if "github.com" in url:
                if '/blob/' in url:
                    # URL points to a single file, convert to raw URL
                    raw_url = url.replace('github.com', 'raw.githubusercontent.com').replace('/blob/', '/')
                    # Proceed to download the single file using the raw URL
                    file_name = raw_url.split('/')[-1]
                    download_path = os.path.join(target_folder, file_name)
                    with TqdmUpTo(unit='B', unit_scale=True, miniters=1, desc=file_name) as t:  
                        urllib.request.urlretrieve(raw_url, download_path, reporthook=t.update_to)
                else:
                    # Special case for GitHub URLs
                    url = url.rstrip("/")
                    if 'Karthik' in db_name:
                        #Hard-code this one because the directory structure is different compared to the other github repos 
                        repo_api_url = 'https://api.github.com/repos/kanantharaman/metabolic-hmms/contents/'
                    else:
                        repo_api_url = url.replace("github.com", "api.github.com/repos").replace("/tree/master", "/contents")
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
                
                # Download the file with progress bar
                with TqdmUpTo(unit='B', unit_scale=True, miniters=1, desc=file_name) as t:  
                    urllib.request.urlretrieve(url, download_path, reporthook=t.update_to)
                
                # Extract the file if it's a tar archive
                if tarfile.is_tarfile(download_path):
                    print(f"Extracting {file_name}...")
                    with tarfile.open(download_path, 'r') as tar_ref:
                        tar_ref.extractall(target_folder)
                    os.remove(download_path)  # Remove the original tar file
                    
                    # Check if a single directory was extracted
                    extracted_files = os.listdir(target_folder)
                    if len(extracted_files) == 1 and os.path.isdir(os.path.join(target_folder, extracted_files[0])):
                        single_dir = os.path.join(target_folder, extracted_files[0])
                        for file_to_move in os.listdir(single_dir):
                            shutil.move(os.path.join(single_dir, file_to_move), target_folder)
                        os.rmdir(single_dir)  # Remove the now-empty directory

                # Extract the file if it's a gz archive
                elif file_name.endswith('.gz'):
                    print(f"Decompressing {file_name}...")
                    with gzip.open(download_path, 'rb') as f_in:
                        with open(download_path[:-3], 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(download_path)  # Remove the original gz file
                    
                # Check for other file types and add .hmm extension if necessary
                elif not file_name.endswith(('.hmm', '.HMM')):
                    if db_name == 'dbCAN':
                        #dbCAN provides a single file called 'dbCAN-fam-HMMs.txt.v11'
                        new_file_name = file_name + '.hmm'
                        os.rename(download_path, os.path.join(target_folder, new_file_name))
                    else:
                        os.system('rm ' + os.path.join(target_folder, file_name))

            
            print(f"{db_name} successfully downloaded and extracted.")
            parsed_json['db_urls'][index]['installed'] = True
            parsed_json['db_urls'][index]['installation_dir'] = target_folder
            update_required = True
            break  # Break after updating the relevant database


    if update_required:
        print(db_path)
        # Update the JSON file outside the loop to reflect the installed status and installation_dir
        json_path = os.path.join(db_path, 'hmm_databases.json')
        with open(json_path, 'w') as f:
            json.dump(parsed_json, f, indent=4)



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

    # Load configuration to get the database path
    parsed_json = load_config()

    
    
    db_path = parsed_json['db_path']
    
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
                    install_databases(db['name'], parsed_json, db_path)
                    
        # Case for installing all nucleotide databases
        elif db_name == 'all_nuc':
            for db in parsed_json['db_urls']:
                if db['molecule_type'] == 'nucleotide' and not db['installed']:
                    install_databases(db['name'], parsed_json, db_path)
                    
        # Case for installing a specific database
        else:
            install_databases(db_name, parsed_json, db_path)
