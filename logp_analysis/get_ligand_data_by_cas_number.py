#example command: python get_ligand_data_by_cas_numper.py 50-00-0
#this is currently designe specifically for me (Ari) by hard-coding the path to my api key, but this can be changed by the use and not have to pass as a command line argument on every single call

import requests
import os,sys

#input cas number as single command line argument
cas = sys.argv[1]

#file with single line only containing valid CAS common chemistry API key, looks something like (elipses added to represent obscured intermediate characters in key) q1pyZq...hDqAPU
api_file = open("/home/ari.ginsparg-umw/cas_api_key.txt", "r")

for line in api_file.readlines():
    API_KEY = line.strip()



def get_cas_data(cas_rn):
    url = f"https://commonchemistry.cas.org/api/detail?cas_rn={cas_rn}"
    
    headers = {
        "X-API-KEY": API_KEY
    }

    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        data = response.json()
        return data
    else:
        print(f"Error {response.status_code}: {response.text}")
        return None

# Example
#cas = "50-00-0"  # formaldehyde
data = get_cas_data(cas)

if data:
    print("Name:", data.get("name"))
    print("SMILES:", data.get("smile"))   # sometimes 'smile' or 'canonicalSmile'
    print("InChI:", data.get("inchi"))