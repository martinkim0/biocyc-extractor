#!/usr/bin/python

"""
Author: Martin Kim martinkim0
BIOE 134 (Genetic Design Automation), Fall 2021

A command-line script for extracting and imputing data from BioCyc database files
using the BioCyc and PubChem APIs.
"""

import os
import sys
import re
from typing import Dict, List

import pubchempy as pcp
from biocyc import biocyc

class BioCycData:

    input_root = 'inputs'
    output_root = 'outputs'
    chem_links_name = 'compound-links.dat'
    rxn_links_name = 'reaction-links.dat'
    chems_name = 'compounds.dat'
    rxns_name = 'reactions.dat'

    def __init__(self, org: str):
        self.org = org
        self.chem_links = self.load_chem_links()
        self.rxn_links = self.load_rxn_links()

    def load_chem_links(self) -> Dict:
        f = open(
            f'{self.input_root}/{self.org}/{self.chem_links_name}'
        )
        exclude_pattern = re.compile('<[^<>]*>')
        exclude_chars = ['&', ';', '[', ']']

        # Dictionary of form {biocyc_id: {common_name, inchi, smiles}}
        chem_links = {}
        for line in f:
            if line.startswith('#'):
                continue
            entries = line.split('\t')
            id = entries[0].strip()
            inchi = entries[1].strip()
            smiles = entries[2].strip()
            name = entries[3].strip()
            
            # Exclude unwanted characters from common name
            for char in exclude_chars:
                name = name.replace(char, '')
            name = re.sub(exclude_pattern, '', name)
            if name.startswith('a '):
                name = name[2:]
            elif name.startswith('an '):
                name = name[3:]

            values = {
                'name': name,
                'inchi': inchi,
                'smiles': smiles
            }
            chem_links[id] = values

        f.close()
        return chem_links

    def load_rxn_links(self) -> Dict:
        f = open(
            f'{self.input_root}/{self.org}/{self.rxn_links_name}'
        )

        # Dictionary of form {rxn_id: [ec_nums]}
        rxn_links = {}
        for line in f:
            if line.startswith('#'):
                continue
            entries = line.split('\t')
            if entries[1].isspace():
                continue
            id = entries[0].strip()
            ecs = [ec[3:].strip() for ec in entries[1:]]
            rxn_links[id] = ecs
        
        f.close()
        return rxn_links

    def extract_data(self):
        self.chems = self.extract_chems()
        self.rxns = self.extract_rxns()
        self.write_data() 

    def extract_chems(self) -> Dict:
        f = open(
            f'{self.input_root}/{self.org}/{self.chems_name}',
            encoding='ISO-8859-1'
        )
        exclude_pattern = re.compile('<[^<>]*>')
        exclude_chars = ['&', ';', '[', ']']

        chems = {}
        id, name, inchi, smiles = (None, None, None, None)
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith('//'):
                chems[id] = self.impute_chem(
                    id, 
                    name, 
                    inchi, 
                    smiles
                )
                id, name, inchi, smiles = (None, None, None, None)
                continue
            
            entries = line.split(' - ')
            entry_type = entries[0]
            if entry_type == 'UNIQUE-ID':
                id = entries[1].strip()
            elif entry_type == 'COMMON-NAME':
                name = entries[1].strip()
                for char in exclude_chars:
                    name = name.replace(char, '')
                name = re.sub(exclude_pattern, '', name)
                if name.startswith('a '):
                    name = name[2:]
                elif name.startswith('an '):
                    name = name[3:]
            elif entry_type == 'NON-STANDARD-INCHI':
                inchi = entries[1].strip()
            elif entry_type == 'INCHI':
                inchi = entries[1].strip()
            elif entry_type == 'SMILES':
                smiles = entries[1].strip()

        f.close()
        return chems

    def impute_chem(self, id: str, name: str, inchi: str, smiles: str) -> Dict:
        # Check if in compound links
        try:
            entry = self.chem_links[id]
            in_links = True
        except:
            in_links = False
        # Check if in BioCyc database
        try:
            compound = biocyc.get(id)
            in_biocyc = not compound is None
        except:
            in_biocyc = False

        # Pass 1: Use links and BioCyc
        if name is None:
            if in_links:
                name = entry['name']
            elif in_biocyc:
                name = compound.name
        if inchi is None:
            if in_links and entry['inchi'] != '':
                inchi = entry['inchi']
            elif in_biocyc and compound.inchi != '':
                inchi = compound.inchi
        if smiles is None and in_links:
            smiles = entry['smiles']

        # Pass 2: Use PubChem
        if name is None:
            try:
                c = pcp.get_compounds(inchi, 'inchi')
                name = c[0].iupac_name
            except: 
                try:
                    c = pcp.get_compounds(smiles, 'smiles')
                    name = c[0].iupac_name
                except:
                    name = ''
                
        if inchi is None:
            try:
                c = pcp.get_compounds(name, 'name')
                inchi = c[0].inchi
            except:
                try:
                    c = pcp.get_compounds(smiles, 'smiles')
                    inchi = c[0].inchi
                except: 
                    inchi = ''
        
        if smiles is None:
            try:
                c = pcp.get_compounds(name, 'name')
                smiles = c[0].canonical_smiles
            except:
                try:
                    c = pcp.get_compounds(inchi, 'inchi')
                    smiles = c[0].canonical_smiles
                except: 
                    smiles = ''

        return {
            'name': name,
            'inchi': inchi,
            'smiles': smiles
        }

    def extract_rxns(self) -> Dict:
        f = open(
            f'{self.input_root}/{self.org}/{self.rxns_name}',
            encoding='ISO-8859-1'
        )

        rxns = {}
        id, left, right, direction = (None, [], [], None)
        while True:
            line = f.readline()
            if not line: 
                break
            if line.startswith('//'):
                rxns[id] = self.impute_rxn(
                    id,
                    left,
                    right,
                    direction
                )
                id, left, right, direction = (None, [], [], None)
                continue

            entries = line.split(' - ')
            entry_type = entries[0].strip()
            if entry_type == 'UNIQUE-ID':
                id = entries[1].strip()
            elif entry_type == 'LEFT':
                left.append(entries[1].strip())
            elif entry_type == 'RIGHT':
                right.append(entries[1].strip())
            elif entry_type == 'REACTION-DIRECTION':
                dir = entries[1].strip()
                if dir == 'REVERSIBLE':
                    direction = 0
                elif dir.endswith('LEFT-TO-RIGHT'):
                    direction = 1
                elif dir.endswith('RIGHT-TO-LEFT'):
                    direction = 2
                else:
                    direction = -1

        f.close()
        return rxns

    def impute_rxn(self, id: str, left: List, right: List, direction: int) -> Dict:
        reactant_names, product_names = [], []
        for chem_id in left:
            if chem_id in self.chems:
                reactant_names.append(self.chems[chem_id]['name'])
            else:
                reactant_names.append(chem_id)
            
        for chem_id in right:
            if chem_id in self.chems:
                product_names.append(self.chems[chem_id]['name'])
            else:
                product_names.append(chem_id)

        return {
            'reactant_names': reactant_names,
            'product_names': product_names,
            'reactant_ids': left,
            'product_ids': right,
            'direction': direction
        }

    def write_data(self):
        try:
            os.mkdir(f'{self.output_root}/{self.org}')
        except: pass

        chems_file = open(f'{self.output_root}/{self.org}/chemicals.txt', 'w')
        header = 'id' + '\t' + 'name' + '\t' + 'inchi' + '\t' + 'smiles' + '\n'
        chems_file.write(header)

        for id in self.chems:
            values = self.chems[id]
            line = id
            line += '\t' + values['name']
            if not values['inchi'] is None:
                line += '\t' + values['inchi']
            if not values['smiles'] is None:
                line += '\t' + values['smiles']
            line += '\n'
            chems_file.write(line)
        chems_file.close()

        rxns_file = open(f'{self.output_root}/{self.org}/reactions.txt', 'w')
        header = 'reaction' + '\t' + 'reactants' + '\t' + 'products' + '\n'
        rxns_file.write(header)

        for id in self.rxns:
            values = self.rxns[id]
            line = " + ".join(values['reactant_names'])
            if values['direction'] == 0:
                line += ' <=> '
            elif values['direction'] == 1:
                line += ' -> '
            else:
                line += ' <- '
            line += " + ".join(values['product_names'])
            line += '\t' + " ".join(values['reactant_ids'])
            line += '\t' + " ".join(values['product_ids'])
            line += '\n'
            rxns_file.write(line)
        rxns_file.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print('Incorrect number of arguments')
        sys.exit(1)
    
    root = 'inputs'
    org = sys.argv[1]
    if not os.path.isdir(f'{root}/{org}'):
        print('Invalid path')
        sys.exit(1)

    for fname in os.listdir(f'{root}/{org}'):
        if not os.path.isfile(f'{root}/{org}/{fname}'):
            print(f'Cannot parse path {fname}')
            sys.exit(1)

    # Must change this line for different organisms
    biocyc.set_organism('ECOLI')

    extractor = BioCycData(org)
    extractor.extract_data()

    print(f'Data written to outputs/{org}')