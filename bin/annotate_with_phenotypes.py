import os
import pandas as pd
from helpers import starting_columns
import glob
import sys
from tqdm import tqdm

tqdm.pandas()
phenotype_db_names = ['grasp', 'ebi', 'clinvar', 'phewas', 'finemapping']


def pack(arr):
    return '\t'.join(map(str, arr)) + '\n'


def parse_grasp(filepath):
    phenotypes = {}
    print('Parsing Grasp')
    with open(filepath) as file:
        r = file.readlines()
    for line in tqdm(list(r)):
        a = line.strip('\n').split('\t')
        if a[1] not in phenotypes:
            phenotypes[a[1]] = set()
        rs_id = f"rs{a[0]}"
        # if rs_id not in snps:
        #     continue
        phenotypes[a[1]].add(rs_id)
    return phenotypes


def parse_finemapping(filepath):
    phenotypes = {}
    print('Parsing Finemapping')
    with open(filepath) as file:
        f = file.readlines()
    for k, line in enumerate(tqdm(list(f))):
        if k == 0:
            continue
        a = line.strip('\n').split(',')
        if a[0] not in phenotypes:
            phenotypes[a[0]] = set()
        rs_id = a[2]
        # if rs_id not in snps:
        #     continue
        phenotypes[a[0]].add(rs_id)
    return phenotypes


def parse_ebi(filepath):
    phenotypes = {}
    print('Parsing EBI')
    with open(filepath, 'r') as file:
        f = file.readlines()
    for k, line in enumerate(tqdm(list(f))):
        if k == 0:
            continue
        try:
            a = line.strip('\n').split('\t')
            if a[7] not in phenotypes:
                phenotypes[a[7]] = set()
            rs_id = f"rs{a[23]}"
            # if rs_id not in snps:
            #     continue
            phenotypes[a[7]].add(rs_id)
        except ValueError:
            continue
    return phenotypes


def parse_gtex(qtlfiles, transqtl, snps):

    result = {'trans': {}, 'cis': {}}
    print('Number of cis-eQTL files:', len(qtlfiles))

    for qtlfile in tqdm(list(qtlfiles)):
        with open(qtlfile) as qfile:

            tis = qtlfile[qtlfile.rfind('/') + 1:qtlfile.find('.v8.')]
            print(tis, end=' ')

            for line in qfile:

                if line.startswith('variant_id'):
                    tit = line.strip('\n').split('\t')
                    titlen = len(tit)
                    continue

                a = line.strip('\n').split('\t')
                a = {tit[x]: a[x] for x in range(titlen)}

                chrpos = '_'.join(a['variant_id'].split('_')[:2])

                if chrpos in snps:
                    result['cis'].set_default(chrpos, (set(), set()))[0].add(tis)
                    result['cis'][chrpos][1].add(a['gene_id'])

    print('Starting to parse transqtl')
    with open(transqtl) as trfile:
        for line in trfile:
            if line.startswith('tissue_id'):
                tit = line.strip('\n').split('\t')
                titlen = len(tit)
                continue

            a = line.strip('\n').split('\t')
            a = {tit[x]: a[x] for x in range(titlen)}

            chrpos = '_'.join(a['variant_id'].split('_')[:2])
            tis = a['tissue_id']
            gen = a['gene_id']

            if chrpos in snps:
                result['trans'].set_default(chrpos, (set(), set()))[0].add(tis)
                result['trans'][chrpos][1].add(gen)
    return result
    
    

def parse_phewas(filepath):
    phenotypes = {}
    print('Parsing phewas')
    with open(filepath, 'r') as file:
        f = file.readlines()
    for k, line in enumerate(tqdm(list(f))):
        if k == 0:
            continue
        a = line[line.find('"rs'):].split('",')
        ph = a[1][1:]
        if ph not in phenotypes:
            phenotypes[ph] = set()
        rs_id = f"rs{int(a[0][3:])}"
        # if rs_id not in snps:
        #     continue
        phenotypes[ph].add(rs_id)

    return phenotypes


def parse_clinvar(filepath):
    phenotypes = {}
    print('Parsing ClinVar')
    with open(filepath, 'r') as file:
        f = file.readlines()
    for k, line in enumerate(tqdm(list(f))):
        if k == 0:
            continue
        a = line.strip('\n').split('\t')
        if 'pathogenic' not in a[6].lower() and 'risk factor' not in a[6].lower():
            continue

        for ph in a[13].split(';'):
            if ph in ('not provided', 'not specified'):
                continue
            if ph not in phenotypes:
                phenotypes[ph] = set()
            rs_id = f"rs{a[9]}"
            # if rs_id not in snps:
            #     continue
            phenotypes[ph].add(rs_id)
    return phenotypes

def arr_to_str(arr):
    return '|'.join(sorted([x for x in arr if x is not None]))

def get_phens_by_id(row, all_phenotypes, ids_phenotypes_dict, gtex):
    snp_id = row['ID']
    snp_posid = row.posID   
    return [arr_to_str([ids_phenotypes_dict[y]
                                for y in all_phenotypes.get(snp_id, {}).get(x, [])
                                if y is not None])
                                for x in phenotype_db_names] + [arr_to_str(snps_dict.get(snp_posid, [None, None])[1]) for snps_dict in gtex.values()]


def remove_phen_name_punctuation(phenotype_name):
    return phenotype_name.lower().replace("'", '').replace('_', ' ')


def main(phenotypes_dir, snps_path, out_path):
    print('Reading files')
    snps = pd.read_table(snps_path)
    snps_positions = snps[[*starting_columns, 'fdrp_bh_ref', 'fdrp_bh_alt']]
    del snps
    snp_ids = snps_positions['ID'].tolist()
    snps_positions['posID'] = snps_positions['#chr'] + '_' + snps_positions['end'].astype(str)
    grasp = os.path.join(phenotypes_dir, 'pheno', 'grasp_pheno.tsv')
    ebi = os.path.join(phenotypes_dir, 'pheno', 'gwas_catalog.tsv')
    clinvar = os.path.join(phenotypes_dir, 'pheno', 'variant_summary.txt')
    fm = os.path.join(phenotypes_dir, 'pheno', 'finemapping.csv')
    phewas = os.path.join(phenotypes_dir, 'pheno', 'phewas-catalog.csv')

    qtlfiles = glob.glob(os.path.join(phenotypes_dir, 'eqtl', 'signif', '*.txt'))
    transqtl = os.path.join(phenotypes_dir, 'eqtl', 'GTEx_Analysis_v8_trans_eGenes_fdr05.txt')
    print('Started parsing DBs')
    phenotypes_for_db_list = [parse_grasp(grasp),
                              parse_ebi(ebi),
                              parse_clinvar(clinvar),
                              parse_phewas(phewas),
                              parse_finemapping(fm),
                              ]
    print('Started parsing GTEX')
    #gtex = parse_gtex(qtlfiles, transqtl, snps_positions[snps_positions[['fdrp_bh_ref', 'fdrp_bh_alt']].min(axis=1) <= 0.05].posID)
    print('Parsing finished')
    gtex = {}
    phenotypes_ids_dict = {}
    ids_phenotypes_dict = {}
    phenotype_id = 1

    for db in phenotypes_for_db_list:
        for phenotype in db:
            r_ph = remove_phen_name_punctuation(phenotype)
            if r_ph not in phenotypes_ids_dict:
                phenotypes_ids_dict[r_ph] = phenotype_id
                ids_phenotypes_dict[phenotype_id] = r_ph
                phenotype_id += 1

    all_phenotypes = {}

    for i in tqdm(range(len(phenotypes_for_db_list)), total=len(phenotypes_for_db_list)):
        for phenotype in phenotypes_for_db_list[i]:
            for rs in phenotypes_for_db_list[i][phenotype]:
                if rs not in all_phenotypes:
                    all_phenotypes[rs] = {x: set() for x in phenotype_db_names}
                all_phenotypes[rs][phenotype_db_names[i]].add(
                    phenotypes_ids_dict[remove_phen_name_punctuation(phenotype)])
    
    print('pheno sizes:', len(phenotypes_ids_dict), len(all_phenotypes))

    snps_positions[[*phenotype_db_names, 'QTLgenes_cis', 'QTLgenes_trans']] = snps_positions.progress_apply(
        lambda x: get_phens_by_id(x, all_phenotypes, ids_phenotypes_dict, gtex), axis=1)
    
    snps_positions.to_csv(out_path, sep='\t', index=False)


if __name__ == '__main__':
    main(*sys.argv[1:])