from sys import stdout
import pandas as pd
import click
from sys import stderr

NRANKS = 20


def parse_minerva_file(minerva_file, read_to_taxa):
    bx_tbl, ebx_tbl = {}, {}
    with open(minerva_file) as mf:
        for line in mf:
            tkns = line.strip().split('\t')
            bx = tkns[0]
            ebx = tkns[0] + '$' + tkns[2]
            try:
                taxa = read_to_taxa[tkns[1]]
            except KeyError:
                taxa = ['unclassified']
            if len(taxa) > 8:
                taxa = taxa[:8]
            try:
                bx_tbl[bx].append(taxa)
            except KeyError:
                bx_tbl[bx] = [taxa]
            try:
                ebx_tbl[ebx].append(taxa)
            except KeyError:
                ebx_tbl[ebx] = [taxa]
    return bx_tbl, ebx_tbl


def parse_kraken_file(kraken_file):
    read_to_taxa = {}
    with open(kraken_file) as kf:
        for line in kf:
            tkns = line.strip().split('\t')
            read = tkns[0]
            taxa = tkns[1].split(';')
            read_to_taxa[read] = taxa
    return read_to_taxa


def promote(taxas, taxa_tree=None):
    if taxa_tree is None:
        taxa_tree = build_taxa_tree(taxas)
    promoted_taxas = []
    for taxa in taxas:
        root = taxa_tree
        promoted = []
        for rank in taxa:
            root = root[rank]
            promoted.append(rank)
        while len(root) == 1:
            rank = list(root.keys())[0]
            root = root[rank]
            promoted.append(rank)
        promoted_taxas.append(promoted)
    return promoted_taxas
    

def build_taxa_tree(taxas):
    taxa_tree = {}
    for taxa in taxas:
        root = taxa_tree
        for rank in taxa:
            try:
                root = root[rank]
            except KeyError:
                root[rank] = {}
                root = root[rank]
    return taxa_tree
    

def promote_and_count(bx_tbl):
    before, after  = {}, {}
    for taxas in bx_tbl.values():
        promoted = promote(taxas)
        for taxa in taxas:
            try:
                before[';'.join(taxa)] += 1
            except KeyError:
                before[';'.join(taxa)] = 1
        for taxa in promoted:
            try:
                after[';'.join(taxa)] += 1
            except KeyError:
                after[';'.join(taxa)] = 1
    return before, after


def promote_and_count2_enhanced(bx_tbl, ebx_tbl):
    tbl = {}
    for ebx, taxas in ebx_tbl.items():
        promoted_once = promote(taxas)
        bx_tree = build_taxa_tree(bx_tbl[ebx.split('$')[0]])
        promoted_twice = promote(promoted_once, bx_tree)

        for original_taxa, promoted_taxa in zip(taxas, promoted_twice):
            otaxa, ptaxa = ';'.join(original_taxa), ';'.join(promoted_taxa)
            key = otaxa + ' -> ' + ptaxa 
            try:
                tbl[key] += 1
            except KeyError:
                tbl[key] = 1
    return tbl


def promote_and_count2(bx_tbl):
    tbl = {}
    for taxas in bx_tbl.values():
        promoted = promote(taxas)
        for original_taxa, promoted_taxa in zip(taxas, promoted):
            otaxa, ptaxa = ';'.join(original_taxa), ';'.join(promoted_taxa)
            key = otaxa + ' -> ' + ptaxa 
            try:
                tbl[key] += 1
            except KeyError:
                tbl[key] = 1
    return tbl




@click.command()
@click.argument('minerva_file')
@click.argument('kraken_file')
def main(minerva_file, kraken_file):
    read_to_taxa = parse_kraken_file(kraken_file)
    bx_tbl, ebx_tbl = parse_minerva_file(minerva_file, read_to_taxa)
    '''
    inp_counts, bx_counts = promote_and_count(bx_tbl)
    _, ebx_counts = promote_and_count(ebx_tbl)
    counts = pd.DataFrame.from_dict({'inputs': inp_counts, 'standard': bx_counts, 'enhanced': ebx_counts}, orient='columns')
    '''
    counts = pd.DataFrame.from_dict({
        'standard': promote_and_count2(bx_tbl),
        'enhanced': promote_and_count2_enhanced(bx_tbl, ebx_tbl),
    }, orient='columns')
    stdout.write(counts.to_csv())


if __name__ == '__main__':
    main()
