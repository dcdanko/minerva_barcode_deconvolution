from sys import stdout
import pandas as pd
import click


NRANKS = 9


def parse_minerva_file(minerva_file, read_to_taxa):
    bx_tbl, ebx_tbl = {}, {}
    with open(minerva_file) as mf:
        for line in mf:
            tkns = line.strip().split('\t')
            bx = tkns[0]
            ebx = tkns[0] + tkns[2]
            try:
                taxa = read_to_taxa[tkns[1]]
            except KeyError:
                taxa = []
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
            taxa = tkns[1].split('|')
            read_to_taxa[read] = taxa
    return read_to_taxa


def promote(taxas):
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
    

def promote_and_count(bx_tbl):
    before, after  = [0] * NRANKS, [0] * NRANKS
    for taxas in bx_tbl.values():
        promoted = promote(taxas)
        for taxa in taxas:
            before[len(taxa)] += 1
        for taxa in promoted:
            after[len(taxa)] += 1
    return before, after
    



@click.command()
@click.argument('minerva_file')
@click.argument('kraken_file')
def main(minerva_file, kraken_file):
    read_to_taxa = parse_kraken_file(kraken_file)
    bx_tbl, ebx_tbl = parse_minerva_file(minerva_file, read_to_taxa)
    inp_counts, bx_counts = promote_and_count(bx_tbl)
    _, ebx_counts = promote_and_count(ebx_tbl)
    counts = pd.concat([inp_counts, bx_counts, ebx_counts])
    stdout.write(counts.to_csv())


if __name__ == '__main__':
    main()
