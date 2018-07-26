import click
from math import log2
from sys import stderr, stdin


def aln_purity(aln_counts, size):
    return max(aln_counts) / size


def aln_entropy(aln_counts, size):
    H = 0
    for count in aln_counts:
        p = count / size
        H += p * log2(p)
    if H < 0:
        H = -H
    return H


def eval_cluster(alns):
    size = len(alns) / 2
    counts = {}
    for aln in alns:
        try:
            counts[aln] += 0.5
        except KeyError:
            counts[aln] = 0.5
    aln_counts = list(counts.values())
    return size, aln_purity(aln_counts, size), aln_entropy(aln_counts, size)


def parse_annotated_fastq(fqf, count_all=False):
    for i, line in enumerate(fqf):
        if (i % 4) == 0:
            tkns = line.strip().split()
            al_token = None
            bc_token = None
            for tkn in tkns:
                if ('AL' in tkn) and ('UNAL' not in tkn):
                    al_token = tkn.split(':')[1]
                if ('BX' in tkn) and (count_all or 'E:' in tkn):
                    bc_token = tkn
                    if count_all:
                        bc_token = bc_token.split('E:')[0]
            if al_token and bc_token:
                yield bc_token, al_token

                
def parse_clusters(fastq_filename, count_all=False):
    current_ebx = None
    current_alns = []
    for ebx, aln in parse_annotated_fastq(fastq_filename, count_all=count_all):
        if ebx == current_ebx:
            current_alns.append(aln)
        else:
            if current_ebx:
                yield current_alns
            current_alns = [aln]
            current_ebx = ebx
    if current_ebx:
        yield current_alns


@click.command()
@click.option('--count-all/--count-enhanced', default=False)
def main(count_all):
    total_size, total_purity, total_entropy, N = 0, 0, 0, 0
    print('size,purity,entropy')
    for cluster in parse_clusters(stdin, count_all=count_all):
        size, purity, entropy = eval_cluster(cluster)
        print(f'{size},{purity},{entropy}')
        total_size += size
        total_purity += purity
        total_entropy += entropy
        N += 1
    ave_size = total_size / N
    ave_purity = total_purity / N
    ave_entropy = total_entropy / N
    print('Ave Size, Ave Purity, Ave Entropy, N', file=stderr)
    print(f'{ave_size}, {ave_purity}, {ave_entropy}, {N}', file=stderr)
    



if __name__ == '__main__':
    main()
