import click
from sys import stdin, stdout, stderr


def get_bc_map(bc_tbl):
    out = {}
    with open(bc_tbl) as bct:
        for line in bct:
            tkns = line.strip().split()
            if len(tkns) != 3:
                continue
            bc, rid, cl = tkns[0], tkns[1], tkns[2]
            ebx = f'{bc}:E:{cl}'
            out[rid] = ebx
    return out


def get_aln_map(sam_file):
    if not sam_file:
        return None
    out = {}
    with open(sam_file) as sf:
        for tkns in (line.strip().split() for line in sf):
            out[tkns[0]] = tkns[2]
    return out

def annotate_rid_line(line, aln_map, bc_map):
    line = line.strip()
    rid = line.split()[0][1:]
    if aln_map:
        try:
            aln = 'AL:' + aln_map[rid]
        except KeyError:
            aln = 'AL:UNAL'
        line += ' ' + aln
    try:
        line += ' ' + bc_map[rid]
    except KeyError:
        pass
    line += '\n'
    return line


@click.command()
@click.option('-s', '--sam-file', default=None)
@click.argument('bc_tbl')
def main(sam_file, bc_tbl):
    aln_map, bc_map = get_aln_map(sam_file), get_bc_map(bc_tbl)
    for i, line in enumerate(stdin):
        if (i % 4) == 0:
            stdout.write(annotate_rid_line(line, aln_map, bc_map))
        else:
            stdout.write(line)
        

if __name__ == '__main__':
    main()
