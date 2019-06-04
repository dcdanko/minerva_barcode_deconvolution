
import click


def parse_bc_table(bc_assignment):
    tbl = {}
    for tkns in (line.strip().split() for line in bc_assignment):
        tbl[tkns[1]] = f'{tkns[0]}_{tkns[2]}'
    return tbl


@click.command()
@click.argument('bc_assignment', type=click.File('r'))
@click.argument('fastq', type=click.File('r'))
@click.argument('outfile', type=click.File('w'))
def main(bc_assignment, fastq, outfile):
    bc_tbl = parse_bc_table(bc_assignment)
    for i, line in enumerate(fastq):
        if i % 4 == 0:
            tkns = line[1:].strip().split()
            rid, bc = tkns[:2]
            new_bc = bc_tbl.get(rid, bc)
            new_id_line = f'@{rid} {new_bc} {" ".join(tkns[2:])}\n'
            outfile.write(line)
        else:
            outfile.write(line)


if __name__ == '__main__':
    main()
