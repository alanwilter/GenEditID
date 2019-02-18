import sqlalchemy
from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import Primer
from dnascissors.model import AmpliconSelection
from dnascissors.model import Guide
from dnascissors.model import Target
from dnascissors.model import Project
import requests
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def find_primer_sequences(geid):
    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    dbsession = DBSession()
    amplicons = []
    amplicon_selections = dbsession.query(AmpliconSelection)\
                                   .join(AmpliconSelection.guide)\
                                   .join(Guide.target)\
                                   .join(Target.project)\
                                   .filter(Project.geid == geid)\
                                   .all()
    primers = dbsession.query(Primer).all()
    for amplicon_selection in amplicon_selections:
        amplicon = amplicon_selection.amplicon
        aprimers = [p for p in primers if amplicon in p.amplicons]
        for primer in aprimers:
            if primer.strand == 'forward':
                fprimer_seq = primer.sequence
                tstart = primer.end + 1
            elif primer.strand == 'reverse':
                rprimer_seq = primer.sequence
                tend = primer.start - 1
        strand = '+'
        if amplicon_selection.guide_strand == 'reverse':
            strand = '-'
        acoord = "chr{}\t{}\t{}\t{}\t{}".format(amplicon.chromosome,
                                                amplicon.start,
                                                amplicon.end,
                                                strand,
                                                "{}_chr{}_{}".format(amplicon.genome.assembly, amplicon.chromosome, amplicon.start))

        tcoord = "chr{}\t{}\t{}\t{}\t{}".format(amplicon.chromosome,
                                                tstart,
                                                tend,
                                                strand,
                                                "{}_chr{}_{}".format(amplicon.genome.assembly, amplicon.chromosome, tstart))
        amplicons.append({'gene_id': amplicon_selection.guide.target.gene_id,
                          'fprimer_seq': fprimer_seq,
                          'rprimer_seq': rprimer_seq,
                          'amplicon_coord': acoord,
                          'target_coord': tcoord,
                          'target_len': tend-tstart+1})
    return amplicons


def find_primer_location(sequence, primer_seq):
    primer_seq_ori = primer_seq
    primer_loc = sequence.find(primer_seq)
    if primer_loc == -1:
        primer_seq = str(Seq(primer_seq, IUPAC.unambiguous_dna).reverse_complement())
        primer_loc = sequence.find(primer_seq)
    return primer_loc, primer_seq, primer_seq_ori


def get_primer_locations(sequence, fprimer_seq, rprimer_seq):
    fprimer_loc, fprimer_seq, fprimer_seq_ori = find_primer_location(sequence, fprimer_seq)
    rprimer_loc, rprimer_seq, rprimer_seq_ori = find_primer_location(sequence, rprimer_seq)
    return fprimer_loc, fprimer_seq, rprimer_loc, rprimer_seq


def get_target_sequence(ensembl_gene_id, fprimer_seq, rprimer_seq):
    url = ("http://rest.ensembl.org/sequence/id/{}?").format(ensembl_gene_id)
    r = requests.get(url, headers={"Content-Type": "text/plain"})
    sequence = r.text

    fprimer_loc, fprimer_seq, rprimer_loc, rprimer_seq = get_primer_locations(sequence, fprimer_seq, rprimer_seq)
    target_seq = sequence[fprimer_loc+len(fprimer_seq):rprimer_loc-len(rprimer_seq)]

    if len(target_seq) == 0:
        fprimer_loc, fprimer_seq, rprimer_loc, rprimer_seq = get_primer_locations(sequence, rprimer_seq, fprimer_seq)
        target_seq = sequence[fprimer_loc+len(fprimer_seq):rprimer_loc-len(rprimer_seq)]

    return fprimer_loc, fprimer_seq, rprimer_loc, rprimer_seq, target_seq


def main():
    geid = 'GEP00013'
    for amplicon in find_primer_sequences(geid):

        ensembl_gene_id = amplicon['gene_id']
        fprimer_seq = amplicon['fprimer_seq']
        rprimer_seq = amplicon['rprimer_seq']
        fprimer_loc, fprimer_seq, rprimer_loc, rprimer_seq, wt_seq = get_target_sequence(ensembl_gene_id, fprimer_seq, rprimer_seq)

        print(geid)
        print('Ensembl Gene ID\t{}'.format(ensembl_gene_id))
        print('forward primer\t{} {}'.format(fprimer_seq, fprimer_loc))
        print('reverse primer\t{} {}'.format(rprimer_seq, rprimer_loc))
        if fprimer_loc == -1 or rprimer_loc == -1:
            print('Primers not found in Ensembl Gene ID {}!'.format(ensembl_gene_id))
        else:
            print('target seq\t{}'.format(wt_seq))
            print('target seq len\t{}'.format(len(wt_seq)))
            print('amplicon coord\t{}'.format(amplicon['amplicon_coord']))
            print('target coord\t{}'.format(amplicon['target_coord']))
            print('target len\t{}'.format(amplicon['target_len']))
            if not len(wt_seq) == amplicon['target_len']:
                print('Target coordinates not equal. Do run blat!')


if __name__ == '__main__':
    main()
