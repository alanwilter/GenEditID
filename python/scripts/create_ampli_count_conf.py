import sqlalchemy
from dnascissors.config import cfg
from dnascissors.model import Base
import scripts.log as logger
import scripts.run_ampli_find
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def create_files(session, refgenome, project):
    amplicount_file = "amplicount_config.csv"
    with open(amplicount_file, "w") as amplicount_output:
        amplicount_output.write("id,fprimer,rprimer,amplicon,coord\n")
        amplifind_amplicon_desc_list = []
        i = 0
        for amplicon in scripts.run_ampli_find.get_amplicons(session, project):
            i += 1
            try:
                amplifind_amplicon = scripts.run_ampli_find.find_amplicon_sequence(refgenome, amplicon['guide_loc'], amplicon['chr'], amplicon['strand'], amplicon['fprimer_seq'], amplicon['rprimer_seq'])
                scripts.run_ampli_find.print_amplifind_report(i, amplicon, amplifind_amplicon)
                # remove duplicated amplicons
                if amplifind_amplicon:
                    if not amplifind_amplicon['desc'] in amplifind_amplicon_desc_list:
                        amplifind_amplicon_desc_list.append(amplifind_amplicon['desc'])
                        fprimer = amplifind_amplicon['fprimer_seq']
                        rprimer = amplifind_amplicon['rprimer_seq']
                        seq = amplifind_amplicon['seq']
                        if (project == 'GEP00001') and (amplicon['chr'] == 17):
                            fprimer_ori = fprimer
                            fprimer = str(Seq(rprimer, IUPAC.unambiguous_dna).reverse_complement())
                            rprimer = str(Seq(fprimer_ori, IUPAC.unambiguous_dna).reverse_complement())
                            seq = str(Seq(seq, IUPAC.unambiguous_dna).reverse_complement())
                        amplicount_output.write("chr{}_{},{},{},{},{}\n".format(amplicon['chr'], amplifind_amplicon['start'], fprimer, rprimer, seq, amplifind_amplicon['gcoord']))
            except Exception as e:
                print('--- Amplicon #{}'.format(i))
                print('Target name\t{}'.format(amplicon['target_name']))
                print(e)
                print('---')
                continue


def main():
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--project", dest="project", action="store", help="The project id.", required=True)
    parser.add_argument("--genome", dest="refgenome", action="store", help="The reference genome fasta file e.g. 'Homo_sapiens.GRCh38.dna.toplevel.fa'", required=True)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(os.path.dirname(__file__), 'create_config.log'))

    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()

    try:
        create_files(session, options.refgenome, options.project)
        session.commit()
    except Exception as e:
        log.exception(e)
        session.rollback()
    finally:
        session.close()


if __name__ == '__main__':
    main()
