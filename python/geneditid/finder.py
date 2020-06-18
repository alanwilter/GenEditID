import os
import logging

from geneditid.config import cfg
from geneditid.connect import dbsession

from geneditid.model import Primer
from geneditid.model import Amplicon
from geneditid.model import Guide
from geneditid.model import Target
from geneditid.model import Project

from Bio.Seq import Seq
from pyfaidx import Fasta


class FinderException(Exception):

    def __init__(self, msg=None):
        if msg is None:
            msg = "Loader error."
        super(FinderException, self).__init__(msg)
        self.message = msg

    def __str__(self):
        return self.message


class AmpliconFinder():

    def __init__(self, dbsession, project_geid):
        self.log = logging.getLogger(__name__)
        self.dbsession = dbsession
        self.project = self.dbsession.query(Project).filter(Project.geid == project_geid).first()
        if not self.project:
            raise LoaderException("Project {} not found".format(project_geid))
        self.log.info('Project {} found'.format(self.project.geid))
        self.config_file = os.path.join(self.project.project_folder, "amplicount_config.csv")

    def get_amplicons(self):
        results = []
        amplicons = self.dbsession.query(Amplicon)\
                                  .filter(Amplicon.project == self.project)\
                                  .all()
        for amplicon in amplicons:
            self.log.info('Amplicon {} retrieved'.format(amplicon.name))
            amplicon_result = {'name': amplicon.name,
                               'refgenome': amplicon.guide.genome.fa_file,
                               'fprimer_seq': amplicon.fprimer.sequence,
                               'rprimer_seq': amplicon.rprimer.sequence,
                               'guide_loc': amplicon.guide_location,
                               'chr': int(amplicon.chromosome),
                               'amplicon_len': amplicon.end-amplicon.start+1,
                               'target_name': amplicon.guide.target.name}
            results.append(amplicon_result)
        return results


    def find_primer(self, sequence, primer_seq):
        primer_seq_ori = primer_seq
        primer_loc = sequence.find(primer_seq)
        if primer_loc == -1:
            primer_seq = str(Seq(primer_seq).reverse_complement())
            primer_loc = sequence.find(primer_seq)
        return primer_loc, primer_seq, primer_seq_ori


    def get_primer_pair(self, sequence, fprimer_seq, rprimer_seq):
        fprimer_loc, fprimer_seq, fprimer_seq_ori = self.find_primer(sequence, fprimer_seq)
        rprimer_loc, rprimer_seq, rprimer_seq_ori = self.find_primer(sequence, rprimer_seq)
        return fprimer_loc, fprimer_seq, rprimer_loc, rprimer_seq


    def find_amplicon_sequence(self, refgenome, amplicon_name, guide_loc, chr, fprimer_seq, rprimer_seq):
        self.log.info("Search amplicon sequence +/- 1000bp around guide location {} on chrom {}".format(guide_loc, chr))
        start = guide_loc - 1000
        end = guide_loc + 1000
        if os.path.exists(refgenome + '.fai'):
            self.log.info('fai file for {} already exists, there is no need to rebuild indexes'.format(refgenome))
            sequence = Fasta(refgenome, rebuild=False, build_index=False, read_ahead=10000)['{}'.format(chr)][start:end].seq
        else:
            self.log.info('fai file for {} do not exist, it will take a while to generate it'.format(refgenome))
            sequence = Fasta(refgenome, read_ahead=10000)['{}'.format(chr)][start:end].seq

        submitted_fprimer_seq = fprimer_seq
        submitted_rprimer_seq = rprimer_seq

        fprimer_loc, fprimer_seq, rprimer_loc, rprimer_seq = self.get_primer_pair(sequence, fprimer_seq, rprimer_seq)

        if fprimer_loc > rprimer_loc:
            fprimer_loc, fprimer_seq, rprimer_loc, rprimer_seq = self.get_primer_pair(sequence, rprimer_seq, fprimer_seq)

        amplicon_seq = sequence[fprimer_loc:(rprimer_loc + len(rprimer_seq))]
        amplicon_start = int(start) + fprimer_loc + 1
        amplicon_end = int(start) + (rprimer_loc + len(rprimer_seq))
        amplicon_coord = "chr{}:{}-{}".format(chr, amplicon_start, amplicon_end)

        msg = ''
        if fprimer_loc == -1 or rprimer_loc == -1:
            raise FinderException('Primers (forward_primer: {}, reverse_primer: {}) not found for amplicon {} (forward_primer_start: {}, reverse_primer_start: {}). Check your primer sequences, or try with a guide location different than {} to search within a different genomic interval than [{}:{}]. If the primer sequence are more than 1,000bp from the cut site each way, they will not be found.'.format(submitted_fprimer_seq.upper(), submitted_rprimer_seq.upper(), amplicon_name, fprimer_loc, rprimer_loc, guide_loc, start, end))
        if not submitted_fprimer_seq == fprimer_seq:
            msg += 'Forward primer sequence different than the one submitted! [submitted: {}, found: {}]\n'.format(submitted_fprimer_seq, fprimer_seq)
        if not submitted_rprimer_seq == rprimer_seq:
            msg += 'Reverse primer sequence different than the one submitted! [submitted: {}, found: {}]\n'.format(submitted_rprimer_seq, rprimer_seq)

        amplicon = {'fprimer_loc': fprimer_loc,
                    'fprimer_seq': fprimer_seq,
                    'rprimer_loc': rprimer_loc,
                    'rprimer_seq': rprimer_seq,
                    'seq': amplicon_seq,
                    'start': amplicon_start,
                    'end': amplicon_end,
                    'coord': amplicon_coord,
                    'info': msg}
        self.log.info('Amplicon sequence {} found'.format(amplicon_seq))
        return amplicon


    def write_amplicount_config_file(self):
        with open(self.config_file, "w") as out:
            out.write("id,fprimer,rprimer,amplicon,coord,info\n")
            found_amplicon_unique_list = []
            for amplicon in self.get_amplicons():
                try:
                    self.log.info('Amplicon {}'.format(amplicon['name']))
                    found_amplicon = self.find_amplicon_sequence(amplicon['refgenome'], amplicon['name'], amplicon['guide_loc'], amplicon['chr'], amplicon['fprimer_seq'], amplicon['rprimer_seq'])
                    # remove duplicated amplicons
                    if found_amplicon:
                        if not found_amplicon['coord'] in found_amplicon_unique_list:
                            found_amplicon_unique_list.append(found_amplicon['coord'])
                            fprimer = found_amplicon['fprimer_seq']
                            rprimer = found_amplicon['rprimer_seq']
                            seq = found_amplicon['seq']
                            out.write("chr{}_{},{},{},{},{},{}\n".format(amplicon['chr'], found_amplicon['start'], fprimer, rprimer, seq, found_amplicon['coord'], found_amplicon['info']))
                except FinderException as e:
                    self.log.error('--- Amplicon #{}'.format(amplicon['name']))
                    self.log.error('Target name\t{}'.format(amplicon['target_name']))
                    self.log.error(e)
                    self.log.error('---')
                    raise e
                except Exception as e:
                    self.log.error('--- Amplicon #{}'.format(amplicon['name']))
                    self.log.error('Target name\t{}'.format(amplicon['target_name']))
                    self.log.error(e)
                    self.log.error('---')
                    raise FinderException('Unexpected error for Amplicon {} on target {}'.format(mplicon['name'], amplicon['target_name']))
