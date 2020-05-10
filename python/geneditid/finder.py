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
            acoord = "chr{}\t{}\t{}\t{}\t{}".format(amplicon.chromosome,
                                                    amplicon.start,
                                                    amplicon.end,
                                                    amplicon.strand,
                                                    "chr{}_{}".format(amplicon.chromosome, amplicon.start))

            tcoord = "chr{}\t{}\t{}\t{}\t{}".format(amplicon.chromosome,
                                                    amplicon.fprimer.end + 1,
                                                    amplicon.rprimer.start - 1,
                                                    amplicon.strand,
                                                    "chr{}_{}".format(amplicon.chromosome, amplicon.start))
            amplicon_result = {'name': amplicon.name,
                               'refgenome': amplicon.guide.genome.fa_file,
                               'gene_id': amplicon.guide.target.gene_id,
                               'fprimer_seq': amplicon.fprimer.sequence,
                               'rprimer_seq': amplicon.rprimer.sequence,
                               'guide_loc': amplicon.guide_location,
                               'strand': amplicon.strand,
                               'chr': int(amplicon.chromosome),
                               'amplicon_coord': acoord,
                               'amplicon_len': amplicon.end-amplicon.start+1,
                               'target_coord': tcoord,
                               'target_len': amplicon.rprimer.start-amplicon.fprimer.end+1,
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


    def find_amplicon_sequence(self, refgenome, guide_loc, chr, strand, fprimer_seq, rprimer_seq):
        self.log.info("Search amplicon sequence +/- 1000bp around guide location {} on chrom {}".format(guide_loc, chr))
        start = guide_loc - 1000
        if os.path.exists(refgenome + '.fai'):
            self.log.info('fai file for {} already exists, there is no need to rebuild indexes'.format(refgenome))
            sequence = Fasta(refgenome, rebuild=False, build_index=False, read_ahead=10000)['{}'.format(chr)][start:guide_loc+1000].seq
        else:
            self.log.info('fai file for {} do not exist, it will take a while to generate it'.format(refgenome))
            sequence = Fasta(refgenome, read_ahead=10000)['{}'.format(chr)][start:guide_loc+1000].seq

        fprimer_loc, fprimer_seq, rprimer_loc, rprimer_seq = self.get_primer_pair(sequence, fprimer_seq, rprimer_seq)

        if fprimer_loc > rprimer_loc:
            fprimer_loc, fprimer_seq, rprimer_loc, rprimer_seq = self.get_primer_pair(sequence, rprimer_seq, fprimer_seq)

        amplicon_seq = sequence[fprimer_loc:(rprimer_loc + len(rprimer_seq))]
        amplicon_start = int(start) + fprimer_loc + 1
        amplicon_end = int(start) + (rprimer_loc + len(rprimer_seq))
        amplicon_desc = "chr{}:{}:{}".format(chr, amplicon_start, amplicon_end)
        acoord = "chr{}\t{}\t{}\t{}\t{}".format(chr,
                                                amplicon_start,
                                                amplicon_end,
                                                strand,
                                                "chr{}_{}".format(chr, amplicon_start))
        genomic_acoord = "chr{}:{}-{}".format(chr, amplicon_start, amplicon_end)

        target_start = int(start) + (fprimer_loc + len(fprimer_seq)) + 1
        target_end = int(start) + rprimer_loc
        tcoord = "chr{}\t{}\t{}\t{}\t{}".format(chr,
                                                target_start,
                                                target_end,
                                                strand,
                                                "chr{}_{}".format(chr, amplicon_start))
        amplicon = {'fprimer_loc': fprimer_loc,
                    'fprimer_seq': fprimer_seq,
                    'rprimer_loc': rprimer_loc,
                    'rprimer_seq': rprimer_seq,
                    'seq': amplicon_seq,
                    'start': amplicon_start,
                    'end': amplicon_end,
                    'desc': amplicon_desc,
                    'coord': acoord,
                    'gcoord': genomic_acoord,
                    'target_coord': tcoord}
        self.log.info('Amplicon sequence {} found'.format(amplicon_seq))
        return amplicon


    def get_status(self, amplicon, found_amplicon):
        self.log.info('--- Amplicon #{}'.format(amplicon['name']))
        self.log.info('Target name\t{}'.format(amplicon['target_name']))
        if found_amplicon:
            self.log.info('Forward primer\t{}'.format(found_amplicon['fprimer_seq']))
            self.log.info('Reverse primer\t{}'.format(found_amplicon['rprimer_seq']))
            self.log.info('target coord\t{}'.format(found_amplicon['target_coord']))
            self.log.info('amplicon coord\t{}'.format(found_amplicon['coord']))
            self.log.info('amplicon desc\t{}'.format(found_amplicon['desc']))
            self.log.info('amplicon seq\t{}'.format(found_amplicon['seq']))
            self.log.info('amplicon seq len\t{}'.format(len(found_amplicon['seq'])))
            status = 'OK'
            if found_amplicon['fprimer_loc'] == -1 or found_amplicon['rprimer_loc'] == -1:
                self.log.warning('>>> Primers not found! [fprimer: {}, rprimer: {}]'.format(found_amplicon['fprimer_loc'], found_amplicon['rprimer_loc']))
                status = 'ERR'
            if not amplicon['amplicon_len'] == len(found_amplicon['seq']):
                self.log.info('>>> Amplicon length different than the one submitted. [submitted: {}, found: {}]'.format(amplicon['amplicon_len'], len(found_amplicon['seq'])))
                status = 'DIFF'
            if not amplicon['fprimer_seq'] == found_amplicon['fprimer_seq']:
                self.log.info('>>> Forward primer sequence different than the one submitted! [submitted: {}, found: {}]'.format(amplicon['fprimer_seq'], found_amplicon['fprimer_seq']))
                status = 'DIFF'
            if not amplicon['rprimer_seq'] == found_amplicon['rprimer_seq']:
                self.log.info('>>> Reverse primer sequence different than the one submitted! [submitted: {}, found: {}]'.format(amplicon['rprimer_seq'], found_amplicon['rprimer_seq']))
                status = 'DIFF'
            if not amplicon['amplicon_coord'] == found_amplicon['coord']:
                self.log.info('>>> Amplicon coord different than the one submitted! [submitted: {}, found: {}]'.format(amplicon['amplicon_coord'], found_amplicon['coord']))
                status = 'DIFF'
            if not amplicon['target_coord'] == found_amplicon['target_coord']:
                self.log.info('>>> Target coord different than the one submitted! [submitted: {}, found: {}]'.format(amplicon['target_coord'], found_amplicon['target_coord']))
                status = 'DIFF'
        self.log.info('---')
        return status


    def write_amplicount_config_file(self):
        with open(self.config_file, "w") as out:
            out.write("id,fprimer,rprimer,amplicon,coord,status\n")
            found_amplicon_unique_list = []
            for amplicon in self.get_amplicons():
                try:
                    self.log.info('Amplicon {}'.format(amplicon['name']))
                    found_amplicon = self.find_amplicon_sequence(amplicon['refgenome'], amplicon['guide_loc'], amplicon['chr'], amplicon['strand'], amplicon['fprimer_seq'], amplicon['rprimer_seq'])
                    status = self.get_status(amplicon, found_amplicon)
                    # remove duplicated amplicons
                    if found_amplicon:
                        if not found_amplicon['desc'] in found_amplicon_unique_list:
                            found_amplicon_unique_list.append(found_amplicon['desc'])
                            fprimer = found_amplicon['fprimer_seq']
                            rprimer = found_amplicon['rprimer_seq']
                            seq = found_amplicon['seq']
                            out.write("chr{}_{},{},{},{},{},{}\n".format(amplicon['chr'], found_amplicon['start'], fprimer, rprimer, seq, found_amplicon['gcoord'], status))
                except Exception as e:
                    self.log.error('--- Amplicon #{}'.format(amplicon['name']))
                    self.log.error('Target name\t{}'.format(amplicon['target_name']))
                    self.log.error(e)
                    self.log.error('---')
                    continue
