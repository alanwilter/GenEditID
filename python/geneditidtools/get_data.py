import os
import argparse
import pandas
import subprocess
import shutil
import gzip
from pathlib import Path

import geneditid.log as logger
from geneditid.config import cfg


def run_process(log, cmd, dry_run=True):
    if not dry_run:
        process = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out = process.communicate()[0]
        retcode = process.returncode
        log.info("command '{}' executed".format(" ".join(cmd)))
        if retcode == 0:
            log.debug(out)
            return out
        else:
            raise subprocess.CalledProcessError(retcode, cmd, out)
    else:
        log.info("[dry-run] command '{}' executed".format(" ".join(cmd)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--meta", dest="meta", action="store", help="Metadata file listing accession codes for SRA", required=True)
    parser.add_argument("--fqdump", dest="fqdump", action="store", help="Path to fastq_dump tool", default='~/sratoolkit/bin/fastq-dump', required=False)
    options = parser.parse_args()

    log = logger.get_custom_logger(os.path.join(cfg['PROJECTS_FOLDER'], 'get_data.log'))

    fqdump_cmd = options.fqdump.replace('~', str(Path.home()))

    codes = pandas.read_csv(options.meta, sep='\t')
    for i, row in codes.iterrows():
        # create study folder
        study_folder = os.path.join(cfg['PROJECTS_FOLDER'], row.study)
        if not os.path.exists(study_folder):
            os.makedirs(study_folder)
            log.info('Project folder {} created'.format(study_folder))

        log.info('Downloading {}'.format(row.accession))
        run_process(log, [fqdump_cmd, '--split-files', row.accession], False)

        log.info('Compressing {}_1.fastq into {}'.format(row.accession, os.path.join(study_folder, row.filename)))
        with open("{}_1.fastq".format(row.accession), 'rb') as f_in:
            with gzip.open(os.path.join(study_folder, row.filename), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove("{}_1.fastq".format(row.accession))

        log.info('Compressing {}_2.fastq into {}'.format(row.accession, os.path.join(study_folder, row.filename2)))
        with open("{}_2.fastq".format(row.accession), 'rb') as f_in:
            with gzip.open(os.path.join(study_folder, row.filename2), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove("{}_2.fastq".format(row.accession))

if __name__ == '__main__':
    main()
