import sqlalchemy
from geneditid.config import cfg
from geneditid.model import Base
from geneditid.model import Target
import requests
import json


def get_gene_id(gene_symbol, species='homo_sapiens'):
    species = species.lower().replace(' ', '_')
    url = ("https://rest.ensembl.org/xrefs/symbol/{}/{}?").format(species, gene_symbol)
    r = requests.get(url, headers={"Content-Type": "application/json"})
    gene = json.loads(r.text)
    if gene:
        return gene[0]['id']


def main():
    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    dbsession = DBSession()
    targets = dbsession.query(Target).all()
    print('Update gene ids on target in {} to be Ensembl ones starting with ENSG...'.format(cfg['DATABASE_URI']))
    for target in targets:
        old_gene_id = target.gene_id
        new_gene_id = get_gene_id(target.gene_id, target.genome.species)
        if not old_gene_id == new_gene_id:
            new_gene_id = get_gene_id(target.name, target.genome.species)
            if not new_gene_id and target.name.startswith('LEPR'):
                new_gene_id = get_gene_id('LEPR', target.genome.species)

            if new_gene_id:
                if target.description:
                    if len(target.description) >= 984:
                        target.description = "{}... [submitted_gene_id: {}]".format(target.description[:984], old_gene_id)
                    else:
                        target.description = "{} [submitted_gene_id: {}]".format(target.description, old_gene_id)
                else:
                    target.description = "[submitted_gene_id: {}]".format(old_gene_id)
                target.gene_id = new_gene_id
                dbsession.commit()
                print(">>>UPDATED: {}\t{}\t{}\t{}".format(target.name, old_gene_id, new_gene_id, target.description))
            else:
                print("!!!  NO ID: {}\t{}\t{}".format(target.name, target.gene_id, target.description))
        else:
            print("UP-TO-DATE: {}\t{}\t{}".format(target.name, target.gene_id, target.description))
    print('done')


if __name__ == '__main__':
    main()
