###
# configuration that should be handled manually and will differ between systems
###

# BASE = '/home/hjp/sleuth_paper_analysis'
# BASE = '/Users/hjp/analysis/sleuth_paper'
# TODO: make this a absolute path before submission
from os.path import expanduser, splitext
from os import getenv
from pandas import read_table

HOME = expanduser('~')
BASE = HOME + '/sleuth_paper_analysis'

N_THREADS = 50

###
# software
###
BIN = BASE + '/software/bin'

KALLISTO = BIN + '/kallisto'
SALMON = BIN + '/salmon'

# UPDATED_PATH = 'PATH=' + ':'.join([
#     RSEM_PATH,
#     getenv('PATH')
#     ])
UPDATED_PATH = 'PATH=' + BIN + ':$PATH'

SAMTOOLS = BIN + '/samtools'

PANDOC = BIN + '/pandoc'

# import os

# os.environ['PATH'] = RSEM_PATH + ':' + os.environ['PATH']
# print(RSEM_PATH + ':' + os.environ['PATH'])
# os.environ['PATH'] = RSEM_PATH + ':' + os.environ['PATH']
# print(type(os.environ['PATH']))

###
# simulation wildcards
###
SIM_NAMES = ['run' + format(i, '02d') for i in range(1,31)]
SIM_IDS = ['sample_' + format(i, '02d') for i in range(1,11)]

###
# annotations
###
TRANSCRIPTOME_NAME = 'Homo_sapiens.gencode.v25'
TRANSCRIPTOME_FA = BASE + '/annotation/' + TRANSCRIPTOME_NAME + '.fa'
TRANSCRIPTOME_GTF = BASE + '/annotation/' + TRANSCRIPTOME_NAME + '.gtf'
SPIKEIN_FA = BASE + '/annotation/' + TRANSCRIPTOME_NAME + '.andErcc92.fa'

# NOTE: bowtie creates files with this prefix
BOWTIE_INDEX = BASE + '/index/' + TRANSCRIPTOME_NAME

KALLISTO_INDEX = BASE + '/index/' + TRANSCRIPTOME_NAME + '.kidx'
SALMON_INDEX = BASE + '/index/' + TRANSCRIPTOME_NAME + '.sidx'
KALLISTO_SPIKEIN_INDEX = BASE + '/index/' + TRANSCRIPTOME_NAME + '.andErcc92.kidx'
SALMON_SPIKEIN_INDEX = BASE + '/index/' + TRANSCRIPTOME_NAME + '.andErcc92.sidx'

GENOME_NAME = 'Homo_sapiens.GRCh38.dna.primary_assembly.rel87'
GENOME_FA = BASE + '/annotation/' + GENOME_NAME + '.fa'

# STAR_DIRECTORY = BASE + '/index/star_' + TRANSCRIPTOME_NAME
# STAR_INDEX = STAR_DIRECTORY + '/Genome'
HISAT_INDEX = BASE + '/index/' + GENOME_NAME



RSEM_ANNOTATION_DIR = '/'.join([
    BASE,
    'annotation',
    TRANSCRIPTOME_NAME + '_rsem'])

MOUSE_NAME = 'Mus_musculus.gencode.vM12'
MOUSE_TRANSCRIPTOME_FA = BASE + '/annotation/' + MOUSE_NAME + '.fa'
MOUSE_TRANSCRIPTOME_GTF = BASE + '/annotation/' + MOUSE_NAME + '.gtf'

MOUSE_KALLISTO_INDEX = BASE + '/index/' + MOUSE_NAME + '.kidx'
MOUSE_SALMON_INDEX = BASE + '/index/' + MOUSE_NAME + '.sidx'

MOUSE_GENOME_NAME = 'Mus_musculus.GRCm38.dna.primary_assembly.rel87'
MOUSE_GENOME_FA = BASE + '/annotation/' + MOUSE_GENOME_NAME + '.fa'

MOUSE_HISAT_INDEX = BASE + '/index/' + MOUSE_GENOME_NAME

YEAST_NAME = 'ASM294v2.pombase.all'
YEAST_ANNOS = BASE + '/annotation/' + YEAST_NAME + '_annos.txt'
YEAST_TRANSCRIPTOME_FA = BASE + '/annotation/' + YEAST_NAME + '.fa'

YEAST_KALLISTO_INDEX = BASE + '/index/' + YEAST_NAME + '.kidx'
YEAST_SALMON_INDEX = BASE + '/index/' + YEAST_NAME + '.sidx'

ERCC_FILE = 'ERCC92_data.txt'

###
# functions
###
def source_r(base, fname):
    return 'OMP_NUM_THREADS=1 Rscript --vanilla --default-packages=methods,stats,utils -e \'setwd("{0}")\' -e \'source("{1}")\''.format(base, fname)

def source_rmd(base, file_name, output_name = None):
    if output_name is None:
        output_name = splitext(file_name)[0]
        output_name += '.html'
    return UPDATED_PATH + ' OMP_NUM_THREADS=1 Rscript --vanilla --default-packages=methods,stats,utils,knitr -e \'setwd("{0}")\' -e \'rmarkdown::render("{1}", output_file = "{2}")\''.format(base, file_name, output_name)

def get_acc_sample_dict(fname):
    ret_pd = read_table(fname)
    ret_pd.index = ret_pd.name.tolist()
    ret = ret_pd.drop(columns = 'name')
    RET_DICT = ret.T.to_dict()
    return RET_DICT
