###
# configuration that should be handled manually and will differ between systems
###

from os.path import expanduser, splitext
from os import getenv
from pandas import read_table

HOME = expanduser('~')
BASE = HOME + '/sleuthALR_paper_analysis'
CONDA = HOME + '/miniconda3'

N_THREADS = 50

###
# software
###
BIN = CONDA + '/bin'

KALLISTO = BIN + '/kallisto'
SALMON = BIN + '/salmon'

UPDATED_PATH = 'PATH=' + BIN + ':$PATH'

SAMTOOLS = BIN + '/samtools'
PANDOC = BIN + '/pandoc'
OPENSSL = BIN + '/openssl'
FASTQ_DUMP = BIN + '/fastq-dump'

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
SPIKEIN_FA = BASE + '/annotation/' + TRANSCRIPTOME_NAME + '.andErcc92.fa'

KALLISTO_INDEX = BASE + '/index/' + TRANSCRIPTOME_NAME + '.kidx'
SALMON_INDEX = BASE + '/index/' + TRANSCRIPTOME_NAME + '.sidx'
KALLISTO_SPIKEIN_INDEX = BASE + '/index/' + TRANSCRIPTOME_NAME + '.andErcc92.kidx'
SALMON_SPIKEIN_INDEX = BASE + '/index/' + TRANSCRIPTOME_NAME + '.andErcc92.sidx'

## Mouse annotations
MOUSE_NAME = 'Mus_musculus.gencode.vM12'
MOUSE_TRANSCRIPTOME_FA = BASE + '/annotation/' + MOUSE_NAME + '.fa'
MOUSE_KALLISTO_INDEX = BASE + '/index/' + MOUSE_NAME + '.kidx'
MOUSE_SALMON_INDEX = BASE + '/index/' + MOUSE_NAME + '.sidx'

## Yeast annotations
YEAST_NAME = 'ASM294v2.pombase.all'
YEAST_ANNOS = BASE + '/annotation/' + YEAST_NAME + '_annos.txt'
YEAST_TRANSCRIPTOME_FA = BASE + '/annotation/' + YEAST_NAME + '.fa'

YEAST_KALLISTO_INDEX = BASE + '/index/' + YEAST_NAME + '.kidx'
YEAST_SALMON_INDEX = BASE + '/index/' + YEAST_NAME + '.sidx'

ERCC_FILE = 'ERCC92_data.txt'
ERCC_FA = 'ERCC92.fa'

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
