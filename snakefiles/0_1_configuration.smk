import yaml

# --------------------------------- Config parsing -------------------------- #

# Checks, if config file is specified through --configfile command. If not,
# uses config.yaml file.
if not config:
    with open("config.yaml") as stream:
        CONFIG = yaml.load(stream)
else:
    CONFIG = config

if len(CONFIG["SAMPLES"]) == 0:
    sys.exit('No samples provided in your CONFIG.')
else:
    STEMS = CONFIG["SAMPLES"].split()

# ------------------------------ Directories--------------------------------- #

if CONFIG["INPUT-DIR"][-1] == "/":
    input_dir = CONFIG["INPUT-DIR"][:-1]
else:
    input_dir = CONFIG["INPUT-DIR"]

if CONFIG["TMP-DIR"][-1] == "/":
    tmp = CONFIG["TMP-DIR"][:-1]
else:
    tmp = CONFIG["TMP-DIR"]

if CONFIG["OUT-DIR"][-1] == "/":
    OUT = CONFIG["OUT-DIR"][:-1]
else:
    OUT = CONFIG["OUT-DIR"]

if CONFIG["SCRATCH"]:
    if CONFIG["SCRATCH"][-1] == "/":
        scratch = CONFIG["SCRATCH"][:-1]
    else:
        scratch = CONFIG["SCRATCH"]
else:
    scratch = tmp

LOGS = OUT+"/"+"LOGS"


# ------------------------------ MACHINE ------------------------------------ #

STAR_THREADS = CONFIG["MAPPING"]["STAR"]["threads"]
JULIA_THREADS = CONFIG["MACHINE"]["threads_julia"]
SAMBAMBA_THREADS = CONFIG['MACHINE']['threads_sambamba']
MEMORY_JAVA = CONFIG["MACHINE"]["memory_java"]


# ------------------------------ REFERENCE ---------------------------------- #

REFERENCE = CONFIG["REFERENCE"]
GTF = CONFIG["GTF"]
gff = CONFIG["GFF3"]
TRANSCRIPTS = CONFIG["TRANSCRIPTS"]


# ------------------------------ Loc ---------------------------------------- #
ALIGN_DIR = OUT + '/ALIGNMENT'
BENCHMARKS = CONFIG['BENCHMARKS']
MULTIQC_DIR = OUT + '/MULTIQC'
LOG_PATH = {
    'WORKFLOW_LOG_FILE': LOGS + '/SNAKEMAKE/workflow.log'
    }

LOG_CONFIG = 'LOG_CONFIG.yaml'
