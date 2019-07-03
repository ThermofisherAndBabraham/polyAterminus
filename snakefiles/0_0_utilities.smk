import logging
import logging.config
import os.path
import re
import sys
import yaml
import os


def init_log():
    try:
        with open(LOG_CONFIG, 'rt') as f:
            log_config = yaml.safe_load(f.read())
        for handler in log_config['handlers']:
            fl = log_config['handlers'][handler]['filename']
            log_config['handlers'][handler]['filename'] = LOG_PATH[fl]
            log_dir = LOG_PATH[fl].split('/')[:-1]
            log_dir = '/'.join(log_dir)
            if not os.path.isdir(log_dir):
                os.makedirs(log_dir)

        logging.config.dictConfig(log_config)

    except Exception as e:
        print('LOGS creation failed! ' + str(e))
        sys.exit()


def install_slimfastq():

    try:
        subprocess.call(['dependencies/slimfastq/slimfastq', '-h'], stdout = subprocess.DEVNULL)
    except:
        if not os.path.exists('dependencies'):
            os.mkdir('dependencies')

        wlogger.info(f'slimfastq was not found in dependencies. Pulling.]')
        p = subprocess.Popen('cd dependencies; ' +
                              'git clone https://github.com/Infinidat/slimfastq.git; ' +
                              'cd slimfastq ' +
                              'git checkout v1.032 && make; ' +
                              'cd ../../snakefiles',
                              stdout=subprocess.PIPE,
                              shell=True
                              )
        p_status = p.wait()
        (p_out, p_err) = p.communicate()
        wlogger.info(f'slimfastq was not found in dependencies. Pulling. DONE]')
        if p_status == 0:
            logger.info(f'slimfastq installed to dependencies: {p_out}')
        else:
            logger.error(f'Installing slimfastq exited with 0: {p_err}')

    return True


def shstring(l, sfq):
    slimfq = 'dependencies/slimfastq/slimfastq '
    out = ''
    for i in l:
        if i.endswith('.sfq') and sfq:
            out = out + slimfq + i + ' | gzip -f; \n'
        else:
            out = out + "cat " + i + "; \n"
    return out


def get_sample_files_named_pipe_script(sample, debug=True):

    input_dir = CONFIG['INPUT-DIR']
    r1_files = []
    r1_stems = []
    r2_files = []
    r2_stems = []
    r1_testlist = []
    r2_testlist = []
    sfq = False


    for root, dirnames, filenames in os.walk(input_dir):
        for f in filenames:
            if ((f.endswith(".fastq.gz") or f.endswith(".fq.gz") or f.endswith(".sfq"))
                    and f.find(sample) > -1):
                if f.find(".sfq") > -1:
                    sfq = True

                lane_splits = re.split(CONFIG["LANE-REGEX"], f)
                wlogger.info(f'File name split by LANE-REGEX to: [{lane_splits}]')

                if (f.find("_R1_") > -1 or f.endswith("_1.fq.gz")
                        or f.endswith("_1.fastq.gz") or f.endswith("_1.sfq")):

                    if debug:
                        wlogger.info(f'Found R1 file: [{root}, {f}]')
                    r1_files.append(root+os.sep+f)

                    if lane_splits[0].find(sample) > -1:
                        stem = lane_splits[0]
                    else:
                        wlogger.error(f'SAMPLE not found in split filename: [{sample}, {lane_splits}]')
                        wlogger.error(f'Check LANE-REGEX and SAMPLES in CONFIG!')

                        raise Exception(
                            "Substring " + sample + " is not found in split filename string "
                            "by LANE-REGEX: " + CONFIG["LANE-REGEX"] + "\n"
                            "Filename: " + f + "\n"
                            "SAMPLE: " + sample
                            )

                    r1_testlist.append(root + os.sep + stem)
                    r1_stems.append(stem)

                elif (f.find("_R2_") > -1 or f.endswith("_2.fq.gz")
                        or f.endswith("_2.fastq.gz") or f.endswith("_2.sfq")):

                    if debug:
                        wlogger.info(f'Found R2 file: [{root}, {f}]')
                    r2_files.append(root + os.sep + f)

                    if lane_splits[0].find(sample) > -1:
                        stem = lane_splits[0]
                    else:
                        wlogger.error(f'SAMPLE not found in split filename: [{sample}, {lane_splits}]')
                        wlogger.error(f'Check LANE-REGEX and SAMPLES in CONFIG!')

                        raise Exception(
                            "Substring " + sample + " is not found in split filename string "
                            "by LANE-REGEX: " + CONFIG["LANE-REGEX"] + "\n"
                            "Filename: " + f + "\n"
                            "SAMPLE: " + sample
                            )

                    r2_testlist.append(root + os.sep + stem)
                    r2_stems.append(stem)

    all_stems = []
    all_stems.extend(r1_stems)
    all_stems.extend(r2_stems)

    if len(list(set(all_stems))) > 1:
        wlogger.warning(f'Split filename by LANE-REGEX is matching to more than one SAMPLES: [{sample}, {all_stems}]')

    if str(r1_testlist) != str(r2_testlist):
        wlogger.error(f'R1 and R2 found lists of files does not mathch, pairs are invalid!')
        r1_testlist.sort()
        wlogger.error(f'R1 file begins: [{r1_testlist}]')

        r2_testlist.sort()
        wlogger.error(f'R2 file begins: [{r2_testlist}]')

        raise Exception(
            "Key " + sample + " is not properly differentiating, "
            "pairs are invalid \n"
            "Check lists of R1: " + str(r1_testlist) + "\n"
            "and R2: " + str(r2_testlist)
            )

    if len(r1_files) == 0 or len(r2_files) == 0:
        if not os.path.exists(input_dir):
            wlogger.error(f'INPUT-DIR does not exists: [{input_dir}]')
            print ("{0} directory do not exit!".format(input_dir))
        else:
            wlogger.error(f'No suitable files found: [INPUT-DIR:{input_dir}, SAMPLE:{sample}]')
            wlogger.error(os.system("ls {0}".format(input_dir)))

        raise Exception(
            "No matching R1 files for stem {0} in dir: {1}".format(sample, input_dir))

    r1_files.sort()
    r2_files.sort()
    r1ft = os.path.join(tmp+'/{0}R1_input_list.sh'.format(sample))
    r2ft = os.path.join(tmp+'/{0}R2_input_list.sh'.format(sample))

    f1 = open(r1ft, 'w')

    wlogger.info(f'Writting {sample} R1 list of found files for `cat`: [{r1_files}]')
    f1.write(shstring(r1_files, sfq))
    f1.close()
    os.system('chmod a+x %s' % (r1ft))

    f2 = open(r2ft, 'w')
    out = shstring(r2_files, sfq)
    wlogger.info(f'Writting {sample} R2 list of found files for `cat`: [{r2_files}]')
    f2.write(out)
    f2.close()
    os.system('chmod a+x %s' % (r2ft))


def PolyAAnalysis_trimm_input(wildcards):

    if CONFIG["SUBSAMPLING"]["run"]:
        input = [
            tmp + "/" + wildcards.stem + "_R1_001subs.fastq.gz",
            tmp + "/" + wildcards.stem + "_R2_001subs.fastq.gz"
        ]

    else:
        input = [
            tmp + "/" + wildcards.stem + "_R1_001Trimmed.fastq.gz",
            tmp + "/" + wildcards.stem + "_R2_001Trimmed.fastq.gz"
        ]

    return input


def get_annotate_ts_strandedness(wildcards):

    for i in CONFIG["ANNOTATE-TS"]["strandedness"].keys():
        if i == wildcards.stem:
            wlogger.info(f'Set strandedness: {wildcards.stem}: {CONFIG["ANNOTATE-TS"]["strandedness"][i]}')
            return CONFIG["ANNOTATE-TS"]["strandedness"][i]
