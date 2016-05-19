import argparse
import glob
import os
import subprocess
import pandas as pd

# Generating parameters for MiGec assembly
def assemble_param(minimal_overseq):
    global output_dir
    samples_overseq = {}
    with open("migec/histogram/estimates.txt") as threshold:
        for line in threshold:
            if minimal_overseq is None:
                samples_overseq[line.split()[0]] = line.split()[4]
                output_dir = "assemble"
            else:
                samples_overseq[line.split()[0]] = minimal_overseq
                output_dir = "assemble_t" + str(minimal_overseq)
    return samples_overseq, output_dir

# Creating metadata file for VDJtools
def metadata_creator():
    label_list = []
    for file in glob.glob("mixcr/*.vdjca"):
        file_label_list = []
        file_id = os.path.splitext(os.path.basename(file))[0]
        file_label_list.append(file_id + ".txt")
        file_label_list.append(file_id)
        file_label_list.extend(file_id.split("_"))
        label_list.append(file_label_list)
    maxLen = max(len(l) for l in label_list)
    metadata = pd.DataFrame(label_list)
    col_names = ["#file name", "sample_id", "label_1"]
    for i in range(3, maxLen):
        col_names.append("label_" + str(i - 1))
    metadata.columns = col_names
    metadata.to_csv("mixcr/metadata.txt", sep='\t', index=False, na_rep="NA")


def migec_checkout(barcodesFile):
    FNULL = open(os.devnull, 'w')
    demultiplexing = subprocess.Popen(
        ['migec', 'CheckoutBatch', '-cute', '--skip-undef', barcodesFile, 'migec/checkout/'],
        stdout=FNULL, stderr=FNULL)
    demultiplexing.wait()


def migec_histogram():
    FNULL = open(os.devnull, 'w')
    hist = subprocess.Popen(['migec', 'Histogram', 'migec/checkout/', 'migec/histogram/'], stdout=FNULL, stderr=FNULL)
    hist.wait()


def migec_assemble(file_R1, file_R2, overseq, output_dir):
    FNULL = open(os.devnull, 'w')
    assemble = subprocess.Popen(['migec', 'Assemble', '-m', overseq, '--filter-collisions', file_R1, file_R2,
                                 "migec/" + output_dir + "/"], stdout=FNULL, stderr=FNULL)
    assemble.wait()


def mixcr(species, file_R1, file_R2):
    print("Starting MiXCR alignment for " + os.path.splitext(os.path.basename(file_R1))[0].split("_R1")[0])
    FNULL = open(os.devnull, 'w')
    mixcr_alignment = subprocess.Popen(['mixcr', 'align', '-r mixcr/alignmentReport.txt', '-f', '-s', species,
                                        file_R1, file_R2,
                                        'mixcr/' + os.path.splitext(os.path.basename(file_R1))[0].split("_R1")[
                                            0] + '.vdjca'],
                                       stdout=FNULL, stderr=FNULL)
    mixcr_alignment.wait()
    print("Starting MiXCR assemble for " + os.path.splitext(os.path.basename(file_R1))[0].split("_R1")[0])

    mixcr_assemble = subprocess.Popen(['mixcr', 'assemble', '-r mixcr/assembleReport.txt', '-f', 'mixcr/' +
                                       os.path.splitext(os.path.basename(file_R1))[0].split("_R1")[0] + '.vdjca',
                                       'mixcr/' + os.path.splitext(os.path.basename(file_R1))[0].split("_R1")[
                                           0] + '.clns'],
                                      stdout=FNULL, stderr=FNULL)
    mixcr_assemble.wait()
    print("Exporting clones for " + os.path.splitext(os.path.basename(file_R1))[0].split("_R1")[0])
    mixcr_export = subprocess.Popen(
        ['mixcr', 'exportClones', 'mixcr/' + os.path.splitext(os.path.basename(file_R1))[0].split("_R1")[0] + '.clns',
         'mixcr/' + os.path.splitext(os.path.basename(file_R1))[0].split("_R1")[0] + '.txt'],
        stdout=FNULL, stderr=FNULL)
    mixcr_export.wait()

# Converting mixcr output for VDJTools, calc basic stats
def vdjtools():
    print("Converting files to vdgtools format")
    FNULL = open(os.devnull, 'w')
    vdjtools_convert = subprocess.Popen(['vdjtools', 'Convert', '-S', 'MiXCR', '-m', 'mixcr/metadata.txt',
                                        'vdjtools/'], stdout=FNULL, stderr=FNULL)
    vdjtools_convert.wait()
    print("Calculating basic statistics")
    vdjtools_basicstats = subprocess.Popen(['vdjtools', 'CalcBasicStats', '-m', 'vdjtools/metadata.txt', 'vdjtools/'],
                                           stdout=FNULL, stderr=FNULL)
    vdjtools_basicstats.wait()


def pipeline(barcodesFile, species, minimal_overseq):
    print("\033[1;36;40mMiBuddy will take care of your data\033[0m")
    print("Starting demultiplexing")
    migec_checkout(barcodesFile)
    print("Demultiplexing is complete")
    print("Collecting MIG statistics")
    migec_histogram()
    print("MIG statistics has been calculated")
    samples_overseq = assemble_param(minimal_overseq)[0]
    assemble_path = assemble_param(minimal_overseq)[1]
    for file in glob.glob("migec/checkout/*_R1.fastq.gz"):
        filename = os.path.splitext(os.path.basename(file))[0]
        if filename.split("_R1")[0] in samples_overseq.keys():
            print("Assembling MIGs for " + filename.split("_R1")[0] + ". Minimal number of reeds per MIG: " +
                  samples_overseq[filename.split("_R1")[0]])
            file_1_path = "migec/checkout/" + filename.split("_R1")[0] + "_R1" + ".fastq.gz"
            file_2_path = "migec/checkout/" + filename.split("_R1")[0] + "_R2" + ".fastq.gz"
            overseq = samples_overseq[filename.split("_R1")[0]]
            migec_assemble(file_1_path, file_2_path, str(overseq), assemble_path)
            mixcr(species, glob.glob("migec/" + assemble_path + "/" + filename.split("_R1")[0] + "_R1*.fastq")[0],
                  glob.glob("migec/" + assemble_path + "/" + filename.split("_R1")[0] + "_R2*.fastq")[0])
    print("Creating metadata file")
    metadata_creator()
    vdjtools()


def main(args):
    dirs = ["migec", "mixcr", "vdjtools"]
    for item in dirs:
        if not os.path.exists(item):
            os.makedirs(item)
    pipeline(args.file_with_barcodes, args.s, args.overseq)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("file_with_barcodes", help="Specify barcodes file", required=True)
    parser.add_argument("-s", help="Specify species: mmu for Mus musculus, hsa - Homo sapiens", required=True)
    parser.add_argument("--overseq", "-minimal_overseq", type=int, default=None,
                        help="Force minimal overseq value for all samples")
    args = parser.parse_args()
    main(args)
