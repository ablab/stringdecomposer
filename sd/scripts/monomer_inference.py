#!/usr/bin/env python

# requirements: clustalw2
import os
import sys
import csv
import argparse
import shutil
import edlib

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO

from threading import Thread

#Log class, use it, not print
class Log:
    text = ""

    def log(self, s):
        self.text += s + "\n"
        print(s)

    def warn(self, s):
        msg = "WARNING: " + s
        self.text += msg + "\n"
        sys.stdout.write(msg)
        sys.stdout.flush()

    def err(self, s):
        msg = "ERROR: " + s + "\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def print_log(self):
        print(self.text)

    def get_log(self):
        return self.text

log = Log()


def process_readline(line, is_python3=sys.version.startswith("3.")):
    if is_python3:
        return str(line, "utf-8").rstrip()
    return line.rstrip()


def sys_call(cmd):
    import shlex
    import subprocess

    if isinstance(cmd, list):
        cmd_list = cmd
    else:
        cmd_list = shlex.split(cmd)

    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    while not proc.poll():
        line = process_readline(proc.stdout.readline())
        if line:
            log.log(line)
        if proc.returncode is not None:
            break

    for line in proc.stdout.readlines():
        line = process_readline(line)
        if line:
            log.log(line)

    if proc.returncode:
        log.err("system call for: \"%s\" finished abnormally, OS return value: %d" % (cmd, proc.returncode))


def parse_args():
    parser = argparse.ArgumentParser(description='Monomer Inference Problem: complement monomers set')
    parser.add_argument('-seq', '--sequences', dest="sequences", help='fasta-file with long reads or genomic sequences')
    parser.add_argument('-mon', '--monomers', dest='monomers', help='fasta-file with monomers')
    parser.add_argument('-o', '--out-dir', dest="outdir", help='output directory (default=\'.\')', default=".", required=False)
    parser.add_argument("--continue", dest="restart", help="continue run from output dir", required=False, action='store_true')
    parser.add_argument('-t', '--threads', dest="threads", help='number of threads (default=1)', default=1, type=int)
    parser.add_argument('--len', help='the monomer length (default=171)', type=int, default=171, required=False)
    parser.add_argument('--lenDiv', '--max-length-divergence',
                        help='the maximum differ from length (default= 0.02*len)',
                        type=int, default=-1, required=False)
    parser.add_argument('--resDiv', '--max-resolved-divergence', help='max divergence in identity for resolve block (default=5%)',
                        type=float, default=5, required=False)
    parser.add_argument('--maxDiv','--max-divergence', help='max divergence in identity for monomeric-block (default=25%)',
                        type=float, default=40, required=False)
    return parser.parse_args()


class MonomericBlock(object):
    read_name = ""
    lft = 0
    rgh = 0
    seq = ""

    def __init__(self, read_name="", lft=0, rgh=0, seq=""):
        self.read_name = read_name
        self.lft = lft
        self.rgh = rgh
        self.seq = seq


def load_fasta(filename, tp = "list"):
    if tp == "map":
        records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
        for r in records:
            records[r] = records[r].upper()
    else:
        records = list(SeqIO.parse(filename, "fasta"))
        for i in range(len(records)):
            records[i] = records[i].upper()
    return records


def set_blocks_seq(sequences, blocks):
    log.log("Get block seqs")
    records = load_fasta(sequences, tp="map")
    for i in range(len(blocks)):
        #print("Block read name:", blocks[i].read_name)
        #print(records[blocks[i].read_name])
        #print(records[blocks[i].read_name][2:6])
        blocks[i].seq = records[blocks[i].read_name][blocks[i].lft:blocks[i].rgh + 1]
        #print("Set seq: " + blocks[i].seq)


def seq_identity(a, b):
    result = edlib.align(a, b, mode="NW", task="locations")
    if result["editDistance"] == -1:
        return 10**9
    return result["editDistance"] * 100 / max(len(a), len(b))


def clustering(blocks, args):
    from scipy.cluster.hierarchy import linkage

    y = []
    for i in range(len(blocks)):
        for j in range(i + 1, len(blocks)):
            y.append(seq_identity(str(blocks[i].seq.seq), str(blocks[j].seq.seq)))

    z = linkage(y, 'single')
    #[[cluster1, cluster2, dist, size]]
    #choose the biggest cluster with dist <= args.resDiv/2
    bst_cluster_id = 0
    cluster_size = 0
    for i in range(len(z)):
        if z[i][2] > args.resDiv/2:
            continue

        if z[i][3] > cluster_size:
            cluster_size = z[i][3]
            bst_cluster_id = i + len(blocks)

    #print(z)
    #generate list of blocks in bigest cluster
    max_cluster = []

    def add_blocks_to_cluster(cid):
        if (cid < len(blocks)):
            max_cluster.append(blocks[cid])
        else:
            add_blocks_to_cluster(int(z[cid - len(blocks)][0]))
            add_blocks_to_cluster(int(z[cid - len(blocks)][1]))

    add_blocks_to_cluster(int(bst_cluster_id))
    log.log("Max cluster is found! Cluster size: " + str(len(max_cluster)))
    return max_cluster


def save_seqs(max_cluster, cluster_seqs_path):
    with open(cluster_seqs_path, "w") as fa:
        for i in range(len(max_cluster)):
            name = "block" + str(i) + "_" + str(max_cluster[i].lft) + "_" + str(max_cluster[i].rgh)
            new_record = SeqRecord(max_cluster[i].seq.seq, id=name, description="")
            SeqIO.write(new_record, fa, "fasta")
    log.log("Seqs from biggest cluster saved to " + cluster_seqs_path)


def get_consensus_seq(cluster_seqs_path, seq_records):
    if (len(seq_records) == 1):
        return str(seq_records[0].seq.seq)

    from Bio.Align.Applications import ClustalwCommandline
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from Bio.Align import MultipleSeqAlignment

    cmd = ClustalwCommandline("clustalw2", infile=cluster_seqs_path)
    log.log("Run clustalw2: " + str(cmd))
    stdout, stderr = cmd()
    aln_file = '.'.join(cluster_seqs_path.split('.')[:-1]) + ".aln"
    log.log("Get Multiply alignment: " + aln_file)

    log.log("Start search for consensus monomer")
    align = AlignIO.read(aln_file, "clustal")
    #print(align.format("fasta"))

    summary_align = AlignInfo.SummaryInfo(align)
    consensus = summary_align.gap_consensus(threshold=0, ambiguous='N')
    consensus = str(consensus).replace('-', '')
    log.log("New consensus monomer: " + str(consensus))
    return consensus


def get_rc(seq):
    res_seq = ""
    for i in range(len(seq) - 1, -1, -1):
        if (seq[i] == 'A' or seq[i] == 'a'):
            res_seq += 'T'
        if (seq[i] == 'T' or seq[i] == 't'):
            res_seq += 'A'
        if (seq[i] == 'C' or seq[i] == 'c'):
            res_seq += 'G'
        if (seq[i] == 'G' or seq[i] == 'g'):
            res_seq += 'C'
    return res_seq


def get_dist_to_exists_monomers(monomers_list, new_monomer):
    mdist = 100
    for i in range(len(monomers_list)):
        mdist = min(mdist, seq_identity(str(monomers_list[i].seq), new_monomer))
        mdist = min(mdist, seq_identity(get_rc(str(monomers_list[i].seq)), new_monomer))
    return mdist


def init_first_run(args):
    if args.lenDiv == -1:
        args.lenDiv = int(args.len * 0.02)
    args.sequences = os.path.abspath(args.sequences)
    args.monomers = os.path.abspath(args.monomers)

    #switch working directory to output dir
    os.chdir(args.outdir)
    local_monmers_path = os.path.join(args.outdir, "monomers.fa")
    shutil.copyfile(args.monomers, local_monmers_path)
    args.monomers = local_monmers_path

    local_seq_path = os.path.join(args.outdir, 'sequence.fa')
    shutil.copyfile(args.sequences, local_seq_path)
    args.sequences = local_seq_path

    summary_path = os.path.join(args.outdir, "summary.csv")
    summary_fw = open(summary_path, "w")
    summary_writer = csv.writer(summary_fw)
    summary_writer.writerow(["IterationId", "Number of resolved blocks", "Number of unresolved blocks", "Number of non-monomeric blocks",
                              "Size of the largest cluster", "Number of deleted monomers at this iteration"])
    iter_id = 0
    return iter_id, summary_fw, summary_writer


def get_lst_iter(outdir):
    iter_id = 0
    while (os.path.isdir(os.path.join(outdir, "iter_" + str(iter_id)))):
        iter_id += 1
    return iter_id - 1


def init_continue(args):
    #switch working directory to output dir
    os.chdir(args.outdir)

    iter_id = get_lst_iter(args.outdir)

    args.monomers = os.path.join(args.outdir, "monomers.fa")
    args.sequences = os.path.join(args.outdir, 'sequence.fa')

    summary_path = os.path.join(args.outdir, "summary.csv")
    summary_fw = open(summary_path, "a")
    summary_writer = csv.writer(summary_fw)

    return iter_id, summary_fw, summary_writer


def init_monomers(monomers_path):
    #save info about monomers
    monomers_list = []
    for record in SeqIO.parse(monomers_path, "fasta"):
        monomers_list.append(record)
    return monomers_list


need_reverse_monomers = -1


def reverse_monomer(args):
    global need_reverse_monomers
    if need_reverse_monomers != -1:
        return need_reverse_monomers
    else:
        res_tsv = os.path.join(args.outdir, "iter_0", "final_decomposition.tsv")
        cnt_need_rev = 0
        cnt_dont_need_rev = 0
        with open(res_tsv, "r") as f:
            csv_reader = csv.reader(f, delimiter='\t')
            for row in csv_reader:
                if row[2] == "start":
                    continue
                identity = float(row[4])

                if identity >= 100 - args.maxDiv:
                    if row[1][-1] == "'":
                        cnt_need_rev += 1
                    else:
                        cnt_dont_need_rev += 1
        if cnt_need_rev > cnt_dont_need_rev:
            need_reverse_monomers = 1
        else:
            need_reverse_monomers = 0
        return need_reverse_monomers


def rc(seq):
    res_seq = ""
    for i in range(len(seq)):
        pos = len(seq) - i - 1
        if seq[pos] == 'A':
            res_seq += "T"
        elif seq[pos] == 'T':
            res_seq += "A"
        elif seq[pos] == 'C':
            res_seq += "G"
        elif seq[pos] == "G":
            res_seq += "C"
        else:
            res_seq += seq[pos]
    return res_seq


def need_update(args, monomer_record, monomer_resolved_block, iter_outdir):
    iter_id = int(iter_outdir.split('_')[-1])
    if iter_id == 0:
        return True
    else:
        prev_dir = '_'.join(iter_outdir.split('_')[:-1]) + '_' + str(iter_id - 1)
        res_tsv = os.path.join(prev_dir, "final_decomposition.tsv")
        prevresolved_blocks = []

        with open(res_tsv, "r") as f:
            csv_reader = csv.reader(f, delimiter='\t')
            for row in csv_reader:
                if row[2] == "start":
                    continue
                identity = float(row[4])
                if identity >= 100 - args.resDiv:
                    if row[1][-1] == "'":
                        row[1] = row[1][:-1]
                    if row[1] == monomer_record.id:
                        prevresolved_blocks.append(MonomericBlock(row[0], int(row[2]), int(row[3])))
        if len(prevresolved_blocks) != len(monomer_resolved_block):
            return True

        for i in range(len(prevresolved_blocks)):
            if (prevresolved_blocks[i].read_name != monomer_resolved_block[i].read_name) or\
                    (prevresolved_blocks[i].lft != monomer_resolved_block[i].lft) or\
                    (prevresolved_blocks[i].rgh != monomer_resolved_block[i].rgh):
                return True
    return False


def update_monomer(args, monomer_record, monomer_resolved_block, iter_outdir):
    if need_update(args, monomer_record, monomer_resolved_block, iter_outdir):
        cluster_seqs_path = os.path.join(iter_outdir, monomer_record.id.split('/')[0] + "_seqs.fa")
        save_seqs(monomer_resolved_block, cluster_seqs_path)
        new_monomer = get_consensus_seq(cluster_seqs_path, monomer_resolved_block)
        if reverse_monomer(args):
            new_monomer = rc(new_monomer)
        monomer_record.seq = Seq(new_monomer)

    return monomer_record


def update_monomers_range(lft, rgh, args, monomer_resolved, monomers_list, iter_outdir):
    for i in range(lft, rgh):
        log.log("==== Update monomer " + str(monomers_list[i].id))
        set_blocks_seq(args.sequences, monomer_resolved[monomers_list[i].id])
        monomers_list[i] = update_monomer(args, monomers_list[i], monomer_resolved[monomers_list[i].id], iter_outdir)


def final_iteration(args, sd_script_path, monomers_list):
    log.log("====== Start final iteration======")
    iter_outdir = os.path.join(args.outdir, "final")
    if not os.path.exists(iter_outdir):
        os.makedirs(iter_outdir)
    os.chdir(iter_outdir)

    local_monmers_path = os.path.join(iter_outdir, "monomers.fa")
    shutil.copyfile(args.monomers, local_monmers_path)

    # run string decomposer
    sys_call(["python3", sd_script_path, args.sequences, local_monmers_path, "-t", str(args.threads)])
    log.log("String decomposer is complete. Results save in: " + iter_outdir)

    G = {}
    gw = {}
    for i in range(len(monomers_list)):
        G[monomers_list[i].id] = []
        for j in range(len(monomers_list)):
            gw[(monomers_list[i].id, monomers_list[j].id)] = 0

    # parse output csv file
    res_tsv = os.path.join(iter_outdir, "final_decomposition.tsv")

    with open(res_tsv, "r") as f:
        csv_reader = csv.reader(f, delimiter='\t')
        prev_row = []
        for row in csv_reader:
            if row[2] == "start":
                continue
            if prev_row != []:
                pident = float(prev_row[4])
                pmon = prev_row[1]
                if pmon[-1] == "'":
                    pmon = pmon[:-1]
                identity = float(row[4])
                mon = row[1]
                if mon[-1] == "'":
                    mon = mon[:-1]

                if identity >= 100 - args.resDiv and pident >= 100 - args.resDiv:
                    gw[(pmon, mon)] += 1
            prev_row = row

    dotst = "digraph MonomerGraph {\n"

    for i in range(len(monomers_list)):
        dotst += monomers_list[i].id + ";\n"
        G[monomers_list[i].id] = []
        for j in range(len(monomers_list)):
            if gw[(monomers_list[i].id, monomers_list[j].id)] > 0:
                G[monomers_list[i].id].append((monomers_list[j].id, gw[(monomers_list[i].id, monomers_list[j].id)]))

    simple_gr = os.path.join(iter_outdir, "monomer_gr")
    with open(simple_gr, "w") as fw:
        for i in range(len(monomers_list)):
            curm = monomers_list[i].id
            for edg in G[curm]:
                fw.write(curm + " " + edg[0] + " "  + str(edg[1]) + "\n")
                dotst += curm + " -> " + edg[0] + " [label=\"" + str(edg[1]) + "\"];\n"

    dotst += "}\n"

    dotpath = os.path.join(iter_outdir, "monomer_graph.dot")
    with open(dotpath, "w") as fw:
        fw.write(dotst)


def main():
    log.log("Start Monomer Inference")
    args = parse_args()
    # get path to current script
    current_script_path = os.path.abspath(os.path.dirname(__file__))

    # path to the run_decomposer script
    sd_script_path = os.path.join(current_script_path, "..", "run_decomposer.py")
    log.log("Path to run_decomposer: " + sd_script_path)

    # create output_dir if not exists
    args.outdir = os.path.abspath(args.outdir)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if args.restart == False:
        iter_id, summary_fw, summary_writer = init_first_run(args)
    else:
        iter_id, summary_fw, summary_writer = init_continue(args)

    monomers_list = init_monomers(args.monomers)

    monomer_set_complete = False
    while (not monomer_set_complete):
        log.log("====== Start iteration " + str(iter_id) + "======")
        #create for current iteration
        iter_outdir = os.path.join(args.outdir, "iter_" + str(iter_id))
        if not os.path.exists(iter_outdir):
            os.makedirs(iter_outdir)
        os.chdir(iter_outdir)
        local_monmers_path = os.path.join(iter_outdir, "monomers.fa")
        shutil.copyfile(args.monomers, local_monmers_path)

        # run string decomposer
        sys_call(["python3", sd_script_path, args.sequences, local_monmers_path, "-t", str(args.threads)])
        log.log("String decomposer is complete. Results save in: " + iter_outdir)

        # parse output csv file
        res_tsv = os.path.join(iter_outdir, "final_decomposition.tsv")
        unresolved_blocks = []
        non_monomeric_cnt = 0
        resolved_cnt = 0

        monomer_resolved = {}
        for i in range(len(monomers_list)):
            monomer_resolved[monomers_list[i].id] = []

        with open(res_tsv, "r") as f:
            csv_reader = csv.reader(f, delimiter='\t')
            for row in csv_reader:
                if row[2] == "start":
                    continue
                identity = float(row[4])
                if identity >= 100 - args.resDiv:
                    resolved_cnt += 1
                    if row[1][-1] == "'":
                        row[1] = row[1][:-1]
                    monomer_resolved[row[1]].append(MonomericBlock(row[0], int(row[2]), int(row[3])))

                if identity <= 100 - args.maxDiv:
                    non_monomeric_cnt += 1

                if (identity > 100 - args.maxDiv) and (identity < 100 - args.resDiv):
                    unresolved_blocks.append(MonomericBlock(row[0], int(row[2]), int(row[3])))

        deleted_cnt = 0
        i = 0
        while (i < len(monomers_list)):
            if len(monomer_resolved[monomers_list[i].id]) == 0:
                deleted_cnt += 1
                if (i + 1 < len(monomers_list)):
                    monomers_list = monomers_list[:i] + monomers_list[i + 1:]
                else:
                    monomers_list = monomers_list[:i]
            else:
                i += 1

        stp = (1 + (len(monomers_list) - 1)//args.threads)
        lft = 0
        threads = []
        while lft < len(monomers_list):
            threads.append(Thread(target=update_monomers_range, args=(lft, min(lft + stp, len(monomers_list)),
                                                                      args, monomer_resolved, monomers_list, iter_outdir)))
            lft += stp

        for i in range(len(threads)):
            threads[i].start()
        for i in range(len(threads)):
            threads[i].join()

        log.log("Number of unresolved monomer block: " + str(len(unresolved_blocks)))
        log.log("Number of resolved blocks: " + str(resolved_cnt))
        log.log("Number of non-monomeric blocks: " + str(non_monomeric_cnt))
        log.log("Number of deleted monomers: " + str(deleted_cnt))

        if len(unresolved_blocks) > 0:
            set_blocks_seq(args.sequences, unresolved_blocks)
            max_cluster = clustering(unresolved_blocks, args)
            if (len(max_cluster) == 1):
                monomer_set_complete = True
        else:
            monomer_set_complete = True

        if not monomer_set_complete:
            cluster_seqs_path = os.path.join(iter_outdir, "cluster_seq.fa")
            save_seqs(max_cluster, cluster_seqs_path)

            new_monomer = get_consensus_seq(cluster_seqs_path, max_cluster)
            if reverse_monomer(args):
                new_monomer = rc(new_monomer)

            dist_to_monomers = get_dist_to_exists_monomers(monomers_list, new_monomer)
            log.log("Min Distance to exsisting monomers: " + str(dist_to_monomers))

            new_monomer_record = SeqRecord(Seq(new_monomer), id="mn_" + str(iter_id), description="")
            monomers_list.append(new_monomer_record)

        with open(args.monomers, "w") as fa:
            for record in monomers_list:
                SeqIO.write(record, fa, "fasta")

        summary_writer.writerow([str(iter_id), str(resolved_cnt), str(len(unresolved_blocks)), str(non_monomeric_cnt), str(len(max_cluster)), str(deleted_cnt)])
        iter_id += 1

    final_iteration(args, sd_script_path, monomers_list)
    summary_fw.close()

if __name__ == "__main__":
    main()