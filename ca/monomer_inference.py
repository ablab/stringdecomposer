#!/usr/bin/env python

# requirements: clustalw2
import os
import sys
import csv
import argparse
import shutil
import edlib
import multiprocessing

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

FreqCeiling = 40
TotalMonomerBlocks = 0

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
    parser.add_argument('-o', '--out-dir', dest="outdir", help="output directory [default='.']", default=".",
                        required=False)
    parser.add_argument("--continue", dest="restart", help="continue run from output dir", required=False,
                        action='store_true')
    parser.add_argument('-t', '--threads', dest="threads", help='number of threads [default=1]', default=1, type=int)
    parser.add_argument('--resDiv', '--max-resolved-divergence',
                        help='max divergence in identity for resolve block [default=5]',
                        type=float, default=5, required=False)
    parser.add_argument('--maxDiv', '--max-divergence',
                        help='max divergence in identity for monomeric-block [default=25]',
                        type=float, default=40, required=False)
    parser.add_argument("--min-cluster-size", dest="min_clst", type=int, default=1, help=
                                                                    "When maximum size of cluster will be less "\
                                                                    "MIN_CLST the monomer inferense will finished "\
                                                                    "[default=2]")
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

    def __hash__(self):
        return hash((self.lft, self.rgh, self.read_name))

    def __eq__(self, other):
        return self.lft == other.lft and self.rgh == other.rgh and self.read_name == other.read_name


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


def get_distance(a, b):
    result = edlib.align(a, b, mode="NW", task="locations")
    return result["editDistance"]


def seq_identity(a, b):
    result = edlib.align(a, b, mode="NW", task="locations")
    if result["editDistance"] == -1:
        return 10**9
    return result["editDistance"] * 100 / max(len(a), len(b))


def get_clusters_list_id(z, args, n, separation=15):
    def inner_get_cluster_list(vid, mx_not_add_cluster=0):
        if vid < n:
            return [vid], mx_not_add_cluster

        zid = vid - n
        if z[zid][2]*2 <= args.resDiv:
            return [vid], mx_not_add_cluster
        elif z[zid][2] >= separation:
            vid1 = int(z[zid][0])
            vid2 = int(z[zid][1])
            clst1, mx1 = inner_get_cluster_list(vid1, mx_not_add_cluster)
            clst2, mx2 = inner_get_cluster_list(vid2, mx_not_add_cluster)
            return clst1 + clst2, max(mx1, mx2, mx_not_add_cluster)
        else:
            vid1 = int(z[zid][0])
            vid2 = int(z[zid][1])
            clst1, clst1_mx = inner_get_cluster_list(vid1, mx_not_add_cluster)
            clst2, clst2_mx = inner_get_cluster_list(vid2, mx_not_add_cluster)
            mx1 = (0, 0)
            mx2 = (0, 0)
            for cid in clst1:
                if cid < n:
                    mx1 = max((1, -cid), mx1)
                else:
                    mx1 = max((z[cid - n][3], -cid), mx1)

            for cid in clst2:
                if cid < n:
                    mx2 = max((1, -cid), mx2)
                else:
                    mx2 = max((z[cid - n][3], -cid), mx2)
            if mx1 > mx2:
                mx_not_add_cluster = max(mx_not_add_cluster, mx2[0], clst1_mx, clst2_mx)
                return clst1, mx_not_add_cluster
            else:
                mx_not_add_cluster = max(mx_not_add_cluster, mx1[0], clst1_mx, clst2_mx)
                return clst2, mx_not_add_cluster
    clstid, mx_nt_add = inner_get_cluster_list(n + len(z) - 1)


    res_clusters = []
    for cl in clstid:
        sz = 1
        if cl >= n:
            sz = int(z[cl - n][3])
        if sz >= mx_nt_add:
            res_clusters.append(cl)
    res_clusters.sort(key=lambda x: (-1,x) if x < n else (int(-z[x - n][3]), x))
    return res_clusters


def calc_dists_between_blocks(blocks, args):
    log.log("Calculate distance between monomers blocks")

    id_get_pairs = []
    for i in range(len(blocks)):
        for j in range(i + 1, len(blocks)):
            id_get_pairs.append((i, j))

    y = multiprocessing.Array('f', len(id_get_pairs), lock=False)

    def calc_dist_for_range(lft, rgh, y):
        for i in range(lft, rgh):
            if lft == 0 and i%100000 == 0:
                print(i, lft, rgh)
            y[i] = seq_identity(str(blocks[id_get_pairs[i][0]].seq.seq), str(blocks[id_get_pairs[i][1]].seq.seq))

    stp = (1 + (len(id_get_pairs) - 1) // args.threads)
    lft = 0

    threads = []
    while lft < len(id_get_pairs):
        threads.append(multiprocessing.Process(target=calc_dist_for_range,
                                               args=(lft, min(lft + stp, len(id_get_pairs)), y)))
        lft += stp

    for i in range(len(threads)):
        threads[i].start()
    for i in range(len(threads)):
        threads[i].join()
    y = y[:]
    return y


def init_small_ids(v, init_dist, small_ids, origin_ids, cluserts_id, blocks, separation, args):
    y = multiprocessing.Array('f', len(blocks), lock=False)

    def calc_dist_for_range(lft, rgh, y):
        for i in range(lft, rgh):
            if cluserts_id[i] == -1:
                y[i] = seq_identity(str(blocks[v].seq.seq), str(blocks[i].seq.seq))

    stp = (1 + (len(blocks) - 1) // args.threads)
    lft = 0
    threads = []
    while lft < len(blocks):
        rgt = min(lft + stp, len(blocks))
        threads.append(multiprocessing.Process(target=calc_dist_for_range,
                                               args=(lft, rgt, y)))
        lft = rgt

    for i in range(len(threads)):
        threads[i].start()
    for i in range(len(threads)):
        threads[i].join()

    cl_size = 0
    for j in range(len(blocks)):
        if cluserts_id[j] != -1:
            continue
        cur_dist = y[j]
        init_dist[j] = cur_dist
        if cluserts_id[j] == -1 and cur_dist < separation:
            origin_ids.append(j)
            small_ids[j] = len(origin_ids) - 1
        if cur_dist*2 <= args.resDiv:
            cl_size += 1
    return cl_size


def find_all_dists(small_ids, origin_ids, blocks, args):
    y = multiprocessing.Array('f', len(origin_ids) * len(origin_ids), lock=False)

    def calc_dist_for_range(lft, rgh, y):
        for i in range(lft, rgh):
            xx = origin_ids[i//len(origin_ids)]
            yy = origin_ids[i%len(origin_ids)]
            y[i] = seq_identity(str(blocks[xx].seq.seq), str(blocks[yy].seq.seq))

    stp = (1 + (len(origin_ids) * len(origin_ids) - 1) // args.threads)
    lft = 0

    threads = []
    while lft < len(origin_ids) * len(origin_ids):
        threads.append(multiprocessing.Process(target=calc_dist_for_range,
                                               args=(lft, min(lft + stp, len(origin_ids) * len(origin_ids)), y)))
        lft += stp

    for i in range(len(threads)):
        threads[i].start()
    for i in range(len(threads)):
        threads[i].join()
    return y[:]


cnt_resolved = 0

isolated_clusters = set()

def init_one_cluster(i, cluserts_id, blocks, args, cid):
    global isolated_clusters
    global cnt_resolved
    separation = 20
    origin_ids = []
    small_ids = [-1]*len(blocks)
    init_dist = [-1]*len(blocks)

    if blocks[i] in isolated_clusters:
        cluserts_id[i] = cid
        return

    cl_size = init_small_ids(i, init_dist, small_ids, origin_ids, cluserts_id, blocks, separation, args)
    print("CLuster size: ", cl_size, " len(origin_ids)=", str(len(origin_ids)))

    #log.log("Cnt blocks in area: " + str(len(origin_ids)))
    #print(all_dists)
    #log.log("All dists are calculated")

    cluserts_id[i] = cid
    queue = [i]
    bg = 0

    while (bg < len(queue)):
        cnt_resolved += 1
        if cnt_resolved % 100 == 0:
            print(str(cnt_resolved) + "/" + str(len(blocks)))
            print("CLose monomers: ", len(origin_ids))
            print("Cluster id: ", cid)
        v = queue[bg]
        bg += 1

        for j in range(0, len(origin_ids)):
            if cluserts_id[origin_ids[j]] != -1:
                continue
            if (init_dist[origin_ids[j]] - init_dist[v]) * 2 > args.resDiv:
                continue
            cur_dist = seq_identity(str(blocks[v].seq.seq), str(blocks[origin_ids[j]].seq.seq))
            if cur_dist*2 <= args.resDiv:
                cluserts_id[origin_ids[j]] = cid
                queue.append(origin_ids[j])

    if len(queue) == 1:
        isolated_clusters.add(blocks[i])
    log.log("Cluster id: " + str(cid) + "; area size: " + str(len(origin_ids)) + "; cluster size: " + str(len(queue)))


def get_clusters(blocks, args):
    global cnt_resolved
    cnt_resolved = 0
    cluserts_id = [-1]*len(blocks)
    cid = 0
    for i in range(len(blocks)):
        if cluserts_id[i] == -1:
            init_one_cluster(i, cluserts_id, blocks, args, cid)
            cid += 1

    clusters = [[] for i in range(cid)]
    for i in range(len(blocks)):
        clusters[cluserts_id[i]].append(blocks[i])

    clusters.sort(key=lambda x: -len(x))
    return clusters


def get_separation(clst1, clst2, args):
    y = multiprocessing.Value('f', 100, lock=True)

    def calc_dist_for_range(lft, rgh, y):
        mn_dist = 100
        for i in range(lft, rgh):
            for j in range(len(clst2)):
                mn_dist = min(mn_dist, seq_identity(str(clst1[i].seq.seq), str(clst2[j].seq.seq)))
        y.value = min(y.value, mn_dist)

    stp = (1 + (len(clst1) - 1) // args.threads)
    lft = 0

    threads = []
    while lft < len(clst1):
        threads.append(multiprocessing.Process(target=calc_dist_for_range,
                                               args=(lft, min(lft + stp, len(clst1)), y)))
        lft += stp

    for i in range(len(threads)):
        threads[i].start()
    for i in range(len(threads)):
        threads[i].join()
    return y.value


def clustering(blocks, args):
    log.log("===CLUSTERING===")

    clusters = get_clusters(blocks, args)
    max_clst = [clusters[0]]
    for i in range(1, len(clusters)):
        mn_sep = 100
        for j in range(0, i):
            mn_sep = min(mn_sep, get_separation(clusters[i], clusters[j], args))
        if mn_sep < 11:
            break
        else:
            max_clst.append(clusters[i])

    log.log("Max cluster is found! Cluster size: " + str(len(max_clst[0])))
    return max_clst


def save_seqs(max_cluster, cluster_seqs_path):
    with open(cluster_seqs_path, "w") as fa:
        for i in range(len(max_cluster)):
            name = "block" + str(i) + "_" + str(max_cluster[i].lft) + "_" + str(max_cluster[i].rgh)
            new_record = SeqRecord(max_cluster[i].seq.seq, id=name, description="")
            SeqIO.write(new_record, fa, "fasta")
    log.log("Seqs from biggest cluster saved to " + cluster_seqs_path)


def get_consensus_seq(cluster_seqs_path, seq_records, arg_threads):
    if (len(seq_records) == 1):
        return str(seq_records[0].seq.seq)

    from Bio.Align.Applications import ClustalwCommandline
    from Bio.Align.Applications import ClustalOmegaCommandline
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from Bio.Align import MultipleSeqAlignment

    aln_file = '.'.join(cluster_seqs_path.split('.')[:-1]) + "_aln.fasta"
    cmd = ClustalOmegaCommandline(infile=cluster_seqs_path, outfile=aln_file, force=True, threads=arg_threads)
    log.log("Run clustalOmega: " + str(cmd))
    stdout, stderr = cmd()
    log.log("Get Multiply alignment: " + aln_file)

    log.log("Start search for consensus monomer")
    align = AlignIO.read(aln_file, "fasta")
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
        mdist = min(mdist, get_distance(str(monomers_list[i].seq), new_monomer))
        mdist = min(mdist, get_distance(get_rc(str(monomers_list[i].seq)), new_monomer))
    return mdist


def init_first_run(args):
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
                              "Size of the largest cluster", "Number of deleted monomers at this iteration", "Radius", "Separation"])
    iter_id = 0
    return iter_id, summary_fw, summary_writer, ""


def get_lst_iter(outdir):
    iter_id = 1000
    while (not os.path.isdir(os.path.join(outdir, "iter_" + str(iter_id)))) and iter_id > -1:
        iter_id -= 1

    prev_id = iter_id - 1
    while (not os.path.isdir(os.path.join(outdir, "iter_" + str(prev_id)))) and prev_id > -1:
        prev_id -= 1

    return iter_id,  os.path.join(outdir, "iter_" + str(prev_id))


def init_continue(args):
    #switch working directory to output dir
    os.chdir(args.outdir)

    iter_id, prev_dir = get_lst_iter(args.outdir)

    args.monomers = os.path.join(args.outdir, "monomers.fa")
    args.sequences = os.path.join(args.outdir, 'sequence.fa')

    summary_path = os.path.join(args.outdir, "summary.csv")
    summary_fw = open(summary_path, "a")
    summary_writer = csv.writer(summary_fw)

    return iter_id, summary_fw, summary_writer, prev_dir


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


def need_update(args, monomer_record, monomer_resolved_block, iter_outdir, prev_dir):
    iter_id = int(iter_outdir.split('_')[-1])
    if iter_id == 0:
        return True
    else:
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

updated_monomers = []

def update_monomer(args, monomer_record, monomer_resolved_block, iter_outdir, prev_dir):
    global updated_monomers
    if need_update(args, monomer_record, monomer_resolved_block, iter_outdir, prev_dir):
        cluster_seqs_path = os.path.join(iter_outdir, monomer_record.id.split('/')[0] + "_seqs.fa")
        save_seqs(monomer_resolved_block, cluster_seqs_path)
        new_monomer = get_consensus_seq(cluster_seqs_path, monomer_resolved_block, args.threads)
        if reverse_monomer(args):
            new_monomer = rc(new_monomer)
        monomer_record.seq = Seq(new_monomer)
        updated_monomers.append(monomer_record.id)

    return monomer_record


def get_radius(new_monomer, max_cluster):
    mx_dist = 0
    for seq in max_cluster:
        mx_dist = max(mx_dist, get_distance(str(new_monomer), str(seq.seq.seq)))
    return mx_dist


def get_hybrid_len(main_mn, mn1, mn2, args):
    mn_identity = 100
    bst_res = (0, 0)
    for prfx in range(len(mn1.seq)):
        suffix = len(str(main_mn.seq)) - prfx
        hbr = str(mn1.seq)[:prfx] + str(mn2.seq)[len(mn2.seq) - suffix:]
        cur_identity = seq_identity(hbr, str(main_mn.seq))
        if mn_identity > cur_identity:
            mn_identity = cur_identity
            bst_res =  (prfx, suffix)
    if (mn_identity * 2 <= args.resDiv):
        return (bst_res[0], bst_res[1], mn_identity)
    return (0, 0, 100)

def detect_hybrid_mn(monomers_list, res_tsv, args):
    global TotalMonomerBlocks
    global FreqCeiling

    cnt_mn = {}
    TotalMonomerBlocks = 0

    with open(res_tsv, "r") as f:
        csv_reader = csv.reader(f, delimiter='\t')
        for row in csv_reader:
            if row[2] == "start":
                continue
            identity = float(row[4])
            if row[-1] == '?':
                continue
            TotalMonomerBlocks += 1
            # print("Identity: ", identity)
            if identity >= 100 - args.resDiv:
                if row[1][-1] == "'":
                    row[1] = row[1][:-1]

                if row[1] not in cnt_mn:
                    cnt_mn[row[1]] = 0
                cnt_mn[row[1]] += 1

    for i in range(len(monomers_list)):
        log.log("hybrid detection for monomer#" + str(i) + " out of " + str(len(monomers_list)))
        bst_hyber_score = 100
        bst_hyber = (0, 0)

        for j in range(len(monomers_list)):
            for g in range(len(monomers_list)):
                if i == j or i == g:
                    continue

                if cnt_mn[monomers_list[j].id] < TotalMonomerBlocks/FreqCeiling:
                    continue

                if cnt_mn[monomers_list[g].id] < TotalMonomerBlocks/FreqCeiling:
                    continue

                hybrid_res = get_hybrid_len(monomers_list[i], monomers_list[j], monomers_list[g], args)
                if (hybrid_res[0] != 0 and hybrid_res[2] < bst_hyber_score):
                    bst_hyber_score = hybrid_res[2]
                    bst_hyber = (j, g)

        if bst_hyber_score < 100:
            cnt1, cnt2, scr = get_hybrid_len(monomers_list[i], monomers_list[bst_hyber[0]], monomers_list[bst_hyber[1]], args)
            mn1 = monomers_list[bst_hyber[0]].id.split('_')[1]
            mn2 = monomers_list[bst_hyber[1]].id.split('_')[1]

            old_name = monomers_list[i].id
            if (mn1 != mn2):
                monomers_list[i].id += "_hybrid_" + mn1 + "(" + str(cnt1) + ")_" + mn2 + "(" + str(cnt2) + ")"
            else:
                monomers_list[i].id += "_variant_" + mn1

            #monomers_list[i].id += "_hybrid_" + mn1 + "_" + mn2
            cnt_mn[monomers_list[i].id] = cnt_mn[old_name]


    return monomers_list


def final_iteration(args, sd_script_path, monomers_list):
    log.log("====== Start final iteration======")

    iter_outdir = os.path.join(args.outdir, "final")
    if not os.path.exists(iter_outdir):
        os.makedirs(iter_outdir)
    os.chdir(iter_outdir)

    local_monmers_path = os.path.join(iter_outdir, "monomers.fa")
    shutil.copyfile(args.monomers, local_monmers_path)

    # run string decomposer
    sys_call(["python3", sd_script_path, args.sequences, local_monmers_path, "-t", str(args.threads), "--fast", "--ed_thr", "20"])
    log.log("String decomposer is complete. Results save in: " + iter_outdir)

    # parse output csv file
    res_tsv = os.path.join(iter_outdir, "final_decomposition.tsv")

    log.log("Start detect hybrids")
    monomers_list = detect_hybrid_mn(monomers_list, res_tsv, args)

    with open(args.monomers, "w") as fa:
        for record in monomers_list:
            SeqIO.write(record, fa, "fasta")


def get_unresolved_blocks(res_tsv, monomers_list, args):
    log.log("==Get unresolved blocls==")
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

    return unresolved_blocks, monomer_resolved, non_monomeric_cnt, resolved_cnt


def delete_unused_monomers(monomers_list, monomer_resolved):
    log.log("==Delte unused monomers==")
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
    return deleted_cnt, monomers_list


def update_all_monomers(monomers_list, args, monomer_resolved, iter_outdir, prev_dir):
    global updated_monomers
    updated_monomers = []
    for i in range(len(monomers_list)):
        log.log("==== Update monomer " + str(monomers_list[i].id))
        set_blocks_seq(args.sequences, monomer_resolved[monomers_list[i].id])
        monomers_list[i] = update_monomer(args, monomers_list[i], monomer_resolved[monomers_list[i].id], iter_outdir, prev_dir)


def chr_statistic(blocks):
    cnt = {}
    for block in blocks:
        if block.read_name not in cnt:
            cnt[block.read_name] = 0
        cnt[block.read_name] += 1

    inclds = []
    for key in cnt:
        inclds.append((-cnt[key], key.split(':')[0]))

    inclds.sort()

    cnts = []
    crms = []
    for incl in inclds:
        cnts.append(-incl[0])
        crms.append(incl[1])
    return cnts, crms


def main():
    global updated_monomers
    log.log("Start Monomer Inference")
    args = parse_args()
    # get path to current script
    current_script_path = os.path.abspath(os.path.dirname(__file__))

    # path to the run_decomposer script
    sd_script_path = os.path.join(current_script_path, "..", "sd", "run_decomposer.py")
    log.log("Path to run_decomposer: " + sd_script_path)

    # create output_dir if not exists
    args.outdir = os.path.abspath(args.outdir)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if args.restart == False:
        iter_id, summary_fw, summary_writer, prev_dir = init_first_run(args)
    else:
        iter_id, summary_fw, summary_writer, prev_dir = init_continue(args)

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

        if not os.path.exists(os.path.join(iter_outdir, "final_decomposition.tsv")):
            extra_options = []
            if updated_monomers != []:
                upmn_path = os.path.join(iter_outdir, "updated_monomers")
                with open(upmn_path, "w") as f:
                    for u_mn in updated_monomers:
                        f.write(u_mn + "\n")
                extra_options = ["--new-monomers", upmn_path]
            # run string decomposer
            sys_call(["python3", sd_script_path, args.sequences, local_monmers_path, "-t", str(args.threads), "--fast"])
            log.log("String decomposer is complete. Results save in: " + iter_outdir)

        # parse output csv file
        res_tsv = os.path.join(iter_outdir, "final_decomposition.tsv")
        dist_to_monomers = 0
        radius = 0

        unresolved_blocks, monomer_resolved, \
        non_monomeric_cnt, resolved_cnt = get_unresolved_blocks(res_tsv, monomers_list, args)

        deleted_cnt, monomers_list = delete_unused_monomers(monomers_list, monomer_resolved)

        update_all_monomers(monomers_list, args, monomer_resolved, iter_outdir, prev_dir)

        log.log("Number of unresolved monomer block: " + str(len(unresolved_blocks)))
        log.log("Number of resolved blocks: " + str(resolved_cnt))
        log.log("Number of non-monomeric blocks: " + str(non_monomeric_cnt))
        log.log("Number of deleted monomers: " + str(deleted_cnt))

        if len(unresolved_blocks) > 0:
            set_blocks_seq(args.sequences, unresolved_blocks)
            max_cluster = clustering(unresolved_blocks, args)
            if (len(max_cluster[0]) < args.min_clst):
                monomer_set_complete = True
        else:
            monomer_set_complete = True

        if not monomer_set_complete:
            for i in range(len(max_cluster)):
                log.log("===DETECT CONSENSUS FOR CLUSTER#" + str(i) + "===")
                cluster_seqs_path = os.path.join(iter_outdir, "cluster_seq_" + str(i) + ".fa")
                save_seqs(max_cluster[i], cluster_seqs_path)

                new_monomer = get_consensus_seq(cluster_seqs_path, max_cluster[i], args.threads)
                radius = get_radius(new_monomer, max_cluster[i])

                if reverse_monomer(args):
                    new_monomer = rc(new_monomer)

                dist_to_monomers = get_dist_to_exists_monomers(monomers_list, new_monomer)
                log.log("Min Distance to exsisting monomers: " + str(dist_to_monomers))

                new_monomer_record = SeqRecord(Seq(new_monomer), id="mn_" + str(iter_id + i), description="")
                monomers_list.append(new_monomer_record)
                updated_monomers.append(new_monomer_record.id)
                ch_cnt, crms = chr_statistic(max_cluster[i])
                summary_writer.writerow(
                    [str(iter_id + i), str(resolved_cnt), str(len(unresolved_blocks)), str(non_monomeric_cnt),
                     str(len(max_cluster[i])), str(deleted_cnt), str(radius),
                     str(dist_to_monomers), str(ch_cnt), str(crms)])
        else:
            summary_writer.writerow(
                [str(iter_id), str(resolved_cnt), str(len(unresolved_blocks)), str(non_monomeric_cnt),
                 str(len(max_cluster[0])), str(deleted_cnt), str(radius), str(dist_to_monomers), "-", "-"])

        with open(args.monomers, "w") as fa:
            for record in monomers_list:
                SeqIO.write(record, fa, "fasta")

        summary_fw.flush()
        prev_dir = iter_outdir
        iter_id += len(max_cluster)

    final_iteration(args, sd_script_path, monomers_list)
    summary_fw.close()

if __name__ == "__main__":
    main()
