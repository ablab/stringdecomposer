from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import edlib

def cnt_identity(lst, cur_mode = "NW"):
    if len(str(lst[0])) == 0:
        return 0
    if len(str(lst[1])) == 0:
        return 0
    result = edlib.align(str(lst[0]), str(lst[1]), mode=cur_mode, task="distance")
    return 100 - result["editDistance"]*100//max(len(str(lst[0])), len(str(lst[1])))

def load_fasta(filename, tp = "list"):
    if tp == "map":
        records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    else:
        records = list(SeqIO.parse(filename, "fasta"))
    return records

def load_decomposition(filename, reads, monomers, identity_threshold):
    reads_mapping = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, sstart, send, idnt = ln.strip().split("\t")[:5]
            if sseqid not in reads_mapping:
                    reads_mapping[sseqid] = []
            s, e, idnt = int(sstart), int(send), float(idnt)
            rev = False
            if qseqid.endswith("'"):
                rev = True
                qseqid = qseqid[:-1]
            if rev:
                idnt = cnt_identity([reads[sseqid].seq[s: e + 1], monomers[qseqid].seq.reverse_complement()])
            else:
                idnt = cnt_identity([reads[sseqid].seq[s: e + 1], monomers[qseqid].seq])
            if idnt > identity_threshold:
                reads_mapping[sseqid].append({"qid": qseqid, "s": s, "e": e, "rev": rev, "idnt": idnt})

    new_mapping = {}
    bad_reads_num = 0
    for r in reads_mapping:
        if len(reads_mapping[r]) > 36:
            cnt_rev = 0
            for i in range(len(reads_mapping[r])):
                h1 = reads_mapping[r][i]
                if h1["rev"]:
                    cnt_rev += 1
            if cnt_rev == 0 or cnt_rev == len(reads_mapping[r]):
                new_mapping[r] = reads_mapping[r]
            elif cnt_rev < 0.1*len(reads_mapping[r]):
                new_mapping[r] = []
                for i in range(len(reads_mapping[r])):
                    h1 = reads_mapping[r][i]
                    if not h1["rev"]:
                        new_mapping[r].append(h1)
            elif cnt_rev > 0.9*len(reads_mapping[r]):
                new_mapping[r] = []
                for i in range(len(reads_mapping[r])):
                    h1 = reads_mapping[r][i]
                    if h1["rev"]:
                        new_mapping[r].append(h1)
            else:
                bad_reads_num += 1
        else:
            bad_reads_num += 1

    print("Bad reads " + str(bad_reads_num))
    for r in new_mapping:
        if len(new_mapping[r]) > 0 and new_mapping[r][0]["rev"]:
            new_mapping[r] = sorted(new_mapping[r], key=lambda x: (-x["e"], x["s"]))
        else:
            new_mapping[r] = sorted(new_mapping[r], key=lambda x: (x["s"], -x["e"]))
    return new_mapping

def identify_monomer_order(filename):
    monomers_lst = load_fasta(filename)
    ordered = []
    for m in monomers_lst:
        ordered.append(m.name)
    order_map = {}
    for i in range(len(ordered)):
        order_map[ordered[i]] = i
    return ordered, order_map

class StatisticsCounter:
    def __init__(self, filename_alns, filename_reads, filename_monomers, identity_threshold):
        self.reads = load_fasta(filename_reads, "map")
        self.monomers = load_fasta(filename_monomers, "map")
        self.hits = load_decomposition(filename_alns, self.reads, self.monomers, identity_threshold)
        self.NOGAP = 10
        self.monomers_ordered, self.monomer_order_map = identify_monomer_order(filename_monomers)
        self.identity_threshold = identity_threshold

    def transition_matrix(self):
        order_lst = {}
        for m in self.monomers:
            order_lst[m] = {} 
            for m2 in self.monomers:
                order_lst[m][m2] = 0
        for r in sorted(self.hits.keys()):
            if len(self.reads) == 0 or r in self.reads:  
                for i in range(len(self.hits[r]) - 1):
                    h1 = self.hits[r][i]
                    h2 = self.hits[r][i + 1]
                    if h1["rev"] == h2["rev"]:
                        rev = h1["rev"]
                        if  (not rev and abs(h2["s"] - h1["e"]) < self.NOGAP) or (rev and abs(h2["e"] - h1["s"]) < self.NOGAP):
                            order_lst[h1["qid"]][h2["qid"]] += 1


        res = ""
        for m in self.monomers_ordered:
            ans = []
            for m2 in self.monomers_ordered:
                ans.append(str(order_lst[m][m2]))
            res += "\t".join(ans) + "\n"
        return res

    def identity_histogram(self):
        pass

    def covered(self):
        sum_len = 0
        sum_filled = 0
        avg_mapping_identity, mapping_num = 0, 0
        for r in self.hits:
            read_len = len(self.reads[r].seq)
            mappings = self.hits[r]
            ar = [0 for _ in range(read_len)]
            for m in mappings:
                s, e = m["s"], m["e"]
                for i in range(s-1, e):
                    ar[i] = 1
                avg_mapping_identity += m["idnt"]
                mapping_num += 1
            cnt = sum(ar)
            sum_filled += cnt
        for r in self.reads:
            sum_len += len(self.reads[r].seq)
        return str(sum_filled*100//sum_len)+"%", str(int(avg_mapping_identity//mapping_num))+"%"

    def abnormal_parts(self):
        gaps = 0
        abnormal_order = 0
        for r in self.hits:
            for i in range(len(self.hits[r]) - 1):
                h1 = self.hits[r][i]
                h2 = self.hits[r][i + 1]
                rev = h1["rev"]
                if (not rev and abs(h2["s"] - h1["e"]) < self.NOGAP) or (rev and abs(h2["e"] - h1["s"]) < self.NOGAP):
                    s1, s2 = self.monomer_order_map[h1["qid"]], self.monomer_order_map[h2["qid"]]
                    if s1 + 1 != s2 and not (s1 == len(self.monomers) - 1 and s2 == 0):
                        abnormal_order += 1 
                else:
                    gaps += 1
        return gaps, abnormal_order
            
    def fast_consensus(self):
        pass

    def generate_stats_table(self):
        columns = ["                 ", "Covered", "Avg identity", "Gaps", "Abnormal mapping pairs"]
        rows = ["Total set of mappings"]
        columns_mp = {}
        for c in columns:
            columns_mp[c] = {}
            for r in rows:
                columns_mp[c][r] = ""
        columns_mp["                 "]["Total set of mappings"]= rows[0]
        columns_mp["Covered"]["Total set of mappings"], columns_mp["Avg identity"]["Total set of mappings"] = self.covered()
        columns_mp["Gaps"]["Total set of mappings"], columns_mp["Abnormal mapping pairs"]["Total set of mappings"] = self.abnormal_parts()
        res = "\t".join(columns) + "\n"
        print(columns_mp)
        for r in rows:
            res_lst = []
            for c in columns:
                res_lst.append(str(columns_mp[c][r]))
            res += "\t".join(res_lst) + "\n"
        return res

    def reads_len(self):
        res = 0
        for r in self.reads:
            res += len(self.reads[r].seq)
        return res

    def save_stats(self, filename_out):
        print("Calculating resulting alignments statistics...")
        with open(filename_out + "_stats" + str(self.identity_threshold) + ".txt", "w") as fout:
            fout.write("Reads number: " + str(len(self.reads)) + "\n")
            fout.write("Reads length: " + str(self.reads_len()) + "\n")
            fout.write("General alignment statistics\n")
            fout.write(self.generate_stats_table())
            fout.write("\n")
            fout.write("Monomers transition matrix\n")
            fout.write(self.transition_matrix())
        print("Statistics saved to " + filename_out + "_stats" + str(self.identity_threshold) + ".txt")

