from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def load_fasta(filename, tp = "list"):
    if tp == "map":
        records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    else:
        records = list(SeqIO.parse(filename, "fasta"))
    return records

def load_decomposition(filename):
    reads_mapping = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, sstart, send, idnt, q = ln.strip().split("\t")[:6]
            if sseqid not in reads_mapping:
                    reads_mapping[sseqid] = []
            s, e, idnt = int(sstart), int(send), float(idnt)
            rev = False
            if qseqid.endswith("'"):
                rev = True
                qseqid = qseqid[:-1]
            reads_mapping[sseqid].append({"qid": qseqid, "s": s, "e": e, "rev": rev, "idnt": idnt, "q": q})
    for r in reads_mapping:
        if len(reads_mapping[r]) > 0 and reads_mapping[r][0]["rev"]:
            reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: (-x["e"], x["s"]))
        else:
            reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: (x["s"], -x["e"]))

    return reads_mapping

class StatisticsCounter:
    def __init__(self, filename_alns, filename_reads, filename_monomers):
        self.reads = load_fasta(filename_reads, "map")
        self.monomers = load_fasta(filename_monomers, "map")
        self.hits = load_decomposition(filename_alns)
        self.NOGAP = 10

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

        symbol = {}
        for m in self.monomers:
            symbol[m] = int(m.split("_")[1])

        s_r = {}
        for s in symbol:
            s_r[symbol[s]] = s

        res = ""
        for m in sorted([symbol[k] for k in order_lst.keys()]):
            ans = []
            for m2 in sorted([symbol[kk] for kk in order_lst[s_r[m]].keys()]):
                ans.append(str(order_lst[s_r[m]][s_r[m2]]))
            res += "\t".join(ans) + "\n"
        return res

    def identity_histogram(self):
        pass

    def covered(self, high_qual = False):
        sum_len = 0
        sum_filled = 0
        avg_mapping_identity, mapping_num = 0, 0
        for r in self.hits:
            read_len = len(self.reads[r].seq)
            mappings = self.hits[r]
            ar = [0 for _ in range(read_len)]
            for m in mappings:
                if not high_qual or m["q"] == "+":
                    s, e = m["s"], m["e"]
                    for i in range(s-1, e):
                        ar[i] = 1
                    avg_mapping_identity += m["idnt"]
                    mapping_num += 1
            cnt = sum(ar)
            sum_len += read_len
            sum_filled += cnt
        return str(sum_filled*100//sum_len)+"%", str(int(avg_mapping_identity//mapping_num))+"%"

    def abnormal_parts(self, high_qual = False):
        gaps = 0
        abnormal_order = 0
        hits = self.hits
        if high_qual:
            high_qual_mp = {}
            for r in self.hits:
                high_qual_mp[r] = []
                for h in self.hits[r]:
                    if h["q"] == "+":
                        high_qual_mp[r].append(h)
            hits = high_qual_mp
        for r in hits:
            for i in range(len(hits[r]) - 1):
                h1 = hits[r][i]
                h2 = hits[r][i + 1]
                if h1["rev"] == h2["rev"]:
                    rev = h1["rev"]
                    if (not rev and abs(h2["s"] - h1["e"]) < self.NOGAP) or (rev and abs(h2["e"] - h1["s"]) < self.NOGAP):
                        s1, s2 = int(h1["qid"].split("_")[1]), int(h2["qid"].split("_")[1])
                        if s1 + 1 != s2 and not (s1 == len(self.monomers) and s2 == 1):
                            abnormal_order += 1 
                else:
                    gaps += 1
        return gaps, abnormal_order
            
    def fast_consensus(self):
        pass

    def generate_stats_table(self):
        columns = ["                 ", "Covered", "Avg identity", "Gaps", "Abnormal mapping pairs"]
        rows = ["Total set of mappings", "High quality mappings"]
        columns_mp = {}
        for c in columns:
            columns_mp[c] = {}
            for r in rows:
                columns_mp[c][r] = ""
        columns_mp["                 "]["Total set of mappings"], columns_mp["                 "]["High quality mappings"] = rows[0], rows[1]
        columns_mp["Covered"]["Total set of mappings"], columns_mp["Avg identity"]["Total set of mappings"] = self.covered()
        columns_mp["Covered"]["High quality mappings"], columns_mp["Avg identity"]["High quality mappings"] = self.covered(True)
        columns_mp["Gaps"]["Total set of mappings"], columns_mp["Abnormal mapping pairs"]["Total set of mappings"] = self.abnormal_parts()
        columns_mp["Gaps"]["High quality mappings"], columns_mp["Abnormal mapping pairs"]["High quality mappings"] = self.abnormal_parts(True)
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
        with open(filename_out + "_stats.txt", "w") as fout:
            fout.write("Reads number: " + str(len(self.reads)) + "\n")
            fout.write("Reads length: " + str(self.reads_len()) + "\n")
            fout.write("General alignment statistics\n")
            fout.write(self.generate_stats_table())
            fout.write("\n")
            fout.write("Monomers transition matrix\n")
            fout.write(self.transition_matrix())
        print("Statistics saved to " + filename_out + "_stats.txt")

