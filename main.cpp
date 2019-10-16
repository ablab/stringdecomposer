#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cstring>
#include <omp.h>

using namespace std;

struct ReadId {
    string name;
    int id;

    ReadId(string name_): name(name_), id(-1) {}
    ReadId(string name_, int id_): name(name_), id(id_) {}
};

struct Seq {
    ReadId read_id;
    string seq;

    Seq(string name_, string seq_): seq(seq_), read_id(ReadId(name_)) {}

    Seq(string name_, string seq_, int id_): seq(seq_), read_id(ReadId(name_, id_)) {}

    size_t size() { return seq.size();}
};

struct MonomerAlignment {
    string monomer_name;
    string read_name;
    int start_pos;
    int end_pos;
    float identity;
    bool best;

    MonomerAlignment() {}

    MonomerAlignment(const MonomerAlignment & m_aln){
        monomer_name = m_aln.monomer_name;
        read_name = m_aln.read_name;
        start_pos = m_aln.start_pos;
        end_pos = m_aln.end_pos;
        identity = m_aln.identity;
        best = m_aln.best;
    }

    MonomerAlignment(string monomer_name_, string read_name_, int start_pos_, int end_pos_, float identity_, bool best_)
    : monomer_name(monomer_name_), read_name(read_name_), start_pos(start_pos_), end_pos(end_pos_), identity(identity_), best(best_) {}
};

bool sortby1(const pair<int, vector<MonomerAlignment>> &a
             , const pair<int, vector<MonomerAlignment>> &b) {
    return (a.first < b.first);
}

class MonomersAligner {

public:
    MonomersAligner(vector<Seq> &monomers, string output_file_name):monomers_(monomers) {
        output_file_best_.open(output_file_name + ".tsv", std::ofstream::out);
        output_file_.open(output_file_name + "_alt.tsv", std::ofstream::out);
    }

    void AlignReadsSet(vector<Seq> &reads) {
        vector<Seq> new_reads;
        vector<int> save_steps;
        for (auto r: reads) {
            int cnt = 0;
            //cout << r.seq.size() << endl;
            for (int i = 0; i < r.seq.size(); i += PART_SZ) {
                if ((int) r.seq.size() - i >= 200) {  
                    Seq seq = Seq(r.read_id.name, r.seq.substr(i, min(PART_SZ + 200, (int) r.seq.size() - i) ), i );
                    new_reads.push_back(seq);
                    ++ cnt;
                }
            }
            save_steps.push_back(cnt);
        }
        cout << "Prepared reads\n";
        
        int start = 0;
        for (int p = 0; p < save_steps.size(); ++ p){
            vector<MonomerAlignment> batch;
            vector<pair<int, vector<MonomerAlignment>>> subbatches;
            #pragma omp parallel for schedule(guided, 50) num_threads(16)
            for (int j = start; j < start + save_steps[p]; ++ j) {
                vector<MonomerAlignment> aln = AlignPartClassicDP(new_reads[j]);
                #pragma omp critical(aligner) 
                {
                    subbatches.push_back(pair<int, vector<MonomerAlignment>> (j, aln));
                }
            }
            sort(subbatches.begin(), subbatches.end(), sortby1);
            for (int j = start; j < start + save_steps[p]; ++ j) {
                int read_index = subbatches[j - start].first;
                for (auto a: subbatches[j - start].second) {
                    MonomerAlignment new_m_aln(a.monomer_name, a.read_name, 
                                                new_reads[read_index].read_id.id + a.start_pos, new_reads[read_index].read_id.id + a.end_pos, 
                                                a.identity, a.best);
                    batch.push_back(new_m_aln);
                }
            }
            cout << "Aligned " << batch[0].read_name << endl;
            batch = PostProcessing(batch);
            SaveBatch(batch);
            start += save_steps[p];
        }           
    }

    ~MonomersAligner() {
        output_file_.close();
        output_file_best_.close();
    }

private:

    void PrecalculateMonomerAlignment(){

    }

    void PrecalculateMonomerEdlib(){

    }

    vector<MonomerAlignment> AlignPartFitting(Seq &read) {

    }

    vector<MonomerAlignment> AlignPartClassicDP(Seq &read) {
        int ins = -1;
        int del = -1;
        int match = 1;
        int mismatch = -1;
        int INF = -1000000;
        int monomers_num = (int) monomers_.size();
        vector<vector<vector<int>>> dp(read.seq.size());
        //cout << dp.size() << endl;
        for (int i = 0; i < read.seq.size(); ++ i) {
            for (auto m: monomers_) {
                dp[i].push_back(vector<int>(m.seq.size()));
                for (int k = 0; k < m.seq.size(); ++ k) {
                    dp[i][dp[i].size() - 1][k] = INF; 
                }
            }
            dp[i].push_back(vector<int>(1));
            dp[i][monomers_num][0] = INF;
        }

        for (int j = 0; j < monomers_.size(); ++ j) {
            Seq m = monomers_[j];
            if (m.seq[0] == read.seq[0]) {
                dp[0][j][0] = match;
            } else {
                dp[0][j][0] = mismatch;
            }
            for (int k = 1; k < m.seq.size(); ++ k) {
                dp[0][j][k] = dp[0][j][k-1] + del;
            }
        }
        //dp[0][monomers_num][0] = 0;
        //cout << "Init dp " << read.seq.size() << endl;
        for (int i = 1; i < read.seq.size(); ++ i) {
            //dp[i][monomers_num][0] = dp[i-1][monomers_num][0];
            for (int j = 0; j < monomers_.size(); ++ j) {
                dp[i][monomers_num][0] = max(dp[i][monomers_num][0], dp[i-1][j][monomers_[j].size() - 1]);
            }
            for (int j = 0; j < monomers_.size(); ++ j) {
                for (int k = 0; k < monomers_[j].size(); ++ k) {
                    int score = INF;
                    int mm_score = monomers_[j].seq[k] == read.seq[i] ? match: mismatch;
                    if (dp[i][monomers_num][0] > INF) {
                        score = max(score, dp[i][monomers_num][0] + mm_score + k*del);
                    }
                    if (k > 0) {
                        if (dp[i-1][j][k-1] > INF) {
                            score = max(score, dp[i-1][j][k-1] + mm_score);
                        }
                        if (dp[i-1][j][k] > INF) {
                            score = max(score, dp[i-1][j][k] + ins);
                        }
                        if (dp[i][j][k-1] > INF) {
                            score = max(score, dp[i][j][k-1] + del);
                        }
                    }
                    dp[i][j][k] = score;
                }
            }
        }
        //cout << "Counted dp \n";
        int max_score = INF; //dp[read.seq.size()-1][monomers_num][0];
        int best_m = monomers_num;
        for (int j = 0; j < monomers_.size(); ++ j) {
            if (max_score < dp[read.seq.size()-1][j][monomers_[j].size() -1] ) {
                max_score = dp[read.seq.size()-1][j][monomers_[j].size() -1];
                best_m = j;
            }
        }
        //cout << max_score << endl;
        vector<MonomerAlignment> ans;
        int i = read.seq.size() - 1;
        int j = best_m;
        int k = dp[i][j].size() - 1;
        MonomerAlignment cur_aln;
        while (i >= 0) {
            if (k == dp[i][j].size() - 1 && j != monomers_num) {
                cur_aln = MonomerAlignment(monomers_[j].read_id.name, read.read_id.name, i, i, dp[i][j][k], true); 
            } 
            if (j == monomers_num) {
                if (i != 0) {
                    for (int p = 0; p < dp[i - 1].size(); ++ p) {
                        if (dp[i - 1][p][dp[i - 1][p].size() - 1] == dp[i][j][k]) {
                            -- i; 
                            j = p;
                            k = dp[i][j].size() - 1;
                            break;
                        } 
                    }
                } else {
                    -- i;
                }
            } else {
                if (k != 0 && dp[i][j][k] == dp[i][j][k-1] + del) {
                    --k;
                } else {
                    if (i != 0 && dp[i][j][k] == dp[i-1][j][k] + ins) {
                        --i;
                    } else{
                        int mm_score = monomers_[j].seq[k] == read.seq[i] ? match: mismatch;
                        if (i != 0 && k != 0 && dp[i][j][k] == dp[i-1][j][k-1] + mm_score) {
                            --i; --k;
                        } else {
                            if (i != 0 && dp[i][monomers_num][0] + k*del + mm_score ==  dp[i][j][k]) {
                                cur_aln.start_pos = i;
                                cur_aln.identity = cur_aln.identity - dp[i][monomers_num][0];
                                ans.push_back(cur_aln);
                                j = monomers_num; k = 0;
                            } else {
                                cur_aln.start_pos = i;
                                ans.push_back(cur_aln);
                                --i;
                            }
                        }
                    }
                }
            }
        }
        reverse(ans.begin(), ans.end());
        return ans;
    }

    void SaveBatch(vector<MonomerAlignment> &batch) {
        for (auto a: batch) {
            string s = a.read_name + "\t" 
                       + a.monomer_name + "\t" 
                       + to_string(a.start_pos) + "\t"
                       + to_string(a.end_pos) + "\t"
                       + to_string(a.identity);
            if (a.best) {
                output_file_best_ << s << "\n";
                s += "\t*";    
            }
            output_file_ << s << "\n";
        }
    }

    vector<MonomerAlignment> PostProcessing(vector<MonomerAlignment> &batch) {
        vector<MonomerAlignment> res;
        res.push_back(batch[0]);
        for (size_t i = 1; i < batch.size(); ++ i) {
            bool add = true;
            for (size_t j = (size_t) max((int) 0, (int) i - 3); j < i; ++ j) {
                if ((batch[i].end_pos - batch[j].end_pos < 50 || batch[i].start_pos - batch[j].start_pos < 50) && batch[i].identity < batch[j].identity) {
                    add = false;
                }
            }
            for (size_t j = i + 1; j < (size_t) min((int) i + 4, (int) batch.size()); ++ j) {
                if ((batch[j].end_pos - batch[i].end_pos < 50 || batch[j].start_pos - batch[i].start_pos < 50) && batch[i].identity < batch[j].identity) {
                    add = false;
                }
            }
            if (add) { res.push_back(batch[i]);}
            // if (batch[i].monomer_name == batch[i - 1].monomer_name && batch[i].end_pos - batch[i-1].end_pos < 150) {
            //     if (batch[i].identity > batch[i-1].identity) {
            //         res[res.size() - 1] = batch[i];
            //     }
            // } else {
            //     res.push_back(batch[i]);
            // }
        }
        return res;
    }

    vector<Seq> monomers_;
    const int SAVE_STEP = 1;
    const int PART_SZ = 5000;
    ofstream output_file_;
    ofstream output_file_best_;
};

vector<Seq> load_fasta(string filename) {
    std::ifstream input_file;
    input_file.open(filename, std::ifstream::in);
    string s;
    vector<Seq> seqs;
    while (getline(input_file, s)) {
        if (s[0] == '>') {
            seqs.push_back(Seq(s.substr(1, s.size() - 1), ""));
        } else {
            seqs[seqs.size()-1].seq += s;
        }
    }
    return seqs;
}

string reverse_complement(string &s){
    string res = "";
    map<char, char> rc = {{'A', 'T'}, {'T', 'A'}, {'G','C'}, {'C','G'}};
    for (int i = (int) s.size() - 1; i >= 0; --i){
        res += rc[s[i]];
    }
    return res;
}

void add_reverse_complement(vector<Seq> &monomers) {
    vector<Seq> rev_c_monomers;
    for (auto s: monomers) {
        rev_c_monomers.push_back(Seq(s.read_id.name + "_rc", reverse_complement(s.seq)));
    }
    monomers.insert(monomers.end(), rev_c_monomers.begin(), rev_c_monomers.end());
    return;
}


int main() {
    vector<Seq> reads = load_fasta("../chrX/results/centromeric_reads/centromeric_reads.fasta");
    vector<Seq> monomers = load_fasta("../chrX/DXZ1_rc_star_monomers.fasta");
    add_reverse_complement(monomers);
    MonomersAligner monomers_aligner(monomers, "cpp_version");
    monomers_aligner.AlignReadsSet(reads);
}