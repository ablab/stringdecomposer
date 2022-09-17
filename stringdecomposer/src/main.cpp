#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <iterator>
#include <omp.h>

#include "edlib.h"

using namespace std;

struct ReadId {
    string name;
    int id = -1;

    ReadId(string name_): name(name_) {}
    ReadId(string name_, int id_): name(name_), id(id_) {}
};

struct Seq {
    ReadId read_id;
    string seq;

    Seq(string name_, string seq_): read_id(ReadId(name_)), seq(seq_) {
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    }

    Seq(string name_, string seq_, int id_): read_id(ReadId(name_, id_)), seq(seq_)  {}

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

    MonomerAlignment(string monomer_name_, string read_name_, int start_pos_, int end_pos_, float identity_, bool best_)
    : monomer_name(monomer_name_), read_name(read_name_), start_pos(start_pos_), end_pos(end_pos_), identity(identity_), best(best_) {}
};

bool sortby1(const pair<int, vector<MonomerAlignment>> &a
             , const pair<int, vector<MonomerAlignment>> &b) {
    return (a.first < b.first);
}

class MonomersAligner {

public:
    MonomersAligner(vector<Seq> &monomers, int ins = -1, int del = -1, int mismatch = -1, int match = 1)
      : monomers_(monomers),
        ins_(ins),
        del_(del),
        mismatch_(mismatch),
        match_(match) {
    }

    void AlignReadsSet(vector<Seq> &reads, int threads, int part_size, int ed_thr, int overlap = 500) {
        vector<Seq> new_reads;
        vector<int> save_steps;
        for (const auto & r: reads) {
            int cnt = 0;
            //cout << r.seq.size() << endl;
            for (size_t i = 0; i < r.seq.size(); i += part_size) {
                if ((int) r.seq.size() - i >= overlap || r.seq.size() < overlap) {
                    Seq seq = Seq(r.read_id.name, r.seq.substr(i, min(part_size + overlap, static_cast<int>(r.seq.size() - i)) ), i );
                    new_reads.push_back(seq);
                    ++ cnt;
                }
            }
            save_steps.push_back(cnt);
        }
        cerr << "Prepared reads\n";
        
        size_t start = 0, p = 0;
        int step = threads*2;
        vector<pair<int, vector<MonomerAlignment>>> subbatches;
        for (size_t i = 0; i < new_reads.size(); i += step) {
            #pragma omp parallel for num_threads(threads)
            for (size_t j = i; j < min(i + step, new_reads.size()); ++ j) {
                std::vector<MonomerAlignment> aln;
                if (ed_thr > -1) {
                    std::vector<Seq> filter_monomers = FilterMonomersForRead(new_reads[j], ed_thr);
                    aln = AlignPartClassicDP(new_reads[j], filter_monomers);
                } else {
                    aln = AlignPartClassicDP(new_reads[j], monomers_);
                }
                
                #pragma omp critical(aligner)
                {
                    subbatches.push_back(pair<int, vector<MonomerAlignment>> (j, aln));
                }
            }
            sort(subbatches.begin() + i, subbatches.begin() + min(i + step, new_reads.size()), sortby1);
            while (p < save_steps.size() && start + save_steps[p] <= subbatches.size()) {
                vector<MonomerAlignment> batch;
                for (size_t j = start; j < start + save_steps[p]; ++ j) {
                    int read_index = subbatches[j].first;
                    for (auto a: subbatches[j].second) {
                        MonomerAlignment new_m_aln(a.monomer_name, a.read_name,
                                                    new_reads[read_index].read_id.id + a.start_pos, new_reads[read_index].read_id.id + a.end_pos,
                                                    a.identity, a.best);
                        batch.push_back(new_m_aln);
                    }
                }
                cerr << (p + 1) * 100/save_steps.size() << "%: Aligned " << batch[0].read_name << endl;
                batch = PostProcessing(batch);
                SaveBatch(batch);
                start += save_steps[p];
                ++ p;
            }
        }           
    }

    ~MonomersAligner() {
    }

private:
    double MonomerEditDistance(Seq& monomer, Seq& read) {
        EdlibAlignResult result = edlibAlign(monomer.seq.c_str(), monomer.seq.size(), read.seq.c_str(), read.seq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
        double res = result.editDistance;
        edlibFreeAlignResult(result);
        return res;
    }

    std::vector<Seq> FilterMonomersForRead(Seq& read, int ed_thr) {
        std::vector<Seq> monomers_for_read;
        std::vector<std::pair<double, int>> mn_edit;
        for (size_t i = 0; i < monomers_.size(); ++i) {
            mn_edit.push_back(std::make_pair(MonomerEditDistance(monomers_[i], read), i));
        }
        std::sort(mn_edit.begin(), mn_edit.end());
        monomers_for_read.push_back(monomers_[mn_edit[0].second]);
        for (size_t i = 1; i < mn_edit.size(); ++i) {
            if (mn_edit[i].first <= ed_thr) {
                monomers_for_read.push_back(monomers_[mn_edit[i].second]);
            }
        }
        return monomers_for_read;
    }

    vector<MonomerAlignment> AlignPartClassicDP(Seq &read, std::vector<Seq>& monomers) {
        int ins = ins_;
        int del = del_;
        int match = match_;
        int mismatch = mismatch_;
        int INF = -1000000;
        int monomers_num = (int) monomers.size();
        vector<vector<vector<long long>>> dp(read.seq.size());
        //cout << dp.size() << endl;
        for (size_t i = 0; i < read.seq.size(); ++ i) {
            for (const auto & m: monomers) {
                dp[i].push_back(vector<long long>(m.seq.size()));
                for (size_t k = 0; k < m.seq.size(); ++ k) {
                    dp[i][dp[i].size() - 1][k] = INF;
                }
            }
            dp[i].push_back(vector<long long>(1));
            dp[i][monomers_num][0] = INF;
        }

        for (size_t j = 0; j < monomers.size(); ++ j) {
            Seq m = monomers[j];
            if (m.seq[0] == read.seq[0]) {
                dp[0][j][0] = match;
            } else {
                dp[0][j][0] = mismatch;
            }
            for (size_t k = 1; k < m.seq.size(); ++ k) {
                long long mm_score = monomers[j].seq[k] == read.seq[0] ? match: mismatch;
                dp[0][j][k] = max(dp[0][j][k-1] + del, (long long)(del*(k-1) + mm_score));
            }
        }
        for (size_t i = 1; i < read.seq.size(); ++ i) {
            for (size_t j = 0; j < monomers.size(); ++ j) {
                dp[i][monomers_num][0] = max(dp[i][monomers_num][0], dp[i-1][j][monomers[j].size() - 1]);
            }
            for (size_t j = 0; j < monomers.size(); ++ j) {
                for (size_t k = 0; k < monomers[j].size(); ++ k) {
                    long long score = INF;
                    int mm_score = monomers[j].seq[k] == read.seq[i] ? match: mismatch;
                    if (dp[i][monomers_num][0] > INF) {
                        score = max(score, dp[i][monomers_num][0] + mm_score + static_cast<long long>(k*del));
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
        int max_score = INF;
        int best_m = monomers_num;
        for (size_t j = 0; j < monomers.size(); ++ j) {
            if (max_score < dp[read.seq.size()-1][j][monomers[j].size() -1] ) {
                max_score = dp[read.seq.size()-1][j][monomers[j].size() -1];
                best_m = j;
            }
        }
        vector<MonomerAlignment> ans;
        long long i = read.seq.size() - 1;
        long long j = best_m;
        long long k = dp[i][j].size() - 1;
        bool monomer_changed = true;
        MonomerAlignment cur_aln;
        while (i >= 0) {
            if (k == static_cast<long long>(dp[i][j].size() - 1) && j != monomers_num && monomer_changed) {
                cur_aln = MonomerAlignment(monomers[j].read_id.name, read.read_id.name, i, i, dp[i][j][k], true);
                monomer_changed = false;
            } 
            if (j == monomers_num) {
                if (i != 0) {
                    for (size_t p = 0; p < dp[i - 1].size(); ++ p) {
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
                        int mm_score = monomers[j].seq[k] == read.seq[i] ? match: mismatch;
                        if (i != 0 && k != 0 && dp[i][j][k] == dp[i-1][j][k-1] + mm_score) {
                            --i; --k;
                        } else {
                            monomer_changed = true;
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
        int prev_end = 0;
        for (auto a: batch) {
            string s = a.read_name + "\t" 
                       + a.monomer_name + "\t" 
                       + to_string(a.start_pos) + "\t"
                       + to_string(a.end_pos) + "\t"
                       + to_string(a.identity) + "\t"
                       + to_string(a.start_pos - prev_end) + "\t"
                       + to_string(a.end_pos - a.start_pos);
            prev_end = a.end_pos;
            cout << s << "\n";
        }
    }

    vector<MonomerAlignment> PostProcessing(vector<MonomerAlignment> &batch) {
        vector<MonomerAlignment> res;
        size_t i = 0;
        while (i < batch.size()) {
            for (size_t j = i + 1; j < min(i + 7, batch.size()); ++ j) {
                if ((batch[i].end_pos - batch[j].start_pos)*2 > (batch[j].end_pos - batch[j].start_pos)) {
                    res.push_back(batch[i]);
                    i = j + 1;
                    break;
                }
            }
            if (i < batch.size() ) res.push_back(batch[i]);
            ++ i;
        }
        return res;
    }

    vector<Seq> monomers_;
    const int SAVE_STEP = 1;
    int ins_;
    int del_;
    int mismatch_;
    int match_;
};



vector<Seq> load_fasta(string filename) {
    std::ifstream input_file;
    input_file.open(filename, std::ifstream::in);
    string s;
    vector<Seq> seqs;
    while (getline(input_file, s)) {
        if (s[0] == '>') {
            string header = s.substr(1, s.size() - 1);
            istringstream iss(header);
            vector<std::string> header_v((istream_iterator<string>(iss)),
                                             istream_iterator<string>());
            seqs.push_back(Seq(header_v[0], ""));
        } else {
            seqs[seqs.size()-1].seq += s;
        }
    }
    set<char> nucs = {'A', 'C', 'G', 'T', 'N'};
    bool hasN = false;
    for (const auto & seq: seqs) {
        for (char c: seq.seq) {
            if (nucs.count(c) == 0) {
                cerr << "ERROR: Sequence " << seq.read_id.name <<" contains undefined symbol (not ACGT): " << c << endl;
                exit(-1); 
            } else if (c == 'N') {
                hasN = true;
            }
        }
    }
    if (hasN) {
        cerr << "WARNING: sequences in " << filename  << " contain N symbol. It will be counted as a separate symbol in scoring!" << endl;
    }
    return seqs;
}

string reverse_complement(string &s){
    string res = "";
    map<char, char> rc = {{'A', 'T'}, {'T', 'A'}, {'G','C'}, {'C','G'}, {'N','N'}};
    for (int i = (int) s.size() - 1; i >= 0; --i){
        try {
            res += rc.at(s[i]);
        }
        catch (std::out_of_range& e)
        {
            cerr << e.what() << std::endl;
            exit(-1);
        }
    }
    return res;
}

void add_reverse_complement(vector<Seq> &monomers) {
    vector<Seq> rev_c_monomers;
    for (auto s: monomers) {
        rev_c_monomers.push_back(Seq(s.read_id.name + "'", reverse_complement(s.seq)));
    }
    monomers.insert(monomers.end(), rev_c_monomers.begin(), rev_c_monomers.end());
    return;
}


int main(int argc, char **argv) {
    if (argc < 5) {
        cout << "Failed to process. Number of arguments < 5\n";
        cout << "./decompose <reads> <monomers> <threads> <part-size> <overlap> [<ins-score> <del-score> <mismatch-score> <match-score>]\n";
        return -1;
    }
    int ins = -1, del = -1, mismatch = -1, match = 1;
    if (argc == 10) {
        ins = stoi(argv[6]);
        del = stoi(argv[7]);
        mismatch = stoi(argv[8]);
        match = stoi(argv[9]);
    }

    int ed_thr = -1;
    if (argc == 11) {
        ed_thr = stoi(argv[10]);
    }

    cerr << "Scores: insertion=" << ins << " deletion=" << del << " mismatch=" << mismatch << " match=" << match << endl;
    vector<Seq> reads = load_fasta(argv[1]);
    vector<Seq> monomers = load_fasta(argv[2]);
    add_reverse_complement(monomers);
    MonomersAligner monomers_aligner(monomers, ins, del, mismatch, match);
    int num_threads = stoi(argv[3]);
    int part_size = stoi(argv[4]);
    int overlap = stoi(argv[5]);
    monomers_aligner.AlignReadsSet(reads, num_threads, part_size, ed_thr, overlap);
}
