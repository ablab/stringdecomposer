#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cstring>
#include <omp.h>
#include <ctime>

#include "edlib.h"

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

    Seq(string name_, string seq_): seq(seq_), read_id(ReadId(name_)) {
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    }

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
    MonomersAligner(vector<Seq> &monomers, int ins = -1, int del = -1, int mismatch = -1, int match = 1)
      : monomers_(monomers),
        ins_(ins),
        del_(del),
        mismatch_(mismatch),
        match_(match) {
    }

    void AlignReadsSet(vector<Seq> &reads, int threads, int part_size) {
        vector<Seq> new_reads;
        vector<int> save_steps;
        for (auto r: reads) {
            int cnt = 0;
            //cout << r.seq.size() << endl;
            for (int i = 0; i < r.seq.size(); i += part_size) {
                if ((int) r.seq.size() - i >= 500 || r.seq.size() < 500) {
                    Seq seq = Seq(r.read_id.name, r.seq.substr(i, min(part_size + 500, (int) r.seq.size() - i) ), i );
                    new_reads.push_back(seq);
                    ++ cnt;
                }
            }
        }
        cerr << reads.size() << " reads divided into " << new_reads.size() << " parts\n";

        dp_ = vector<vector<vector<vector<int>>>>(threads, vector<vector<vector<int>>>(part_size + 500, vector<vector<int>>(monomers_.size() + 1)));
        for (int k = 0; k < threads; ++ k) {
            for (int i = 0; i < part_size + 500; ++ i) {
                for (int j = 0; j < monomers_.size(); ++ j) {
                    dp_[k][i][j] = vector<int>(monomers_[j].seq.size());
                }
                dp_[k][i][monomers_.size()] = vector<int>(1);
            }
        }
        cerr << "Prepared memory\n";
        cerr << threads << endl;

        int step = threads*2;
        for (int i = 0; i < new_reads.size(); i += step) {
            vector<pair<int, vector<MonomerAlignment>>> subbatches;
            #pragma omp parallel for num_threads(threads)
            for (int j = i; j < min(i + step, (int) new_reads.size()); ++ j) {
                clock_t begin = clock();
                vector<int> cur_monomers = ChooseBestMonomers(new_reads[j]);
                vector<MonomerAlignment> aln = AlignPartClassicDP(new_reads[j], cur_monomers);
                clock_t end = clock();
                double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                cerr << "Thread " << omp_get_thread_num() << " Time " <<  elapsed_secs << endl;
                #pragma omp critical(aligner)
                {
                    subbatches.push_back(pair<int, vector<MonomerAlignment>> (j, aln));
                }
            }
            sort(subbatches.begin(), subbatches.end(), sortby1);

            #pragma omp parallel for num_threads(threads)
            for (int p = i; p < min(i + step, (int) new_reads.size()); ++ p){
                vector<MonomerAlignment> batch;
                int read_index = subbatches[p - i].first;
                for (auto a: subbatches[p - i].second) {
                    MonomerAlignment new_m_aln(a.monomer_name, a.read_name, 
                                                new_reads[read_index].read_id.id + a.start_pos, new_reads[read_index].read_id.id + a.end_pos, 
                                                a.identity, a.best);
                    batch.push_back(new_m_aln);
                }
                batch = PostProcessing(batch);
                #pragma omp critical(aligner2)
                {
                    SaveBatch(batch);
                    cerr << "Aligned " << batch[0].read_name << endl;
                }
            }
        }         
    }

    ~MonomersAligner() {
    }

private:

    int EditDistance(string &query, string &target) {
        EdlibAlignResult result = edlibAlign(query.c_str(), query.size(), target.c_str(), target.size(), 
                                            edlibNewAlignConfig(query.size()/4, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
        int editdist = query.size() + 100;
        if (result.status == EDLIB_STATUS_OK && result.editDistance >= 0) {
            editdist = result.editDistance;
        }
        edlibFreeAlignResult(result);
        return editdist;
    }

    vector<int> ChooseBestMonomers(Seq &read) {
        vector<pair<int, int>> best_monomers(monomers_.size());
        for (int i = 0; i < monomers_.size(); ++ i) {
            int edist = EditDistance(monomers_[i].seq, read.seq);
            best_monomers[i].first = edist;
            best_monomers[i].second = i;
        }
        sort(best_monomers.begin(), best_monomers.end());
        int cnt = min(100, monomers_.size());
        vector<int> res(cnt);
        for (int i = 0; i < cnt; ++ i) {
            res[i] = best_monomers[i].second;
            //cout << best_monomers[i].first << " " << monomers_[best_monomers[i].second].read_id.name << endl;
        }
        return res;
    }

    void PrecalculateMonomerAlignment(){

    }

    void PrecalculateMonomerEdlib(){

    }

    vector<MonomerAlignment> AlignPartFitting(Seq &read) {

    }

    vector<MonomerAlignment> AlignPartClassicDP(Seq &read, vector<int> &cur_monomers) {
        int ins = ins_;
        int del = del_;
        int match = match_;
        int mismatch = mismatch_;
        int INF = -1000000;
        int monomers_num = (int) monomers_.size();
        vector<vector<vector<int>>> &dp = dp_[omp_get_thread_num()];
        for (int i = 0; i < read.seq.size(); ++ i) {
            for (int p = 0; p < cur_monomers.size(); ++ p) {
                int j = cur_monomers[p];
                for (int k = 0; k < monomers_[j].seq.size(); ++ k) {
                    dp[i][j][k] = INF;
                }
            }
            dp[i][monomers_num][0] = INF;
        }

        for (int p = 0; p < cur_monomers.size(); ++ p) {
            int j = cur_monomers[p];
            Seq m = monomers_[j];
            if (m.seq[0] == read.seq[0]) {
                dp[0][j][0] = match;
            } else {
                dp[0][j][0] = mismatch;
            }
            for (int k = 1; k < m.seq.size(); ++ k) {
                int mm_score = monomers_[j].seq[k] == read.seq[0] ? match: mismatch;
                dp[0][j][k] = max(dp[0][j][k-1] + del, del*(k-1) + mm_score);
            }
        }
        for (int i = 1; i < read.seq.size(); ++ i) {
            for (int p = 0; p < cur_monomers.size(); ++ p) {
                int j = cur_monomers[p];
                dp[i][monomers_num][0] = max(dp[i][monomers_num][0], dp[i-1][j][monomers_[j].size() - 1]);
            }
            for (int p = 0; p < cur_monomers.size(); ++ p) {
                int j = cur_monomers[p];
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
        int max_score = INF;
        int best_m = monomers_num;
        for (int p = 0; p < cur_monomers.size(); ++ p) {
            int j = cur_monomers[p];
            if (max_score < dp[read.seq.size()-1][j][monomers_[j].size() -1] ) {
                max_score = dp[read.seq.size()-1][j][monomers_[j].size() -1];
                best_m = j;
            }
        }
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
                    for (int pp = 0; pp < cur_monomers.size(); ++ pp) {
                        int p = cur_monomers[pp];
                    //for (int p = 0; p < dp[i - 1].size(); ++ p) {
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
            cout << s << "\n";
        }
    }

    vector<MonomerAlignment> PostProcessing(vector<MonomerAlignment> &batch) {
        vector<MonomerAlignment> res;
        res.push_back(batch[0]);
        for (size_t i = 1; i < batch.size(); ++ i) {
            bool add = true;
            for (size_t j = (size_t) max((int) 0, (int) i - 6); j < i; ++ j) {
                if ((abs(batch[i].end_pos - batch[j].end_pos) < 50 || abs(batch[i].start_pos - batch[j].start_pos) < 50) && batch[i].identity < batch[j].identity) {
                    add = false;
                }
            }
            for (size_t j = i + 1; j < (size_t) min((int) i + 7, (int) batch.size()); ++ j) {
                if ((abs(batch[j].end_pos - batch[i].end_pos) < 50 || abs(batch[j].start_pos - batch[i].start_pos) < 50) && batch[i].identity <= batch[j].identity) {
                    add = false;
                }
            }
            if (add) { res.push_back(batch[i]);}
        }
        return res;
    }

    vector<Seq> monomers_;
    const int SAVE_STEP = 1;
    int ins_;
    int del_;
    int mismatch_;
    int match_;
    vector<vector<vector<vector<int>>>> dp_;
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
        rev_c_monomers.push_back(Seq(s.read_id.name + "'", reverse_complement(s.seq)));
    }
    monomers.insert(monomers.end(), rev_c_monomers.begin(), rev_c_monomers.end());
    return;
}


int main(int argc, char **argv) {
    if (argc < 3) {
        cout << "Failed to process. Number of arguments < 4\n";
        cout << "./decompose <reads> <monomers> <threads> [<ins-score> <del-score> <mismatch-score> <match-score>]\n";
        return -1;
    }
    int ins = -1, del = -1, mismatch = -1, match = 1;
    if (argc == 8) {
        ins = stoi(argv[4]);
        del = stoi(argv[5]);
        mismatch = stoi(argv[6]);
        match = stoi(argv[7]);
    }
    cerr << "Scores: insertion=" << ins << " deletion=" << del << " mismatch=" << mismatch << " match=" << match << endl;
    vector<Seq> reads = load_fasta(argv[1]);
    vector<Seq> monomers = load_fasta(argv[2]);
    add_reverse_complement(monomers);
    MonomersAligner monomers_aligner(monomers, ins, del, mismatch, match);
    int num_threads = stoi(argv[3]);
    int part_size = 20000;
    monomers_aligner.AlignReadsSet(reads, num_threads, part_size);
}
