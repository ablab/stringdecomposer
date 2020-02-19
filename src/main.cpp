#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cstring>
#include <omp.h>
#include <ctime>

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
            save_steps.push_back(cnt);
        }
        cerr << "Prepared reads\n";
        dp = vector<vector<vector<int>>>(part_size + 500, vector<vector<int>>(monomers_.size() + 1));
        for (int i = 0; i < part_size + 500; ++ i) {
            for (int j = 0; j < monomers_.size(); ++ j) {
                dp[i][j] = vector<int>(monomers_[j].seq.size());
            }
            dp[i][monomers_.size()] = vector<int>(1);
        }
        dp_local = vector<vector<vector<pair<int, int> >>>(part_size + 500, vector<vector<pair<int, int>>>(monomers_.size() + 1));
        for (int i = 0; i < part_size + 500; ++ i) {
            for (int j = 0; j < monomers_.size(); ++ j) {
                dp_local[i][j] = vector<pair<int, int>>(monomers_[j].seq.size());
            }
            dp_local[i][monomers_.size()] = vector<pair<int, int>>(1);
        }
        int start = 0;
        for (int p = 0; p < save_steps.size(); ++ p){
            vector<MonomerAlignment> batch;
            vector<pair<int, vector<MonomerAlignment>>> subbatches;
            #pragma omp parallel for num_threads(threads)
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
            cerr << "Aligned " << batch[0].read_name << endl;
            batch = PostProcessing(batch);
            SaveBatch(batch);
            start += save_steps[p];
        }           
    }

    ~MonomersAligner() {
    }

private:

    void PrecalculateMonomerAlignment(){

    }

    void PrecalculateMonomerEdlib(){

    }

    vector<MonomerAlignment> AlignPartFitting(Seq &read) {

    }

    vector<MonomerAlignment> AlignPartClassicDP(Seq &read) {
        int ins = ins_;
        int del = del_;
        int match = match_;
        int mismatch = mismatch_;
        int INF = -1000000;
        int monomers_num = (int) monomers_.size();
        vector<vector<int>> queue(monomers_num);
        for (int i = 0; i < read.seq.size(); ++ i) {
            dp[i][monomers_num][0] = 0;
        }
        for (int i = 0; i < monomers_num; ++ i) {
            queue[i].resize(monomers_[i].seq.size() + 1);
            queue[i][0] = 1;
        }
        float per = 0.85;
        int len = 25;
        vector<float> per_read(read.seq.size());
        vector<int> len_read(read.seq.size());
        for (int i = 0; i < min((int) read.seq.size(), 300); ++ i) {
            per_read[i] = per;
            len_read[i] = 200;
            per_read[read.seq.size() - i - 1] = per;
            len_read[read.seq.size() - i - 1] = 200;
        }
        for (int i = 300; i < max((int)read.seq.size() - 300, 0); ++ i) {
            per_read[i] = per;
            len_read[i] = len;
        }
        for (int j = 0; j < monomers_.size(); ++ j) {
            Seq m = monomers_[j];
            if (m.seq[0] == read.seq[0]) {
                dp[0][j][0] = match;
                dp_local[0][j][0].first = 1;
                dp_local[0][j][0].second = 1;
            } else {
                dp[0][j][0] = mismatch;
                dp_local[0][j][0].first = 0;
                dp_local[0][j][0].second = 1;
            }
            queue[j][0] = 2;
            queue[j][1] = 0;
            for (int k = 1; k < m.seq.size(); ++ k) {
                int mm_score = monomers_[j].seq[k] == read.seq[0] ? match: mismatch;
                if (k < len_read[0]) {
                    dp[0][j][k] = max(dp[0][j][k-1] + del, del*(k-1) + mm_score);
                    if (dp[0][j][k] == dp[0][j][k-1] + del) {
                        dp_local[0][j][k].first = dp_local[0][j][k-1].first;
                        dp_local[0][j][k].second = k;       
                    } else {
                        dp_local[0][j][k].first = 1;
                        dp_local[0][j][k].second = k;
                    }
                    queue[j][queue[j][0]] = k;
                    queue[j][0] ++;
                } 
            }
        }
        vector<vector<int>> new_queue(monomers_num);
        for (int i = 0; i < monomers_num; ++ i) {
            new_queue[i].resize(monomers_[i].seq.size() + 1);
        }
        clock_t begin = clock();
        for (int i = 1; i < read.seq.size(); ++ i) {
            for (int j = 0; j < monomers_.size(); ++ j) {
                dp[i][monomers_num][0] = max(dp[i][monomers_num][0], dp[i-1][j][monomers_[j].size() - 1]);
            }
            int cnt = 0;
            for (int j = 0; j < monomers_num; ++ j) {
                new_queue[j][0] = 1;
                int k = 0;
                while (k < len && dp[i][monomers_num][0] > INF) {
                    int mm_score = monomers_[j].seq[0] == read.seq[i] ? match: mismatch;
                    dp[i][j][k] = dp[i][monomers_num][0] + mm_score + k*del;
                    dp_local[i][j][k].first = mm_score == match ? 1: 0;
                    dp_local[i][j][k].second = k; 
                    new_queue[j][new_queue[j][0]] = k;
                    new_queue[j][0] ++;
                    ++ k;
                }
                int q2 = 1;

                for (int q1 = 1; q1 < queue[j][0]; ++ q1) {
                    while (q2 < new_queue[j][0] && new_queue[j][q2] < queue[j][q1] && monomers_[j].size() > new_queue[j][q2] + 1) {
                        k = new_queue[j][q2];
                        if ((dp_local[i][j][k].second + 1) < len_read[i] || dp_local[i][j][k].first/float(dp_local[i][j][k].second+1) > per_read[i]){
                            if (new_queue[j][0] == q2 + 1){
                                new_queue[j][new_queue[j][0]] = k + 1;
                                ++ new_queue[j][0];
                                dp[i][j][k + 1] = dp[i][j][k] + del;
                                dp_local[i][j][k + 1].first = dp_local[i][j][k].first;
                                dp_local[i][j][k + 1].second = dp_local[i][j][k].second + 1;
                            } else {
                                dp[i][j][k + 1] = max(dp[i][j][k + 1], dp[i][j][k] + del);
                                if (dp[i][j][k + 1] == dp[i][j][k] + del) {
                                    dp_local[i][j][k + 1].first = dp_local[i][j][k].first;
                                    dp_local[i][j][k + 1].second = dp_local[i][j][k].second + 1;
                                }
                            }
                        }
                        ++ q2;
                    }
                    if (q2 == new_queue[j][0]) {
                        k = queue[j][q1];
                        if ((dp_local[i-1][j][k].second + 1) < len_read[i] || dp_local[i-1][j][k].first/float(dp_local[i-1][j][k].second+1) > per_read[i]){
                            new_queue[j][new_queue[j][0]] = k;
                            ++ new_queue[j][0];
                            dp[i][j][k] = dp[i-1][j][k] + ins;
                            dp_local[i][j][k].first = dp_local[i-1][j][k].first;
                            dp_local[i][j][k].second = dp_local[i-1][j][k].second + 1; 
                        }
                    } else {
                        k = new_queue[j][q2];
                        if ((dp_local[i-1][j][k].second + 1) < len_read[i] || dp_local[i-1][j][k].first/float(dp_local[i-1][j][k].second+1) > per_read[i]){
                            dp[i][j][k] = max(dp[i-1][j][k] + ins, dp[i][j][k]);
                            if (dp[i][j][k] == dp[i-1][j][k] + ins) {
                                dp_local[i][j][k].first = dp_local[i-1][j][k].first;
                                dp_local[i][j][k].second = dp_local[i-1][j][k].second + 1;
                            }
                        }
                    }
                    if (monomers_[j].size() > queue[j][q1] + 1) {
                        k = queue[j][q1];
                        int mm_score = monomers_[j].seq[k + 1] == read.seq[i] ? match: mismatch; 
                        int local_mm_score = monomers_[j].seq[k + 1] == read.seq[i] ? 1: 0; 
                        if ((dp_local[i][j][k].second + 1) < len_read[i] || (dp_local[i][j][k].first + local_mm_score)/float(dp_local[i][j][k].second+1) > per_read[i]){
                            if (new_queue[j][0] == q2 + 1){
                                new_queue[j][new_queue[j][0]] = k + 1;
                                ++ new_queue[j][0];
                                dp[i][j][k + 1] = dp[i-1][j][k] + mm_score;
                                dp_local[i][j][k + 1].first = dp_local[i-1][j][k].first + local_mm_score;
                                dp_local[i][j][k + 1].second = dp_local[i-1][j][k].second + 1;
                            } else {
                                dp[i][j][k + 1] = max(dp[i][j][k + 1], dp[i-1][j][k] + mm_score);
                                if (dp[i][j][k + 1] == dp[i-1][j][k] + mm_score) {
                                    dp_local[i][j][k + 1].first = dp_local[i-1][j][k].first + local_mm_score;
                                    dp_local[i][j][k + 1].second = dp_local[i-1][j][k].second + 1;
                                }
                            }
                        }
                    }

                }
                while (q2 < new_queue[j][0] && monomers_[j].size() > new_queue[j][q2] + 1) {
                    k = new_queue[j][q2];
                    if ((dp_local[i][j][k].second + 1) < len_read[i] || dp_local[i][j][k].first/float(dp_local[i][j][k].second+1) > per_read[i]){
                        if (new_queue[j][0] == q2 + 1){
                            new_queue[j][new_queue[j][0]] = k + 1;
                            ++ new_queue[j][0];
                            dp[i][j][k + 1] = dp[i][j][k] + del;
                            dp_local[i][j][k + 1].first = dp_local[i][j][k].first;
                            dp_local[i][j][k + 1].second = dp_local[i][j][k].second + 1;
                        } else {
                            dp[i][j][k + 1] = max(dp[i][j][k + 1], dp[i][j][k] + del);
                            if (dp[i][j][k + 1] == dp[i][j][k] + del) {
                                dp_local[i][j][k + 1].first = dp_local[i][j][k].first;
                                dp_local[i][j][k + 1].second = dp_local[i][j][k].second + 1;
                            }
                        }
                    }
                    ++ q2;
                }
                cnt += new_queue[j][0];
            }
            queue.swap(new_queue);
        }

        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cerr << "Time " <<  elapsed_secs << endl;
        int max_score = INF;
        int best_m = monomers_num;
        for (int j = 0; j < monomers_.size(); ++ j) {
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
    vector<vector<vector<int>>> dp;
    vector<vector<vector<pair<int, int> >>> dp_local;
    vector<vector<int>> queue;
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
    if (argc < 4) {
        cout << "Failed to process. Number of arguments < 5\n";
        cout << "./decompose <reads> <monomers> <threads> <part-size> [<ins-score> <del-score> <mismatch-score> <match-score>]\n";
        return -1;
    }
    int ins = -1, del = -1, mismatch = -1, match = 1;
    if (argc == 9) {
        ins = stoi(argv[5]);
        del = stoi(argv[6]);
        mismatch = stoi(argv[7]);
        match = stoi(argv[8]);
    }
    cerr << "Scores: insertion=" << ins << " deletion=" << del << " mismatch=" << mismatch << " match=" << match << endl;
    vector<Seq> reads = load_fasta(argv[1]);
    vector<Seq> monomers = load_fasta(argv[2]);
    add_reverse_complement(monomers);
    MonomersAligner monomers_aligner(monomers, ins, del, mismatch, match);
    int num_threads = stoi(argv[3]);
    int part_size = stoi(argv[4]);
    monomers_aligner.AlignReadsSet(reads, num_threads, part_size);
}
