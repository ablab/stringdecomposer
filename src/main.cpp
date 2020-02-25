#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map> 
#include <set>
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

class MonomerIndex{
public:
    MonomerIndex(vector<Seq> &monomers, size_t kmer_sz): monomers_(monomers), kmer_sz_(kmer_sz) {
        for (int i = 0; i < monomers_.size(); ++ i){
            if (kmer_sz_ > monomers_[i].seq.size()) continue;
            for (int j = 0; j < monomers_[i].seq.size() - kmer_sz_; ++ j) {
                string kmer = monomers_[i].seq.substr(j, kmer_sz_);
                if (index_.count(kmer) == 0) {
                    index_[kmer] = vector<pair<int, int>>();
                }
                index_[kmer].push_back(pair<int, int>(i, j));
            }
        }
        index_[""] = vector<pair<int, int>>();
        cerr << "Index size " <<  index_.size() << endl;
    }

    vector<pair<int, int>> GetKmerMonomers(string kmer){
        if (index_.count(kmer) > 0) {
            return index_[kmer];
        } else {
            return index_[""];
        }
    }

    size_t kmer_sz() {
        return kmer_sz_;
    }

private:
    vector<Seq> &monomers_;
    size_t kmer_sz_;
    unordered_map<string, vector<pair<int, int>> > index_;

};

class MonomersAligner {

public:
    MonomersAligner(vector<Seq> &monomers, int ins = -1, int del = -1, int mismatch = -1, int match = 1)
      : monomers_(monomers),
        ins_(ins),
        del_(del),
        mismatch_(mismatch),
        match_(match),
        index_(monomers_, 10){}

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
            if (save_steps.size() > 1) {
                save_steps[save_steps.size() - 1] += save_steps[save_steps.size() - 2]; 
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

        vector<pair<int, vector<MonomerAlignment>>> subbatches;
        #pragma omp parallel for num_threads(threads)
        for (int j = 0; j < new_reads.size(); ++ j) {
            clock_t begin = clock();
            vector<MonomerAlignment> aln = AlignPartClassicDP(new_reads[j]);
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
        for (int p = 0; p < save_steps.size(); ++ p){
            vector<MonomerAlignment> batch;
            int start = p == 0? 0: save_steps[p-1];
            for (int j = start; j < save_steps[p]; ++ j) {
                int read_index = subbatches[j].first;
                for (auto a: subbatches[j].second) {
                    MonomerAlignment new_m_aln(a.monomer_name, a.read_name, 
                                                new_reads[read_index].read_id.id + a.start_pos, new_reads[read_index].read_id.id + a.end_pos, 
                                                a.identity, a.best);
                    batch.push_back(new_m_aln);
                }
            }
            batch = PostProcessing(batch);
            cerr << "Aligned " << batch[0].read_name << endl;
            #pragma omp critical(aligner2)
            {
                SaveBatch(batch);
            }
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
        vector<vector<vector<int>>> &dp = dp_[omp_get_thread_num()];
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


        map<int, int> read_kmers;
        string s = read.seq.substr(0, 171);
        for (int j = 0; j < s.size() - index_.kmer_sz(); ++ j) {
            string kmer = s.substr(j, index_.kmer_sz());
            vector<pair<int, int>> kmer_mono = index_.GetKmerMonomers(kmer);
            for (auto p: kmer_mono) {
                if (abs(p.second - j) < 10) {
                    if (read_kmers.count(p.first) == 0) {
                        read_kmers[p.first] = 0;
                    }
                    read_kmers[p.first] += 1;
                }
            }
        }
        //vector<int> similar_monomers = index_.GetSimilarMonomers(read.seq, 0);
        //cout << read_kmers.size() << endl;
        unordered_map<int, int> valid_monomers;
        for (auto const &it: read_kmers){
            if (it.second > 20) {
                valid_monomers[it.first] = 0;
            }
        }
        if (valid_monomers.size() == 0){
            for (int j = 0; j < monomers_num; ++ j){
                valid_monomers[j] = 0;
            }
        }
        //cout << "i=0 " <<  valid_monomers.size() << endl;
        for (auto &it: valid_monomers){ 
            Seq m = monomers_[it.first];
            int j = it.first;
            if (m.seq[0] == read.seq[0]) {
                dp[0][j][0] = match;
            } else {
                dp[0][j][0] = mismatch;
            }
            for (int k = 1; k < m.seq.size(); ++ k) {
                int mm_score = m.seq[k] == read.seq[0] ? match: mismatch;
                dp[0][j][k] = max(dp[0][j][k-1] + del, del*(k-1) + mm_score);
            }
            ++ it.second;
        }
        for (int i = 1; i < read.seq.size(); ++ i) {
            //cout << "i=" << i << endl;
            read_kmers.clear();
            string s = read.seq.substr(i, 171);
            if (s.size() > index_.kmer_sz()){
                for (int j = 0; j < s.size() - index_.kmer_sz(); ++ j) {
                    string kmer = s.substr(j, index_.kmer_sz());
                    vector<pair<int, int>> kmer_mono = index_.GetKmerMonomers(kmer);
                    for (auto p: kmer_mono) {
                        if (abs(p.second - j) < 20) {
                            if (read_kmers.count(p.first) == 0) {
                                read_kmers[p.first] = 0;
                            }
                            read_kmers[p.first] += 1;
                        }
                    }
                }
            }
            //cout << "read kmers size " << read_kmers.size() << endl;
            set<int> invalid_kmers;
            for (auto const &it: read_kmers){
                if (it.second > 20) {
                    valid_monomers[it.first] = 0;
                } else{
                    if (it.second == 0) {
                        invalid_kmers.insert(it.first);
                    }
                }
            }
            //cout << "valid monomers " << valid_monomers.size() << endl;
            if (valid_monomers.size() == 0){
                for (int j = 0; j < monomers_num; ++ j){
                    valid_monomers[j] = 0;
                }
            }
            for (auto const &it: invalid_kmers) {
                read_kmers.erase(it);
            }

            for (int j = 0; j < monomers_.size(); ++ j) {
                dp[i][monomers_num][0] = max(dp[i][monomers_num][0], dp[i-1][j][monomers_[j].size() - 1]);
            }
            //for (int j = 0; j < monomers_.size(); ++ j) {
            //    for (int k = 0; k < monomers_[j].size(); ++ k) {
            set<int> invalid;
            //cout << "valid monomers " << valid_monomers.size() << endl;
            for (auto &it: valid_monomers){
                int j = it.first;
                //cout << " " << monomers_[j].read_id.name << endl;
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
                ++ it.second;
                if (it.second > 200) {
                    invalid.insert(it.first);
                }
            }
            for (auto const &it: invalid) {
                valid_monomers.erase(it);
            }
        }
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

private:
    vector<Seq> monomers_;
    const int SAVE_STEP = 1;
    int ins_;
    int del_;
    int mismatch_;
    int match_;
    MonomerIndex index_;
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
