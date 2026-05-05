#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>
#include <iomanip>
#include <filesystem>
#include <random>

#define FILE1 "MATa1"
#define FILE2 "MATalpha2"
#define FILE3 "MCM1"
#define FILE4 "MIG1"
#define FILE5 "PHO4"
#define FILE6 "RCS1"
#define FILE7 "ROX1"
#define FILE8 "TAF"
#define PROMOTER_FILE "promoters"
#define file "promoters"


using namespace std;
namespace fs = std::filesystem;

//バックグランド頻度を定数として定義するようにする。
const double BG_total=7519429*2+4637676*2;
const double BG[4]={
    7519429/BG_total,//A
    4637676/BG_total,//C
    4637676/BG_total,//G
    7519429/BG_total//T
};



//ファイル読み込みをしたい。
vector<string> readfile(const string&filename){
        ifstream ifs(filename);
    if (!ifs) {
        cerr << "Error: Cannot open " << filename << endl;
        exit(1);
    }
    string line;
    vector<string>seq;
        while (getline(ifs, line))//get lineで改行までの1行ずつを読み込んでくる。
        if (!line.empty()) seq.push_back(line);
    return seq;
}

//塩基を数値に変換したい。
int B_to_i(char c){
    switch(c){
        case 'A':return 0;
        case 'C':return 1;
        case 'G':return 2;
        case 'T':return 3;
        default:return -1;
    }
}

//頻度表を作りたい。
vector<vector<int>> frequency_table(const vector<string>&seq){
    int L=seq[0].size();
    int N=seq.size();
    vector<vector<int>>freq(4,vector<int>(L,0));
    for(int s=0;s<N;s++){
        const string& seqs=seq[s];
        for(int i=0;i<L;i++){
            int b=B_to_i(seqs[i]);
            if(b>=0){
                freq[b][i]++;
            }
        }
    }
    return freq;
}
//oddsスコア行列を作りたい。
vector<vector<double>> odds_score(const vector<vector<int>>&freq){
    int L=freq[0].size();
    vector<vector<double>>odds_score_table(4,vector<double>(L,0.0));
        for(int b=0;b<4;b++){
            for(int j=0;j<L;j++){
                double P=0.0;
                double p=(freq[b][j]+1.0)/(L+4.0);
                odds_score_table[b][j]=log(p/BG[b]);
            }
        }
        return odds_score_table;
    }


//promoters ファイルを読み込む。
struct Promoter{
    string gene;
    string sequence;
};
vector<Promoter> read_promoter(const string&filename){
    vector<Promoter>promoters;
    ifstream ifs(filename);
    string line;
    string gene;
    while(getline(ifs,line)){
        if(line.empty())continue;
        if(line[0]=='>'){
            gene=line.substr(1);
        }else{
            promoters.push_back({gene,line});
        }
    }
    return promoters;
}

// 閾値をランダム配列法で決定する
double calculate_asikiri(const vector<vector<double>>& odds_score_table, const double BG[4]) {
    int L = odds_score_table[0].size();
    mt19937 rng(42);
    discrete_distribution<int> dist({BG[0], BG[1], BG[2], BG[3]});
    vector<char> base_seq;
    for (int i = 0; i < 100000; i++)
        base_seq.push_back("ACGT"[dist(rng)]);

    vector<double> asikiri_scores;
    for (int j = 0; j <= (int)base_seq.size() - L; j++) {
        double score = 0.0;
        for (int k = 0; k < L; k++) {
            int b = B_to_i(base_seq[j + k]);
            if (b >= 0) score += odds_score_table[b][k];
        }
        asikiri_scores.push_back(score);
    }
    sort(asikiri_scores.begin(), asikiri_scores.end(), greater<double>());
    int idx = (int)(0.0001 * asikiri_scores.size());
    return  asikiri_scores[idx];
}

// スキャン関数
int score_scan(const vector<Promoter>& promoters,
               const vector<vector<double>>& odds_score_table,
               double asikiri) {
    int L = odds_score_table[0].size();
    for (const auto& promoter : promoters) {
        const string& gene_seq = promoter.sequence;
        for (int i = 0; i <= (int)gene_seq.size() - L; i++) {
            double s = 0.0;
            for (int j = 0; j < L; j++) {
                int b = B_to_i(gene_seq[i + j]);
                if (b >= 0) s += odds_score_table[b][j];
            }
            if (s > asikiri) {
                cout << "pro:" << promoter.gene << endl;
                cout << "pos:" << i + 1 << endl;
                cout << "hit(" << gene_seq.substr(i, L) << ")="
                     << fixed << setprecision(2) << s << endl;
                cout << endl;
            }
        }
    }
    return 0;
}


int main() {
    const char* filename[8] = {FILE1, FILE2, FILE3, FILE4,
                                FILE5, FILE6, FILE7, FILE8};
    const char* tfname[8]   = {"MATa1","MATalpha2","MCM1","MIG1",
                                "PHO4","RCS1","ROX1","TAF"};

    // ✅ 修正⑤: promotersをループの外で1回だけ読み込む
    vector<Promoter> promoters = read_promoter(PROMOTER_FILE);

    for (int i = 0; i < 8; i++) {
        vector<string> seqs = readfile(filename[i]);
        int N = seqs.size();  // 配列数
        vector<vector<int>> freq = frequency_table(seqs);

        vector<vector<double>> odds_score_table = odds_score(freq);

        cout << "Motif: " << tfname[i] << endl;
        cout << "Position\tA\tC\tG\tT" << endl;
        for (int j = 0; j < (int)odds_score_table[0].size(); j++) {
            cout << j + 1 << "\t";
            for (int b = 0; b < 4; b++)
                cout << fixed << setprecision(2) << odds_score_table[b][j] << "\t";
            cout << endl;
        }

        double asikiri = calculate_asikiri(odds_score_table, BG);
        cout << "\n閾値: " << fixed << setprecision(5) << asikiri << endl;
        cout << "Scanning promoters for " << tfname[i] << " motif..." << endl;

        score_scan(promoters, odds_score_table, asikiri);
    }
    return 0;
}


 