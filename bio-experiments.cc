// ============================================================
// 転写因子結合部位予測プログラム（閾値なし版）
// 機能: スコア行列(PWM)を用いてプロモーター配列をスキャンし、
//       各位置の一致度スコアを出力する
//
// コンパイル: g++ -std=c++17 -O2 -o tfbs_raw main.cpp
// 実行:       ./tfbs_raw
// ============================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>
#include <iomanip>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

// ============================================================
// バックグラウンド頻度（出芽酵母全ゲノムから算出）
// A:7519429, C:4637676, G:4637676, T:7519429
// ============================================================
const double BG_TOTAL = 7519429.0*2 + 4637676.0*2;
const double BG[4] = {
    7519429.0 / BG_TOTAL,  // A
    4637676.0 / BG_TOTAL,  // C
    4637676.0 / BG_TOTAL,  // G
    7519429.0 / BG_TOTAL,  // T
};

// 塩基 → インデックス (A=0, C=1, G=2, T=3)
int b2i(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return -1;
    }
}

// ============================================================
// STEP 1: motifファイルを読み込む（1行1配列）
// ============================================================
vector<string> read_motif(const string& filename) {
    vector<string> seqs;
    ifstream f(filename);
    string line;
    while (getline(f, line))
        if (!line.empty()) seqs.push_back(line);
    return seqs;
}

// ============================================================
// STEP 2: PWM（対数オッズスコア行列）を構築する
// pwm[塩基][位置] = log2( p_i(x) / q(x) )
//
// 手順:
//   1. 各列の塩基出現頻度を数える
//   2. 擬似頻度(+1補正)を加えて確率 p_i(x) を計算
//   3. バックグラウンド確率 q(x) で割り対数を取る
// ============================================================
vector<vector<double>> build_pwm(const vector<string>& seqs) {
    int L = seqs[0].size();
    int n = seqs.size();

    // 頻度表: freq[塩基][位置]
    vector<vector<int>> freq(4, vector<int>(L, 0));
    for (const string& s : seqs)
        for (int i = 0; i < L; i++) {
            int b = b2i(s[i]);
            if (b >= 0) freq[b][i]++;
        }

    // 擬似頻度(+1補正) → 確率 → 対数オッズスコア
    vector<vector<double>> pwm(4, vector<double>(L, 0.0));
    for (int b = 0; b < 4; b++)
        for (int i = 0; i < L; i++) {
            double p = (freq[b][i] + 1.0) / (n + 4.0);
            pwm[b][i] = log2(p / BG[b]);
        }
    return pwm;
}

// ============================================================
// STEP 3: 一致度スコア hit(x) を計算する
// hit(x) = Σ pwm[x_i][i]
// ============================================================
double calc_hit(const string& seq, int start,
                const vector<vector<double>>& pwm) {
    int L = pwm[0].size();
    double score = 0.0;
    for (int i = 0; i < L; i++) {
        int b = b2i(seq[start + i]);
        if (b < 0) return -1e9;
        score += pwm[b][i];
    }
    return score;
}

// ============================================================
// STEP 4: FASTAファイルを読み込む
// ============================================================
// 順序を保持するため vector<pair> で返す
vector<pair<string,string>> read_fasta(const string& filename) {
    vector<pair<string,string>> result;
    ifstream f(filename);
    string line, gene;
    while (getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            gene = line.substr(1);
            result.push_back({gene, ""});
        } else {
            result.back().second += line;
        }
    }
    return result;
}

// ============================================================
// STEP 5: プロモーター配列をスキャンして全スコアを出力する
//         （閾値なし：全位置のスコアを記録し最高スコア位置を表示）
// ============================================================
struct ScanResult {
    int    position;   // 最高スコアの位置（1始まり）
    double max_score;  // 最高スコア
    string subseq;     // 対応する部分配列
};

ScanResult scan_promoter(const string& seq,
                          const vector<vector<double>>& pwm) {
    int L = pwm[0].size();
    ScanResult best = {0, -1e9, ""};
    for (int i = 0; i <= (int)seq.size() - L; i++) {
        double score = calc_hit(seq, i, pwm);
        if (score > best.max_score) {
            best = {i + 1, score, seq.substr(i, L)};
        }
    }
    return best;
}

// ============================================================
// メイン関数
// ============================================================
int main() {
    cout << "============================================" << endl;
    cout << "  転写因子結合部位予測 (PWM法・閾値なし)" << endl;
    cout << "============================================" << endl << endl;

    // motifファイルを全て読み込む
    vector<string>                 tf_names;
    vector<vector<vector<double>>> pwm_list;

    for (const auto& entry : fs::directory_iterator("data/motif/")) {
        string tf = entry.path().filename().string();
        auto seqs = read_motif(entry.path().string());
        if (seqs.empty()) continue;
        tf_names.push_back(tf);
        pwm_list.push_back(build_pwm(seqs));
    }

    // プロモーター配列を読み込む
    auto promoters = read_fasta("data/seq/promoters");

    // ============================================================
    // 各転写因子 × 各プロモーターの最高スコア位置を出力
    // 閾値を設けず、最もスコアが高い位置を「最有力結合候補」として示す
    // ============================================================
    for (int t = 0; t < (int)tf_names.size(); t++) {
        const string& tf  = tf_names[t];
        const auto&   pwm = pwm_list[t];
        int L = pwm[0].size();

        cout << "--- 転写因子: " << tf
             << " (モチーフ長: " << L << "bp) ---" << endl;
        cout << setw(12) << left  << "遺伝子名"
             << setw(8)  << right << "位置"
             << setw(10)          << "スコア"
             << "  部分配列" << endl;
        cout << string(50, '-') << endl;

        for (const auto& [gene, seq] : promoters) {
            auto res = scan_promoter(seq, pwm);
            cout << setw(12) << left  << gene
                 << setw(8)  << right << res.position
                 << setw(10) << fixed << setprecision(2) << res.max_score
                 << "  " << res.subseq << endl;
        }
        cout << endl;
    }

    return 0;
}