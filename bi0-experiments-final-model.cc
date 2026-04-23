#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

// ===== マクロ定義 =====
#define NUM_FEATURES 53
#define NUM_SEQS     10000

// ===== 決定木ノード構造体 =====
struct TreeNode {
    int    feature_id;    // 分岐に使う特徴量のインデックス
    double threshold;     // 分岐閾値（以下→左、超える→右）
    int    left_class_id; // 左側の予測ラベル
    int    right_class_id;// 右側の予測ラベル
};

// ===== (1) ファイル読み込み =====
// protein_solubility_dataset.txt を読み込み、
// feature_name・dataset・labels に格納する
void LoadSolubilityFile(
    const string& filename,
    vector<string>&            feature_name,
    vector<vector<double>>&    dataset,
    vector<int>&                    labels)
{
    ifstream ifs(filename);
    if (!ifs) {
        cerr << "Error: Cannot open file: " << filename << endl;
        return;
    }

    string line;

    // --- 1行目：ヘッダー（id + 53特徴量名 + label） ---
    getline(ifs, line);
    {
        istringstream ss(line);
        string token;
        ss >> token; // "id" を読み飛ばし
        for (int i = 0; i < NUM_FEATURES; i++) {
            ss >> token;
            feature_name[i] = token;
        }
        // "label" は読み飛ばし
    }

    // --- 2行目以降：データ行 ---
    int idx = 0;
    while (getline(ifs, line) && idx < NUM_SEQS) {
        istringstream ss(line);
        string token;
        ss >> token; // protein id を読み飛ばし

        for (int j = 0; j < NUM_FEATURES; j++) {
            ss >> dataset[idx][j];
        }
        int lbl;
        ss >> lbl;
        labels[idx] = lbl;
        idx++;
    }
}

// ===== (2) データセット分割 =====
// インデックス配列をシャッフルし、先頭 test_ratio 割をテスト用に割り当てる
void DivideDataset(
    const vector<vector<double>>& dataset,
    const vector<int>&                 labels,
    vector<vector<double>>&       training_dataset,
    vector<int>&                       training_labels,
    vector<vector<double>>&       test_dataset,
    vector<int>&                       test_labels,
    double test_ratio)
{
    int n = (int)dataset.size();
    int n_test  = (int)(n * test_ratio);
    int n_train = n - n_test;

    // インデックス配列を作成してシャッフル
    vector<int> idx(n);
    iota(idx.begin(), idx.end(), 0);

    mt19937 rng(42); // seed 固定（再現性確保）
    shuffle(idx.begin(), idx.end(), rng);

    // テストデータ
    test_dataset.resize(n_test);
    test_labels.resize(n_test);
    for (int i = 0; i < n_test; i++) {
        test_dataset[i] = dataset[idx[i]];
        test_labels[i]  = labels[idx[i]];
    }

    // トレーニングデータ
    training_dataset.resize(n_train);
    training_labels.resize(n_train);
    for (int i = 0; i < n_train; i++) {
        training_dataset[i] = dataset[idx[n_test + i]];
        training_labels[i]  = labels[idx[n_test + i]];
    }
}

// ===== (3) ジニ不純度の計算 =====
// ラベル列から Gini = 2p(1-p) を返す
double GiniImpurity(const vector<int>& labs) {
    if (labs.empty()) return 0.0;
    int pos = 0;
    for (int l : labs) pos += l;
    double p = (double)pos / labs.size();
    return 2.0 * p * (1.0 - p);
}

// ===== (4) 重み付きジニ不純度の計算 =====
// 左右のラベル列から (N_L/N)*G_L + (N_R/N)*G_R を返す
double WeightedGini(
    const vector<int>& left_labs,
    const vector<int>& right_labs)
{
    int nl = (int)left_labs.size();
    int nr = (int)right_labs.size();
    int n  = nl + nr;
    if (n == 0) return 0.0;
    return ((double)nl / n) * GiniImpurity(left_labs)
         + ((double)nr / n) * GiniImpurity(right_labs);
}

// ===== (5) 多数決ラベルの取得 =====
int MajorityLabel(const vector<int>& labs) {
    int pos = 0;
    for (int l : labs) pos += l;
    return (pos * 2 >= (int)labs.size()) ? 1 : 0;
}

// ===== (6) 1ノードの学習（深さ1用） =====
// 全特徴量 × 1〜99パーセンタイル閾値を全探索し、
// 重み付きジニ不純度を最小化する分割を TreeNode に書き込む
void TrainDecisionNode(
    const vector<vector<double>>& data,
    const vector<int>&                 labs,
    TreeNode&                               node)
{
    int n = (int)data.size();
    double best_gini = 1e18;
    int    best_feat = 0;
    double best_thr  = 0.0;
    int    best_left = 0, best_right = 1;

    for (int f = 0; f < NUM_FEATURES; f++) {
        // 特徴量 f の値をソートして 1〜99 パーセンタイルを閾値候補にする
        vector<double> vals(n);
        for (int i = 0; i < n; i++) vals[i] = data[i][f];
        sort(vals.begin(), vals.end());

        for (int p = 1; p <= 99; p++) {
            int pos_idx = (int)((double)p / 100.0 * (n - 1));
            double thr = vals[pos_idx];

            vector<int> left_labs, right_labs;
            for (int i = 0; i < n; i++) {
                if (data[i][f] <= thr) left_labs.push_back(labs[i]);
                else                   right_labs.push_back(labs[i]);
            }

            double g = WeightedGini(left_labs, right_labs);
            if (g < best_gini) {
                best_gini  = g;
                best_feat  = f;
                best_thr   = thr;
                best_left  = MajorityLabel(left_labs);
                best_right = MajorityLabel(right_labs);
            }
        }
    }

    node.feature_id     = best_feat;
    node.threshold      = best_thr;
    node.left_class_id  = best_left;
    node.right_class_id = best_right;
}

// ===== (7) 性能評価（深さ1用） =====
void Evaluation(
    const TreeNode&                         node,
    const vector<vector<double>>& test_data,
    const vector<int>&                 test_labs)
{
    int tp = 0, fp = 0, fn = 0, tn = 0;
    for (int i = 0; i < (int)test_data.size(); i++) {
        int pred = (test_data[i][node.feature_id] <= node.threshold)
                   ? node.left_class_id
                   : node.right_class_id;
        int actual = test_labs[i];
        if      (pred == 1 && actual == 1) tp++;
        else if (pred == 1 && actual == 0) fp++;
        else if (pred == 0 && actual == 1) fn++;
        else                               tn++;
    }

    double accuracy  = (double)(tp + tn) / (tp + fp + fn + tn);
    double precision = (tp + fp > 0) ? (double)tp / (tp + fp) : 0.0;
    double recall    = (tp + fn > 0) ? (double)tp / (tp + fn) : 0.0;
    double f_score   = (precision + recall > 0)
                       ? 2.0 * precision * recall / (precision + recall) : 0.0;

    cout << fixed << setprecision(4);
    cout << "Accuracy:  " << accuracy  << "\n"
              << "Precision: " << precision << "\n"
              << "Recall:    " << recall    << "\n"
              << "F-score:   " << f_score   << "\n"
              << "Confusion Matrix\n"
              << "TP: " << tp << "  FP: " << fp << "\n"
              << "FN: " << fn << "  TN: " << tn << "\n";
}

// ===== (8) 性能評価（深さ2用） =====
// decision_tree[0] が親ノード、[1] が左子、[2] が右子
void Evaluation(
    const vector<TreeNode>&            dtree,
    const vector<vector<double>>& test_data,
    const vector<int>&                 test_labs)
{
    int tp = 0, fp = 0, fn = 0, tn = 0;
    for (int i = 0; i < (int)test_data.size(); i++) {
        const TreeNode& root = dtree[0];
        int pred;
        if (test_data[i][root.feature_id] <= root.threshold) {
            // 左子ノード
            const TreeNode& child = dtree[1];
            pred = (test_data[i][child.feature_id] <= child.threshold)
                   ? child.left_class_id : child.right_class_id;
        } else {
            // 右子ノード
            const TreeNode& child = dtree[2];
            pred = (test_data[i][child.feature_id] <= child.threshold)
                   ? child.left_class_id : child.right_class_id;
        }

        int actual = test_labs[i];
        if      (pred == 1 && actual == 1) tp++;
        else if (pred == 1 && actual == 0) fp++;
        else if (pred == 0 && actual == 1) fn++;
        else                               tn++;
    }

    double accuracy  = (double)(tp + tn) / (tp + fp + fn + tn);
    double precision = (tp + fp > 0) ? (double)tp / (tp + fp) : 0.0;
    double recall    = (tp + fn > 0) ? (double)tp / (tp + fn) : 0.0;
    double f_score   = (precision + recall > 0)
                       ? 2.0 * precision * recall / (precision + recall) : 0.0;

    cout << fixed << setprecision(4);
    cout << "Accuracy:  " << accuracy  << "\n"
              << "Precision: " << precision << "\n"
              << "Recall:    " << recall    << "\n"
              << "F-score:   " << f_score   << "\n"
              << "Confusion Matrix\n"
              << "TP: " << tp << "  FP: " << fp << "\n"
              << "FN: " << fn << "  TN: " << tn << "\n";
}

// ===== (9) 深さ2の決定木の学習 =====
void TrainDecisionTree(
    const vector<vector<double>>& data,
    const vector<int>&                 labs,
    vector<TreeNode>&                  dtree)
{
    // --- 親ノード(dtree[0])を学習 ---
    TrainDecisionNode(data, labs, dtree[0]);

    // --- 親ノードでデータを左右に分割 ---
    vector<vector<double>> left_data,  right_data;
    vector<int>                 left_labs,  right_labs;

    for (int i = 0; i < (int)data.size(); i++) {
        if (data[i][dtree[0].feature_id] <= dtree[0].threshold) {
            left_data.push_back(data[i]);
            left_labs.push_back(labs[i]);
        } else {
            right_data.push_back(data[i]);
            right_labs.push_back(labs[i]);
        }
    }

    // --- 左子ノード(dtree[1])を学習 ---
    TrainDecisionNode(left_data, left_labs, dtree[1]);

    // --- 右子ノード(dtree[2])を学習 ---
    TrainDecisionNode(right_data, right_labs, dtree[2]);
}

// ===== main 関数 =====
int main(void) {
    // データ構造の初期化
    vector<string>         feature_name(NUM_FEATURES, "");
    vector<vector<double>> dataset(NUM_SEQS,
                                     vector<double>(NUM_FEATURES, 0.0));
    vector<int> labels(NUM_SEQS);

    // (1) ファイル読み込み
    LoadSolubilityFile("protein_solubility_dataset.txt",
                       feature_name, dataset, labels);

    // (2) データ分割（テスト 20%）
    vector<vector<double>> training_dataset, test_dataset;
    vector<int>                 training_labels,  test_labels;
    double test_ratio = 0.2;

    DivideDataset(dataset, labels,
                  training_dataset, training_labels,
                  test_dataset,     test_labels,
                  test_ratio);

    cout << "Train size: " << training_dataset.size()
              << "  Test size: " << test_dataset.size() << "\n\n";

    // ===== 実装(2)：ベースライン（全員可溶と予測）の評価 =====
    cout << "=== Baseline (all predict soluble=1) ===\n";
    {
        // ダミーノード：feature=0, threshold=1e18 → 全データが左(label=1)へ
        TreeNode dummy;
        dummy.feature_id     = 0;
        dummy.threshold      = 1e18;
        dummy.left_class_id  = 1;
        dummy.right_class_id = 1;
        Evaluation(dummy, test_dataset, test_labels);
    }

    // ===== 実装(3)：深さ1の決定木 =====
    cout << "\n=== Depth-1 Decision Tree ===\n";
    TreeNode decision_tree_d1;
    TrainDecisionNode(training_dataset, training_labels, decision_tree_d1);

    cout << "Best feature: " << feature_name[decision_tree_d1.feature_id]
              << "  Threshold: " << decision_tree_d1.threshold
              << "  Left->label " << decision_tree_d1.left_class_id
              << "  Right->label " << decision_tree_d1.right_class_id << "\n";
    Evaluation(decision_tree_d1, test_dataset, test_labels);

    // ===== 実装(4)：深さ2の決定木 =====
    cout << "\n=== Depth-2 Decision Tree ===\n";
    vector<TreeNode> decision_tree_d2(3);
    TrainDecisionTree(training_dataset, training_labels, decision_tree_d2);

    // 木の構造を表示
    cout << "[Root]   feature=" << feature_name[decision_tree_d2[0].feature_id]
              << "  thr=" << decision_tree_d2[0].threshold << "\n"
              << "[Left ]  feature=" << feature_name[decision_tree_d2[1].feature_id]
              << "  thr=" << decision_tree_d2[1].threshold
              << "  L=" << decision_tree_d2[1].left_class_id
              << "  R=" << decision_tree_d2[1].right_class_id << "\n"
              << "[Right]  feature=" << feature_name[decision_tree_d2[2].feature_id]
              << "  thr=" << decision_tree_d2[2].threshold
              << "  L=" << decision_tree_d2[2].left_class_id
              << "  R=" << decision_tree_d2[2].right_class_id << "\n";
    Evaluation(decision_tree_d2, test_dataset, test_labels);

    return 0;
}