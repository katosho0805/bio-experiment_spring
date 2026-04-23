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

#define NUM_FEATURES 53
#define NUM_SEQS   10000


//決定木モデルの構造体
struct TreeNode{
    int feature_id;
    double threshold;//閾値
    int left_class_id;
    int right_class_id;
};  

//ファイルopen関数のために構造体を作る
void LoadSolubilityFile(const string& filename,
    vector<string>&feature_name,
    vector<vector<double>>&dataset,
    vector<int>&labels){
        ifstream ifs(filename);
        if (!ifs) {
        cerr << "Error: Cannot open file: " << filename << endl;
        return;
    }

    string line;

    getline(ifs, line);          // 1行読み込む（ヘッダー行）
    istringstream ss(line);      // 空白などを飛ばしてくれる
    string token;
    ss >> token;                 // "id" を読み捨てる
     for (int i = 0; i < NUM_FEATURES; i++) {
       ss >> token;
       feature_name[i] = token;
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
 

//データセット分割 
void DivideDataset(const vector<vector<double>>& dataset,
    const vector<int>&                 labels,
    vector<vector<double>>&       training_dataset,
    vector<int>&                       training_labels,
    vector<vector<double>>&       test_dataset,
    vector<int>&                       test_labels,
    double test_ratio){
      vector <int> index;
      for(int i=0;i<(int)dataset.size();i++){
        index.push_back(i);
      }
    // 乱数生成
    mt19937 gen(42);

    // シャッフル実行
    shuffle(index.begin(), index.end(), gen);
    for(int i=0;i<(int)index.size();i++){
        if(i<test_ratio*dataset.size()){
            test_dataset.push_back(dataset[index[i]]);
            test_labels.push_back(labels[index[i]]);//テストデータセットに追加
        }else{
            training_dataset.push_back(dataset[index[i]]);
            training_labels.push_back(labels[index[i]]);//トレーニングデータセットに追加
        }
      } 
    }

// ジニ不純度の計算
double GiniImpurity(const vector<int>& labs) {
    int n = labs.size();
    int l = 0;
    for(int i = 0; i < (int)labs.size(); i++){
        if(labs[i] == 1) l++;
    }
    double p1 = (double)l / n;
    double p0 = 1.0 - p1;
    return 1.0 - p1*p1 - p0*p0;
}

// 重み付きジニ不純度の計算
double WeightedGini(const vector<int>& left_labs, const vector<int>& right_labs){
    int left_size  = left_labs.size();
    int right_size = right_labs.size();
    if(left_size == 0 || right_size == 0) return 1.0;

    double total_size = left_size + right_size;
    return (left_size  / total_size) * GiniImpurity(left_labs)
         + (right_size / total_size) * GiniImpurity(right_labs);
}

// 多数決ラベルの取得
int MajorityLabel(const vector<int>& labs){
    int pos = 0;
    for(int l : labs) pos += l;
    return (pos * 2 >= (int)labs.size()) ? 1 : 0;
}

// 1ノードの学習
void TrainDecisionNode(const vector<vector<double>>& data, const vector<int>& labs, TreeNode& node){
    int n = data.size();
    double best_gini = 1.0;
    int    best_feature_id = 0;
    double best_threshold  = 0.0;
    int    best_left  = 0;
    int    best_right = 1;

    for(int feature_id = 0; feature_id < NUM_FEATURES; feature_id++){
        // 特徴量の値を取り出してソート（パーセンタイル閾値用）
        vector<double> feature_values;
        for(int i = 0; i < n; i++){
            feature_values.push_back(data[i][feature_id]);
        }
        sort(feature_values.begin(), feature_values.end());

        // 1〜99パーセンタイルを閾値候補として評価
        for(int p = 1; p <= 99; p++){
            int pos_idx = (int)((double)p / 100.0 * (n - 1));
            double threshold = feature_values[pos_idx];

            vector<int> left_labs, right_labs;
            for(int j = 0; j < n; j++){
                if(data[j][feature_id] <= threshold){
                    left_labs.push_back(labs[j]);
                }else{
                    right_labs.push_back(labs[j]);
                }
            }

            double gini = WeightedGini(left_labs, right_labs);
            if(gini < best_gini){
                best_gini       = gini;
                best_feature_id = feature_id;
                best_threshold  = threshold;
                best_left       = MajorityLabel(left_labs);
                best_right      = MajorityLabel(right_labs);
            }
        }
    }

    node.feature_id     = best_feature_id;
    node.threshold      = best_threshold;
    node.left_class_id  = best_left;
    node.right_class_id = best_right;
}

// 深さ2の決定木の学習
void TrainDecisionTree(const vector<vector<double>>& data, const vector<int>& labs, vector<TreeNode>& dtree){
    // 親ノード(dtree[0])を学習
    TrainDecisionNode(data, labs, dtree[0]);

    // 親ノードでデータを左右に分割
    vector<vector<double>> left_data, right_data;
    vector<int>            left_labs, right_labs;
    for(int i = 0; i < (int)data.size(); i++){
        if(data[i][dtree[0].feature_id] <= dtree[0].threshold){
            left_data.push_back(data[i]);
            left_labs.push_back(labs[i]);
        }else{
            right_data.push_back(data[i]);
            right_labs.push_back(labs[i]);
        }
    }

    // 左子ノード(dtree[1])・右子ノード(dtree[2])を学習
    TrainDecisionNode(left_data,  left_labs,  dtree[1]);
    TrainDecisionNode(right_data, right_labs, dtree[2]);
}

// 性能指標を表示するヘルパー
void PrintMetrics(int tp, int fp, int fn, int tn){
    double accuracy  = (double)(tp + tn) / (tp + fp + fn + tn);
    double precision = (tp + fp > 0) ? (double)tp / (tp + fp) : 0.0;
    double recall    = (tp + fn > 0) ? (double)tp / (tp + fn) : 0.0;
    double f1_score  = (precision + recall > 0)
                       ? 2.0 * precision * recall / (precision + recall) : 0.0;

    cout << fixed << setprecision(4);
    cout << "Accuracy:  " << accuracy  << "\n"
         << "Precision: " << precision << "\n"
         << "Recall:    " << recall    << "\n"
         << "F1-score:  " << f1_score  << "\n"
         << "Confusion Matrix\n"
         << "TP: " << tp << "  FP: " << fp << "\n"
         << "FN: " << fn << "  TN: " << tn << "\n";
}

// 深さ1の評価
void Evaluation(const TreeNode& node, const vector<vector<double>>& test_data, const vector<int>& test_labels){
    int tp = 0, fp = 0, fn = 0, tn = 0;
    for(int i = 0; i < (int)test_data.size(); i++){
        int pred   = (test_data[i][node.feature_id] <= node.threshold)
                     ? node.left_class_id : node.right_class_id;
        int actual = test_labels[i];
        if      (pred == 1 && actual == 1) tp++;
        else if (pred == 1 && actual == 0) fp++;
        else if (pred == 0 && actual == 1) fn++;
        else                               tn++;
    }
    PrintMetrics(tp, fp, fn, tn);
}

// 深さ2の評価
void Evaluation(const vector<TreeNode>& dtree, const vector<vector<double>>& test_data, const vector<int>& test_labels){
    int tp = 0, fp = 0, fn = 0, tn = 0;
    for(int i = 0; i < (int)test_data.size(); i++){
        int pred;
        if(test_data[i][dtree[0].feature_id] <= dtree[0].threshold){
            const TreeNode& child = dtree[1];
            pred = (test_data[i][child.feature_id] <= child.threshold)
                   ? child.left_class_id : child.right_class_id;
        }else{
            const TreeNode& child = dtree[2];
            pred = (test_data[i][child.feature_id] <= child.threshold)
                   ? child.left_class_id : child.right_class_id;
        }
        int actual = test_labels[i];
        if      (pred == 1 && actual == 1) tp++;
        else if (pred == 1 && actual == 0) fp++;
        else if (pred == 0 && actual == 1) fn++;
        else                               tn++;
    }
    PrintMetrics(tp, fp, fn, tn);
}


int main(void){
    vector<string>feature_name(NUM_FEATURES,"");
    vector<vector<double>>dataset(NUM_SEQS,vector<double>(NUM_FEATURES,0.0));
    vector<int>labels(NUM_SEQS);

    LoadSolubilityFile("protein_solubility_dataset.txt",feature_name,dataset,labels);

    vector<vector<double>>training_dataset;
    vector<int>training_labels;
    vector<vector<double>>test_dataset;
    vector<int>test_labels;
    double test_ratio=0.2;

    DivideDataset(dataset,labels,training_dataset,training_labels,test_dataset,test_labels,test_ratio);

    cout << "Training dataset size: " << training_dataset.size() << "\n";
    cout << "Test dataset size: "     << test_dataset.size()     << "\n\n";

    // 深さ1の決定木
    cout << "=== Depth-1 Decision Tree ===\n";
    TreeNode decision_tree_d1;
    TrainDecisionNode(training_dataset, training_labels, decision_tree_d1);
    cout << "Feature: " << feature_name[decision_tree_d1.feature_id]
         << "  Threshold: " << decision_tree_d1.threshold
         << "  Left=" << decision_tree_d1.left_class_id
         << "  Right=" << decision_tree_d1.right_class_id << "\n";
    Evaluation(decision_tree_d1, test_dataset, test_labels);

    // 深さ2の決定木
    cout << "\n=== Depth-2 Decision Tree ===\n";
    vector<TreeNode> decision_tree_d2(3);
    TrainDecisionTree(training_dataset, training_labels, decision_tree_d2);
    cout << "[Root]  feature=" << feature_name[decision_tree_d2[0].feature_id]
         << "  thr=" << decision_tree_d2[0].threshold << "\n"
         << "[Left]  feature=" << feature_name[decision_tree_d2[1].feature_id]
         << "  thr=" << decision_tree_d2[1].threshold
         << "  L=" << decision_tree_d2[1].left_class_id
         << "  R=" << decision_tree_d2[1].right_class_id << "\n"
         << "[Right] feature=" << feature_name[decision_tree_d2[2].feature_id]
         << "  thr=" << decision_tree_d2[2].threshold
         << "  L=" << decision_tree_d2[2].left_class_id
         << "  R=" << decision_tree_d2[2].right_class_id << "\n";
    Evaluation(decision_tree_d2, test_dataset, test_labels);

    return 0;
}