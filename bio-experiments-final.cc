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
 void DivideDataset(    const vector<vector<double>>& dataset,
    const vector<int>&                 labels,
    vector<vector<double>>&       training_dataset,
    vector<int>&                       training_labels,
    vector<vector<double>>&       test_dataset,
    vector<int>&                       test_labels,
    double test_ratio){
      vector <int> index;
      for(int i=0;i<dataset.size();i++){
        index.push_back(i);
      }
    // 乱数生成
    mt19937 gen(42);

    // シャッフル実行
    shuffle(index.begin(), index.end(), gen);
    for(int i=0;i<index.size();i++){
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
       int pos = 0;
       for (int l : labs) pos += l;
       double p = (double)pos / labs.size();
       return 2.0 * p * (1.0 - p);
}

    //重み付きジニ不純度の計算
    double weight_gini(const vector<int>& left_labs,const vector<int>& right_labs){
        int left_size=left_labs.size();
        int right_size=right_labs.size();
        if(left_size==0 || right_size==0){
             return 1.0;
        }
         
        double left_gini=1.0;
        double right_gini=1.0;
        double gini(const vector<int>& labs) {
          int n = labs.size();
          int l=0;
          for(int i=0;i<labs.size();i++){
            if(labs[i]==1){
                l++;
            }
          }
          double p1 = (double)l / n;
          double p0 = 1.0 - p1;
          return 1.0 - p1*p1 - p0*p0;
      }
        left_gini=gini(left_labs);
        right_gini=gini(right_labs);
        double total_size=left_size+right_size;
        return (left_size/total_size)*left_gini+(right_size/total_size)*right_gini;
    }


    //TrainDecisionNode関数の実装
    void TrainDecisionNode(const vector<vector<double>>& data,const vector<int>& labs,TreeNode& node){
       int n=data.size();
         double best_gini=1.0;
    int best_feature_id=0;
    double best_threshold=0.0;
    for(int feature_id=0;feature_id<NUM_FEATURES;feature_id++){
        vector<double> feature_values;
        for(int i=0;i<n;i++){
            feature_values.push_back(data[i][feature_id]);
        }
        // 特徴量の値を昇順にソート
        sort(feature_values.begin(), feature_values.end());
        // 閾値を決定
        for(int i=0;i<feature_values.size()-1;i++){
            double threshold=(feature_values[i]+feature_values[i+1])/2.0;
            // 閾値でデータを分割
            vector<int> left_labs, right_labs;
            for(int j=0;j<n;j++){
                if(data[j][feature_id]<threshold){
                    left_labs.push_back(labs[j]);
                }else{
                    right_labs.push_back(labs[j]);
                }
            }
            // 重み付きジニ不純度を計算
            double gini=weight_gini(left_labs,right_labs);
            // 最も低いジニ不純度を持つ分割を選択
            if(gini<best_gini){
                best_gini=gini;
                best_feature_id=feature_id;
                best_threshold=threshold;
            }
        }
    }
    // 最適な分割情報をノードに保存
    node.feature_id=best_feature_id;
    node.threshold=best_threshold;
    }
//深さ1の決定木の評価 
void Evaluation(const TreeNode& node,const vector<vector<double>>& test_data,const vector<int>& test_labels){
    int tp=0; int fp=0; int fn=0; int tn=0;
    double accuracy=0.0; 
    double precision=0.0;
    double recall=0.0;
    double f1_score=0.0;
    for(int i=0;i<left_labs.size();i++){
        


        

    }
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

    cout << "Training dataset size: " << training_dataset.size() << endl;
    cout << "Test dataset size: " << test_dataset.size() << endl;

    return 0;

}