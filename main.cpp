#include <bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define ld long double
#define print_lst_1x(lst) for (auto el: lst){cout << el << ' ';} cout << endl;
#define print_lst_2x(lst) for (auto el: lst){for (auto zn: el) {cout << zn << ' ';} cout << '\n';} cout << endl;
using namespace std;



struct TA {
    ld n;

    TA(){};

    TA(ld k){
        n = k;
    }

    // переопределили + как минимум
    TA operator+(const TA a) {
        TA result(*this);
        return a.n <= result.n ? a : result;
    }

    // переопределили * как сложение
    TA operator*(const TA a) {
        TA result(*this);
        return TA(result.n + a.n);
    }
};

std::ostream& operator << (std::ostream& o, const TA& a){
    o << a.n;
    return o;
}



class Matrix;

class SparseMatrix{
public:
    int n;
    vector <TA> values;
    vector <int> columns;
    vector <int> row_begin;

    SparseMatrix();
    SparseMatrix(vector <TA> &val, vector <int> &col, vector <int> &rw_bg);
    SparseMatrix(Matrix &M);
    SparseMatrix square();
};


class Matrix {
public:
    int n;
    vector <vector <TA>> matrix;


    Matrix(){};

    Matrix(vector <vector <TA>> &mx){
        matrix = mx;
        n = matrix.size();
    }

    /// Меняем разреженную структуру матрицы на дефолтную структуру
    Matrix(SparseMatrix &SM){
        n = SM.n;
        matrix.resize(n, vector <TA> (n, INFINITY));
        for (int i = 0; i < n; i++){
            for (int j = SM.row_begin[i]; j < SM.row_begin[i + 1]; j++){
                matrix[i][SM.columns[j]] = SM.values[j];
            }
        }
    }


    void resize(int k){
        n = k;
        matrix = vector <vector <TA>> (n, vector<TA>(n, INFINITY));
    }

    /// переопределили умножение матриц
    Matrix operator*(const Matrix& B){
        int N = this->n;
        Matrix ANS;
        ANS.resize(n);

        for (int i = 0; i < N; ++i){
            const TA* aRow = this->matrix[i].data();
            TA* cRow = ANS.matrix[i].data();

            for (int k = 0; k < N; ++k){
                const TA* bRow = B.matrix[k].data();
                TA aVal = aRow[k];

                for (int j = 0; j < N; ++j){
                    cRow[j] = cRow[j] + (aVal * bRow[j]);
                }
            }
        }
        return ANS;
    }
};

void print(Matrix M){
    for (int i = 0; i < M.n; i++){
        for (int j = 0; j < M.n; j++){
            cout << M.matrix[i][j] << '\t';
        }
        cout << '\n';
    }
    cout << '\n';
}


SparseMatrix::SparseMatrix(){}

SparseMatrix::SparseMatrix(vector <TA> &val, vector <int> &col, vector <int> &rw_bg){
    n = rw_bg.size() - 1;
    values = val;
    columns = col;
    row_begin = rw_bg;
}

/// Меняем дефолтную структуру матрицы на разреженную структуру
SparseMatrix::SparseMatrix(Matrix &M){
    n = M.n;

    row_begin.resize(n + 1);
    int value_size = 0;

    for (int i = 0; i < n; i++){
        bool flag = true;
        for (int j = 0; j < n; j++){
            if (M.matrix[i][j].n != INFINITY){
                if (flag){
                    row_begin[i] = value_size;
                    flag = false;
                }
                values.push_back(M.matrix[i][j]);
                value_size++;
                columns.push_back(j);
            }
        }
    }
    row_begin[n] = value_size;
}

/// Возведение в квадрат
SparseMatrix SparseMatrix::square(){
    vector <TA> new_values;
    vector <int> new_columns;
    vector <int> new_row_begin;

    int new_values_size = 0;
    for (int i = 0; i < n; ++i) {
        new_row_begin.push_back(new_values_size);

        for (int j = 0; j < n; ++j){
            TA sum = INFINITY;

            for (int p = row_begin[i]; p < row_begin[i + 1]; ++p){
                TA val = values[p];

                int l = row_begin[columns[p]], r = row_begin[columns[p] + 1];
                int m;
                while (r - l > 1) {
                    m = (r + l) / 2;
                    if (columns[m] > j) {
                        r = m;
                    } else {
                        l = m;
                    }
                }

                if (columns[l] == j) {
                    sum = sum + (val * values[l]);
                }
            }

            if (sum.n != INFINITY) {
                new_values.push_back(sum);
                new_values_size++;
                new_columns.push_back(j);
            }
        }
    }
    new_row_begin.push_back(new_values_size);

    return SparseMatrix(new_values, new_columns, new_row_begin);
}

void print(SparseMatrix &SM){
    for (TA el: SM.values){
        cout << el << ' ';
    }
    cout << '\n';
    for (int el: SM.columns){
        cout << el << ' ';
    }
    cout << '\n';
    for (int el: SM.row_begin){
        cout << el << ' ';
    }
    cout << '\n';
}



class general_matrix_structure{
public:
    Matrix M;
    SparseMatrix SM;
    int dimention;
    int non_zero_elements_count;

    general_matrix_structure(){};
};



class Graph{
public:
    // количество вершин
    int ver_count;
    // матрица смежности
    general_matrix_structure adj_matrix;
    // матрица минимальных путей
    general_matrix_structure min_path_matrix;


    Graph(const string& file_name){
        ifstream fin(file_name);

        if (fin.is_open()){
            fin >> ver_count;
            adj_matrix.dimention = ver_count;
            min_path_matrix.dimention = ver_count;

            adj_matrix.M.resize(ver_count);
            min_path_matrix.M.resize(ver_count);

            for (int i = 0; i < ver_count; i++){
                adj_matrix.M.matrix[i][i] = 0;
                min_path_matrix.M.matrix[i][i] = 0;
            }

            int start, finish;
            ld value;
            int edge_count = ver_count;

            while(fin >> start >> finish >> value){
                adj_matrix.M.matrix[start][finish] = value;
                min_path_matrix.M.matrix[start][finish] = value;
                edge_count++;
                if (start == finish && value == 0){
                    edge_count--;
                }
            }
            adj_matrix.non_zero_elements_count = edge_count;
            min_path_matrix.non_zero_elements_count = edge_count;
        }
        fin.close();
//        get_min_path_matrix();
    }


    /// Функция для определения перехода с одной структуры на другую
    int transition_function(int x){
        return min(x * (x - 1), int(1.15105900938376493059 * pow(x, 1.94433481614137959603)));
    }

    /// Создание матрицы минимальных путей
    void get_min_path_matrix(){
        int limit_value = transition_function(ver_count);

        if (min_path_matrix.non_zero_elements_count < limit_value){
            min_path_matrix.SM = min_path_matrix.M;
        }

        int degree = 1;
        while (min_path_matrix.SM.values.size() < limit_value){
            min_path_matrix.SM = min_path_matrix.SM.square();
            degree *= 2;
            cout << degree << endl;
            if (degree >= ver_count){
                min_path_matrix.M = min_path_matrix.SM;
                return;
            }
        }
        min_path_matrix.M = min_path_matrix.SM;
        while (degree < ver_count){
            min_path_matrix.M = min_path_matrix.M * min_path_matrix.M;
            degree *= 2;
            cout << degree << endl;
        }
        // не создаю SM, хз мб и не пригодится...
    }

    /// Получение пути
    vector <int> get_path(int s, int t){
        if (min_path_matrix.M.matrix[s][t].n == INFINITY){
            cout << "NO WAY!\n";
            return {};
        }
        vector <int> path = {s};
        int sn = s;
        while(sn != t){
            for (int i = 0; i < ver_count; ++i){
                if (i != sn && adj_matrix.M.matrix[sn][i].n != INFINITY &&
                (adj_matrix.M.matrix[sn][i].n + min_path_matrix.M.matrix[i][t].n) == min_path_matrix.M.matrix[sn][t].n){
                    path.push_back(i);
                    sn = i;
                    break;
                }
            }
        }
        return path;
    }
};

void print(Graph G){
    for (int i = 0; i < G.ver_count; i++){
        for (int j = 0; j < G.ver_count; j++){
            cout << G.adj_matrix.M.matrix[i][j] << '\t';
        }
        cout << '\n';
    }
    cout << '\n';
}



string generate_random_graph(int vertices, int edges){
    // Генератор случайных чисел
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> val_dist(1000, 9999); // Веса рёбер
    uniform_int_distribution<> idx_dist(0, vertices - 1); // Индексы вершин

    // Матрица смежности
    vector<vector<int>> adjacency_matrix(vertices, vector<int>(vertices, 0));

    // Генерация рёбер
    int generated_edges = 0;
    while (generated_edges < edges) {
        int i = idx_dist(gen);
        int j = idx_dist(gen);

        // Проверяем, чтобы не было петли (i != j) и чтобы ребро не было уже создано
        if (i != j && adjacency_matrix[i][j] == 0) {
            adjacency_matrix[i][j] = val_dist(gen); // Присваиваем случайный вес
            generated_edges++;
        }
    }

    string filename = "graph_test.txt";
    ofstream outfile(filename);
    if (outfile.is_open()) {
        outfile << vertices << endl;
        for (int i = 0; i < vertices; ++i) {
            for (int j = 0; j < vertices; ++j) {
                if (adjacency_matrix[i][j] != 0) {
                    outfile << i << " " << j << " " << adjacency_matrix[i][j] << "\n";
                }
            }
        }
        outfile.close();
    }

    return filename;
}



/// Проверка общей работоспособности
void experiment_1(){
    /// эксперимент времени
    auto start1 = std::chrono::high_resolution_clock::now();

    Graph G ("C:\\Users\\Mukharyam\\CLionProjects\\Huawei_practice\\graph_2000_5105.txt");

    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = end1 - start1;
    cout << duration1.count() << endl;



    for (int el: G.get_path(0, 1)){
        cout << el << ' ';
    }
    cout << '\n';


    for (int el: G.get_path(2, 1)){
        cout << el << ' ';
    }
    cout << '\n';
}


/// Поиск момента перехода
void experiment_2(){
    for (int i = 100; i < 300; i+=50){
        int m, l = int(1.15 * pow(i, 1.9)), r = i * (i - 1) + 1;
        ld dur1, dur2;
        while(r - l > 1){
            m = (r + l) / 2;
            Graph g(generate_random_graph(i, m));

            /// эксперимент времени
            auto start1 = std::chrono::high_resolution_clock::now();

            g.adj_matrix.M * g.adj_matrix.M;

            auto end1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration1 = end1 - start1;
            dur1 = duration1.count();


            /// эксперимент времени
            auto start2 = std::chrono::high_resolution_clock::now();

            g.adj_matrix.SM.square();

            auto end2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration2 = end2 - start2;
            dur2 = duration2.count();

            cout << "!!!!\n";

            if (dur1 <= dur2){
                r = m;
            } else {
                l = m;
            }
        }
        cout << i << ' ' << m << endl;
        cout << dur1 << ' ' << dur2 << endl << endl;
    }
}



void experiment_3(){
    Graph g ("C:\\Users\\Mukharyam\\CLionProjects\\Huawei_practice\\graph_3.txt");
//    print(g);

    ld dur;

    cout << "Go!\n";
    /// эксперимент времени
    auto start1 = std::chrono::high_resolution_clock::now();

    g.get_min_path_matrix();

    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = end1 - start1;
    dur = duration1.count();

    cout << dur << endl;
}





int solve() {
//    Graph g ("C:\\Users\\Mukharyam\\CLionProjects\\Huawei_practice\\graph_6.txt");
    experiment_3();

}



int main(){
    solve();
//    generate_random_graph(100, 9000);
}

