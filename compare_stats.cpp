#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

vector<vector<int>> readMatrix(const string& filename) {
    vector<vector<int>> matrix;
    ifstream file(filename);
    string line;

    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }
    matrix.push_back(vector<int>(40, 0));
    while (getline(file, line)) {
        stringstream ss(line);
        vector<int> row;
        int num;
        row.push_back(0);
        while (ss >> num) {
            row.push_back(num);
        }
        matrix.push_back(row);
    }

    file.close();
    return matrix;
}

void compare_to_best(int min, int max, const string& name1, const string& name2, const vector<vector<int>>& matrix1, const vector<vector<int>>& matrix2, ostream& out=cout){
    int num0=0;
    int num1=0;
    int num2=0;
    int num3=0;
    for(int i=min-1;i<max;i++){
        for(int j=min-1;j<=i;j++){
            if(matrix1[i][j]==matrix2[i][j]) num0++;
            else if(matrix1[i][j]+1==matrix2[i][j]) num1++;
            else if(matrix1[i][j]+2==matrix2[i][j]) num2++;
            else if(matrix1[i][j]+3==matrix2[i][j]) num3++;
        }
    }
    cout << name1 << ":\n";
    cout << "Off by zero: " << num0 << "times" <<endl;
    cout << "Off by one: " << num1 << "times" <<endl;
    cout << "Off by two: " << num2 << "times" <<endl;
    cout << "Off by three: " << num3 << "times" <<endl;
}

void print_special_values(int min, int max, const string& name1, const string& name2, const vector<vector<int>>& matrix1, const vector<vector<int>>& matrix2, ostream& out=cout){
    cout << name1 << ":\n";
    cout << "Z(9,30): " << matrix1[29][8] << endl;
    cout << "Z(14,29): " << matrix1[28][13] << endl;
    cout << "Z(20,20): " << matrix1[19][19] << endl;
    cout << "Z(20,37): " << matrix1[36][19] << endl;
    cout << "Z(30,30): " << matrix1[29][29] << endl;
    cout << "Z(40,40): " << matrix1[39][39] << endl;

    cout << endl << name2 << ":\n";
    cout << "Z(9,30): " << matrix2[29][8] << endl;
    cout << "Z(14,29): " << matrix2[28][13] << endl;
    cout << "Z(20,20): " << matrix2[19][19] << endl;
    cout << "Z(20,37): " << matrix2[36][19] << endl;
    cout << "Z(30,30): " << matrix2[29][29] << endl;
    cout << "Z(40,40): " << matrix2[39][39] << endl;
}

void compare_edgenum(int min, int max, const string& name1, const string& name2, const vector<vector<int>>& matrix1, const vector<vector<int>>& matrix2, ostream& out=cout){
    int sum1=0;
    int sum2=0;
    for(int i=min-1;i<max;i++){
        for(int j=min-1;j<max;j++){
            sum1+=matrix1[i][j];
            sum2+=matrix2[i][j];
        }
    }
    cout << name1 << ": "<<sum1 <<endl;
    cout << name2 << ": "<<sum2 <<endl;
    cout << "Result1" << ": " <<(double)sum1/sum2<<endl;
    cout << "Result2" << ": " <<(double)sum2/sum1<<endl;
}

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

int main(int argc, char* argv[]) {

    if (argc !=3 && (argc != 5 || !is_number(argv[3]) || !is_number(argv[4]))) {
        cerr << "Usage: " << argv[0] << " <input1.txt> <input2.txt> [<min size> <max size>]" << endl;
        return 1;
    }
    cout << argc <<endl;
    int min = 2;
    int max = 40;
    if(argc==5){
        min = stoi( argv[3] );
        max = stoi( argv[4] );
        if (min>max || min<2 || max>50) {
            cerr << "third parameter can't be greater than fourth parameter and both parameters must be in the range of: 2-50" << endl;
            return 2;
        }
    }
    string ending = "_results.txt";
    string file1 = argv[1]+ending; // First input file
    string file2 = argv[2]+ending; // Second input file

    auto matrix1 = readMatrix(file1);
    auto matrix2 = readMatrix(file2);
    if(matrix1.size()<max) cerr << "the size of" << argv[1] << "is smaller than max: " << max << endl;
    if(matrix2.size()<max) cerr << "the size of" << argv[2] << "is smaller than max: " << max << endl;

    compare_edgenum(min, max, argv[1], argv[2], matrix1, matrix2);
    cout << endl << "Values for table: " <<endl;
    print_special_values(min, max, argv[1], argv[2], matrix1, matrix2);
    cout << endl << "Some stats:" << endl;
    compare_to_best(min, max, argv[1], argv[2], matrix1, matrix2);

    return 0;
}