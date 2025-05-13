#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

long long sumMatrix(const string& filename) {
    vector<vector<int>> matrix;
    ifstream file(filename);
    string line;

    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }

    long long sum=0;
    while (getline(file, line)) {
        stringstream ss(line);
        int num;
        while (ss >> num) {
            sum +=num;
        }
    }

    file.close();
    return sum;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input1.txt> <input2.txt>" << endl;
        return 1;
    }

    string ending = "_results.txt";
    string file1 = argv[1]+ending; // First input file
    string file2 = argv[2]+ending; // Second input file

    long long sum1 = sumMatrix(file1);
    long long sum2 = sumMatrix(file2);

    cout << argv[1] << ": "<<sum1 <<endl;
    cout << argv[2] << ": "<<sum2 <<endl;
    cout << "Result1" << ": " <<(double)sum1/sum2<<endl;
    cout << "Result2" << ": " <<(double)sum2/sum1<<endl;

    return 0;
}