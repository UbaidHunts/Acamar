#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

const int DIM = 121; // Assuming DIM is your matrix dimension

void readMatrixFromFile(const string& file_path, float A[DIM][DIM]) {
    ifstream file(file_path);

    if (file.is_open()) {
        for (int i = 0; i < DIM; ++i) {
            string line;
            if (!getline(file, line)) {
                cerr << "Error: File contains fewer lines than expected." << endl;
                return;
            }
            istringstream iss(line);
            for (int j = 0; j < DIM; ++j) {
                if (!(iss >> A[i][j])) {
                    cerr << "Error: Invalid data format in file." << endl;
                    return;
                }
            }
        }
        file.close();
    } else {
        cerr << "Unable to open file: " << file_path << endl;
    }
}

int main() {
    string file_path = "A.txt";
    float A[DIM][DIM];
    readMatrixFromFile(file_path, A);
    cout << "Matrix read from file:" << endl;
    
    return 0;
}
