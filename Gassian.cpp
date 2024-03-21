#include <vector>
#include <iostream>
using namespace std;

class Fract{
public:
    static int gcd(int a, int b){
        a = abs(a);
        b = abs(b);
        return b?gcd(b, a % b) : a;
    }
    static int lcm(int a, int b){
        return a * (b / gcd(a, b)); 
    }
    static const Fract Zero;

    int first = 0;
    int second = 0;

    bool IsZero() const {
        return Fract2Float() == 0;
    }
    float Fract2Float() const {
        if (second != 0)
            return (float) first / (float) second;
        else
            return 0;
    }
    const Fract GetReciprocal() const {
        return Fract(second, first);
    }
    void ReduceFraction(){
        int g = gcd(first, second); 
        if (g == 0)
            return;
        first /= g;
        second /= g;
    }

    Fract(){};
    Fract(float num){
        second = 1;
        while(num != (int) num){
            num *= 10;
            second *= 10;
        }

        first = num;

        this->ReduceFraction();
    }
    Fract(int first, int second){
        this->first = first;
        this->second = second;
        this->ReduceFraction();
    }

    Fract operator+(const Fract& other){
        int _lcm = lcm(second, other.second);

        int time1 = _lcm / second;
        int time2 = _lcm / other.second;

        int new_first = first * time1 + other.first * time2;

        Fract f = Fract(new_first, _lcm);
        f.ReduceFraction();
        return f;
    }

    Fract operator-(const Fract& other){
        Fract t = other;
        t.first *= -1;
        return *this + t;
    }
    Fract operator*(const Fract& other){
        int new_first = first * other.first; 
        int new_second = second * other.second;
        int g = gcd(new_first, new_second);

        if (g == 0)
            return Zero;
        
        new_first /= g;
        new_second /= g;

        return Fract(new_first, new_second);
    }
    Fract operator / (const Fract& other){
        if (other.IsZero()){
            cerr << "Divide by zero" << endl;
            return Zero;
        }
        return *this * other.GetReciprocal();
    }
};
const Fract Fract::Zero = Fract(0, 0);



// Swap the "first" and "second" line
template<typename T>
void SwapLine(vector<vector<T>>& matrix,  int first, int second){
    int n = matrix.size();
    if (first >= n || second >= n)
        return;
    
    swap(matrix[first], matrix[second]);
}


vector<vector<float>> GassianElimination(vector<vector<float>> matrix){
    if (matrix.size() <= 0 || matrix[0].size() <= 0)
        return vector<vector<float>>();

    int n = matrix.size(), m = matrix[0].size();

    vector<vector<Fract>> fract(n, vector<Fract>(m));

    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            float num = matrix[i][j];
            fract[i][j] = Fract(num);
        }
    }

    // In fact, true rank may be smaller
    int rank = min(n, m);

    for (int i = 0; i < rank; i++){
        // Find first not zero on pivot
        int j = i; 
        while(j < rank && fract[j][i].IsZero()) 
            j++;
        // Can not find no zero
        if (j >= rank)
            break;   
        SwapLine(fract, i, j);        

        // normalized
        Fract reci = fract[i][i].GetReciprocal();
        for (int k = i; k < m; k++)
            fract[i][k] = fract[i][k] * reci;
        
        //elimination
        for (int k = i + 1; k < n; k++){
            if (!fract[k][i].IsZero()){
                Fract time = fract[i][i] / fract[k][i]; 

                for (int z = i; z < m; z++){
                    fract[k][z] = fract[k][z] * time - fract[i][z];
                }
            }
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            matrix[i][j] = fract[i][j].Fract2Float();
    
    return matrix;
} 


#include <fstream> // Include the necessary header file
#include <sstream>

vector<vector<float>> ReadMatrix(const char* filename = "matrix.txt"){
    ifstream file(filename);
    vector<vector<float>> matrix;
    string line;

    while(getline(file, line)){
        vector<float> row;
        stringstream ss(line);
        float value;

        while(ss >> value)
            row.push_back(value);
        matrix.push_back(row); 
    }

    return matrix;
}

void WriteMatrix(const vector<vector<float>>& matrix, const char* filename = "matrix_gassian.txt"){
    ofstream file(filename);
    for (auto row : matrix){
        for (auto value : row)
            file << value << " ";
        file << endl;
    }
}

int main(){

    WriteMatrix(GassianElimination(ReadMatrix("matrix")));
    return 0;
}

