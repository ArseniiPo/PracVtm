#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using std::vector, std::cout, std::cin, std::ostream, std::istream;

namespace additional {
    constexpr long double EPS = 1e-15;
    template<typename it1, typename it2>
    long double dot_product(it1 fir, it2 sec, int len) {
        long double res = 0;
        for (int i = 0; i < len; ++i, ++fir, ++sec) {
            res += *fir * *sec;
        }
        return res;
    }
    long double
    scalar(const vector<long double> &vec1, const vector<long double> &vec2) {
        long double res = 0;
        for (size_t i = 0; i < vec1.size(); ++i) {
            res += vec1[i] * vec2[i];
        }
        return res;
    }
    // OPEARTIONS ON VECTORS
    vector<long double> operator*(long double a, const vector<long double> &vec) {
        vector<long double> res = vec;
        for (size_t i = 0; i < vec.size(); ++i) {
            res[i] *= a;
        }
        return res;
    }
    vector<long double>
    operator-(const vector<long double> &fir, const vector<long double> &sec) {   
        vector<long double> res(fir.size());
        for (size_t i = 0; i < fir.size(); ++i) {
            res[i] = fir[i] - sec[i];
        }
        return res;
    }

    void operator/=(vector<long double> &vec, long double a) {
        for (size_t i = 0; i < vec.size(); ++i) {
            vec[i] /= a;
        }
    }
    std::ostream& operator<< (std::ostream& out, const vector<long double> &A) {
        for (size_t i = 0; i < A.size(); ++i) {
            out << A[i] << " ";
        }
        return out;
    }
    // OPERATIONS ON MATRICIES
    vector<vector<long double>>
    transpose(const vector<vector<long double>> &A) {
        vector<vector<long double>> B(A[0].size(), vector<long double>(A.size(), 0));
        for (size_t i = 0; i < A.size(); ++i) {
            for (size_t j = 0; j < A[i].size(); ++j) {
                B[j][i] = A[i][j];
            }
        }
        return B;
    }
    std::ostream& operator<< (std::ostream& out, const vector<vector<long double>> &A) {
        for (size_t i = 0; i < A.size(); ++i) {
            for (size_t j = 0; j < A[i].size(); ++j) {
                out << A[i][j] << " ";
            }
            out << std::endl;
        }
        return out;
    }

    vector<vector<long double>>
    get_matrix() {
        int n, m; cin >> n >> m;
        vector<vector<long double>> mat(n, vector<long double>(m, 0)); 
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                cin >> mat[i][j];
            }
        }
        return mat;
    }
};
using namespace additional;


// -----ADDITIONAL CLASSES--------------
class Householder {
    vector<long double> v{};
    size_t n = 0;
public:
    friend class SetOfHouseholder;
    Householder(const vector<long double> &arr) : n(arr.size()) {
        v = arr;
        v /= sqrtl(scalar(v, v));
    };
    vector<long double> get_vec() { return v;  };
    vector<long double> operator* (const vector<long double> &vec) const {
        return vec - 2 * scalar(vec, v) * v; 
    }
};

class SetOfHouseholder {
    size_t n = 0;
    vector<Householder> matrix;
public:
    SetOfHouseholder() : matrix{} {}; 
    void add(Householder &hh) {
        matrix.push_back(hh);
        this->n = hh.n;
    };

    void add(const vector<long double> &vec) {
        matrix.push_back(Householder (vec));
        this->n = vec.size();
    };

    SetOfHouseholder (const SetOfHouseholder &Q) : n(Q.n), matrix(Q.matrix) {};
    vector<long double> operator*(const vector<long double> &vec) const {
        vector<long double> res = vec;
        for (size_t i = 0; i < matrix.size(); ++i) {
            res =  matrix[i] * res;
        }
        return res;
    };
    vector<long double> operator[](int index) const {
        vector<long double> tmp(n, 0);
        tmp[index] = 1; 
        return (*this) * tmp; 
    };
    friend std::ostream& operator<<(std::ostream& out, const SetOfHouseholder &Q) {
        for (size_t i = 0; i < Q.n; ++i) {
            for (size_t j = 0; j < Q.n; ++j) {
                out << Q[i][j] << " ";
            }
            out << std::endl;
        }
        return out;
    }
    void invert() { reverse(matrix.begin(), matrix.end()); };
};


class QR {
vector<vector<long double>> matrix;
public:
    QR(const vector<vector<long double>> &vec) : matrix(vec) {};
    SetOfHouseholder get_unitary() {
        SetOfHouseholder Q;
        int size = std::min(matrix[0].size() - 1, matrix.size());
        for (int i = 0; i < size - 1; ++i) {
            vector<long double> tmp(size, 0);
            copy(matrix[i].begin() + i + 1, matrix[i].end(), tmp.begin() + i);
            Q.add(tmp);
        }
        return Q;
    };
    vector<vector<long double>>
    get_traingular() {
        vector<vector<long double>> result(matrix.size(), vector<long double>(matrix[0].size() - 1, 0));
        size_t size = std::min(matrix[0].size() - 1, matrix.size());
        for (size_t i = 0; i < matrix.size(); ++i) {
            size_t len = std::min(i + 1, size);
            copy(matrix[i].begin(), matrix[i].begin() + len, result[i].begin());
        }
        return transpose(result);
    }
};

//----------ALGORITHM--------------
void
multiply(vector<long double> &vec, const vector<long double> &refl, int i) {
    long double product = dot_product(vec.begin() + i, refl.begin(), vec.size() - i);
    for (size_t j = 0; j < refl.size(); ++j) {
        vec[i + j] -= 2 * product * refl[j];
    }
} 

vector<long double>
step_of_qr(vector<long double> &vec, int i) {
    vector<long double> refl(vec.begin() + i, vec.end());
    long double coef = sqrtl(scalar(refl, refl));
    if (refl[0] < 0) {
        coef *= -1;
    }
    refl[0] += coef;
    refl /= sqrtl(scalar(refl, refl));
    vec[i] = -coef;
    for_each(vec.begin() + i + 1, vec.end(), [](long double &a){ a = 0;});
    return refl;
}


QR
qr_solve(const vector<vector<long double>> &mat) {
    vector<vector<long double>> A = transpose(mat);
    for (size_t i = 0; i < std::min(A[i].size(), A.size()) - 1; ++i) {
        vector<long double> refl = step_of_qr(A[i], i);
        cout << refl << std::endl;
        A[i].push_back(0);
        copy(refl.begin(), refl.end(), A[i].begin() + i + 1);
        for (size_t j = i + 1; j < A.size(); ++j) {
            multiply(A[j], refl, i);
        }
    }
    return QR(A);
}

int main(void)
{
    vector<vector<long double>> mat = get_matrix();
    QR res = qr_solve(mat);
    SetOfHouseholder Q = res.get_unitary();
    vector<vector<long double>> R = res.get_traingular();
    cout << "Q\n" << Q << "\nR\n" << R << std::endl;
    /* //----------CHECK----------------
    vector<vector<long double>> A = transpose(R);
    Q.invert();
    for (int i = 0; i < A.size(); ++i) {
        A[i] = Q * A[i];
        for_each(A[i].begin(), A[i].end(), [](long double &a){if (fabs(a) < EPS) a = 0;});
    }
    A = transpose(A);
    long double result = 0;
    for (int i = 0; i < A.size(); ++i) {
        result += scalar(A[i] - mat[i], A[i] - mat[i]);
    }
    cout << "FROBENIOUS: " << sqrtl(result) << std::endl; */
    return 0; 
}
