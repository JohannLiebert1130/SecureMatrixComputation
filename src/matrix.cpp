#include "matrix.h"
#include <cstring>
#include "FHEContext.h"
#include "EvalMap.h"
#include "powerful.h"
#include "binio.h"
#include <fstream>
#include <algorithm>


long FindM1(long k, long L, long c, long p, long d, long s, long chosen_m, int dim, bool verbose=false)
{
    // get a lower-bound on the parameter N=phi(m):
    // 1. Each level in the modulus chain corresponds to pSize=p2Size/2
    //    bits (where we have one prime of this size, and all the others are of
    //    size p2Size).
    //    When using DoubleCRT, we need 2m to divide q-1 for every prime q.
    // 2. With nBits of ctxt primes, 
    //    the largest modulus for "fresh ciphertexts" has size
    //          Q0 ~ 2^{nBits}
    // 3. We break each ciphertext into upto c digits, do each digit is as large
    //    as    D=2^{nBits/c}
    // 4. The added noise variance term from the key-switching operation is
    //    c*N*sigma^2*D^2, and this must be mod-switched down to w*N (so it is
    //    on par with the added noise from modulus-switching). Hence the ratio
    //    P that we use for mod-switching must satisfy c*N*sigma^2*D^2/P^2<w*N,
    //    or    P > sqrt(c/w) * sigma * 2^{(L+1)*pSize/c}
    // 5. With this extra P factor, the key-switching matrices are defined
    //    relative to a modulus of size
    //          Q0 = q0*P ~ sqrt{c/w} sigma 2^{nBits*(1+1/c)}
    // 6. To get k-bit security we need N>log(Q0/sigma)(k+110)/7.2, i.e. roughly
    //          N > nBits*(1+1/c)(k+110) / 7.2

    // Compute a bound on m, and make sure that it is not too large
    double cc = 1.0+(1.0/(double)c);

    double dN = ceil((L+1)*FHE_pSize*cc*(k+110)/7.2);
    // FIXME: the bound for dN is not conservative enough...
    // this should be re-worked.
  
    long N = NTL_SP_BOUND;
    if (N > dN) N = dN;
    else {
        cerr << "Cannot support a bound of " << dN;
        Error(", aborting.\n");
    }

    long m = 0;
    size_t i=0;

    // find the first m satisfying phi(m)>=N and d | ord(p) in Z_m^*
    // and phi(m)/ord(p) >= s
    if (chosen_m) {
        if (GCD(p, chosen_m) == 1) {
        long ordP = multOrd(p, chosen_m);
        if (d == 0 || ordP % d == 0) {
            // chosen_m is OK
            m = chosen_m;
        }
        }
    }
    else if (p==2) { // use pre-computed table, divisors of 2^n-1 for some n's

        static long ms[][4] = {  // pre-computed values of [phi(m),m,d]
        //phi(m), m, ord(2),c_m*1000 (not used anymore)
        { 1176,  1247, 28,  3736}, // gens=5(42)
        { 2880,  3133, 24,  3254}, // gens=6(60), 7(!2)
        { 4050,  4051, 50, 0},     // gens=130(81)
        { 4096,  4369, 16,  3422}, // gens=129(16),3(!16)
        { 4704,  4859, 28, 0},     // gens=7(42),3(!4)
        { 5292,  5461, 14,  4160}, // gens=3(126),509(3)
        { 5760,  8435, 24,  8935}, // gens=58(60),1686(2),11(!2)
        { 7500,  7781, 50, 0},     // gens=353(30),3(!5)
        { 8190,  8191, 13,  1273}, // gens=39(630)
        { 9900, 10261, 30, 0},     // gens=3(330)
        {10752, 11441, 48,  3607}, // gens=7(112),5(!2)
        {10800, 11023, 45, 0},     // gens=270(24),2264(2),3(!5)
        {12000, 13981, 20,  2467}, // gens=10(30),23(10),3(!2)
        {11520, 15665, 24, 14916}, // gens=6(60),177(4),7(!2)
        {14112, 14351, 18, 0},     // gens=7(126),3(!4)
        {15004, 15709, 22,  3867}, // gens=5(682)
        {18000, 18631, 25,  4208}, // gens=17(120),1177(6)
        {15360, 20485, 24, 12767}, // gens=6(80),242(4),7(2)
        {16384, 21845, 16, 12798}, // gens=129(16),273(4),3(!16)
        {17280 ,21931, 24, 18387}, // gens=6(60),467(6),11(!2)
        {19200, 21607, 40, 35633}, // gens=13(120),2789(2),3(!2)
        {21168, 27305, 28, 15407}, // gens=6(126),781(6)
        {23040, 23377, 48,  5292}, // gens=35(240),5(!2)
        {23310, 23311, 45, 0},     // gens=489(518)
        {24576, 24929, 48,  5612}, // gens=12(256),5(!2)
        {27000, 32767, 15, 20021}, // gens=3(150),873(6),6945(2)
        {31104, 31609, 72,  5149}, // gens=11(216),5(!2)
        {43690, 43691, 34, 0},     // gens=69(1285)
        {49500, 49981, 30, 0},     // gens=3(1650)
        {46080, 53261, 24, 33409}, // gens=3(240),242(4),7(!2)
        {54000, 55831, 25, 0},     // gens=22(360),3529(6)
        {49140, 57337, 39,  2608}, // gens=39(630),40956(2)
        {51840, 59527, 72, 21128}, // gens=58(60),1912(6),7(!2)
        {61680, 61681, 40,  1273}, // gens=33(771),17(!2)
        {65536, 65537, 32,  1273}, // gens=2(32),3(!2048)
        {75264, 82603, 56, 36484}, // gens=3(336),24294(2),7(!2)
        {84672, 92837, 56, 38520}  // gens=18(126),1886(6),3(!2)
        };
        for (i=0; i<sizeof(ms)/sizeof(long[4]); i++) { 
            if (ms[i][0] < N || GCD(p, ms[i][1]) != 1) continue;
            long ordP = multOrd(p, ms[i][1]);
            long nSlots = ms[i][0]/ordP;
            if (d != 0 && ordP % d != 0) continue;
            if (nSlots < s) continue;

            m = ms[i][1];
            break;
        }
    }

  // If m is not set yet, just set it close to N. This may be a lousy
  // choice of m for this p, since you will get a small number of slots.

    if (m==0) {
        // search only for odd values of m, to make phi(m) a little closer to m
        for (long candidate=N|1; candidate<10*N; candidate+=2) {
            if (GCD(p,candidate)!=1) continue;

            long ordP = multOrd(p,candidate); // the multiplicative order of p mod m
            
            if (d>1 && ordP%d!=0 ) continue;
            if (ordP > 100) continue;  // order too big, we will get very few slots

            long n = phi_N(candidate); // compute phi(m)
            if (n < N) continue;       // phi(m) too small



            /*************/
            long slots = n / ordP;
            
            if (slots < dim * dim) continue;
            if (slots > dim * dim) return -1;
            /*************/
            //if (slots > dim * dim) return -1;


            cout << "phi(m): " << n << endl;
            m = candidate;  // all tests passed, return this value of m
            break;
        }
    }

    if (verbose) {
        cerr << "*** Bound N="<<N<<", choosing m="<<m <<", phi(m)="<< phi_N(m)
            << endl;
    }

    return m;
}


/* ----- some helping functions ------ */
long myModulu(long num, long mod)
//the regular modulo operator on minus numbers make trouble when trying to calculate enteries in the diagonal matrices.
//for example, -7%3=-1, but I want it to be 2 (-7+3=-4, -4+3=-1, -1+3= 2)
{
    //adding "mod" to "num" until it positive
    while(num<0)
        num+=mod;
    return num%mod;
}

unsigned int numDigit(long num){ return (num/10)==0 ? 1:numDigit(num/10)+1; } //How much digits in the input number. Uses for a "nice" printing

unsigned int largestNumInMatrixDigits(vector<vector<long> > matrix) //Find the longest number in a matrix and returns it length (in digits). Uses for a "nice" printing
{
    int largest = 0;
    long temp;
    for(unsigned int i=0; i < matrix.size(); i++)
        for(unsigned int j=0; j<matrix[i].size(); j++)
            if((temp = numDigit(matrix[i][j])) > largest)
                largest = temp;
    return largest;
}
void printNum(long num, int size){
    int len = numDigit(num);
    cout << num;
    for(int i=len; i < size; i++)
        cout << " ";
}

#define VEC_SIZE 4

// Simple class to measure time for each method
class Timer
{
public:
    void start() { m_start = my_clock(); }
    void stop() { m_stop = my_clock(); }
    double elapsed_time() const {
        return m_stop - m_start;
    }

private:
    double m_start, m_stop;
    double my_clock() const {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
};


/* --------------------- PTMatrix (Plain Text Matrix) class --------------------------*/

vector<vector<long> > PTMatrix::Row2Diagonal(const vector<vector<long> >& rowMatrix)
{
    vector<vector<long> > diagonalMatrix = vector<vector<long> >(rowMatrix[0].size(), vector<long>(rowMatrix.size(),0));
    for(unsigned int i=0, sz1 = diagonalMatrix.size(); i < sz1; i++)
        for(unsigned int j=0, sz2 = diagonalMatrix[i].size(); j < sz2; j++)
            diagonalMatrix[i][j] = rowMatrix[j][(i+j)%sz1];
    return diagonalMatrix;
}

vector<vector<long> > PTMatrix::Diagonal2Row(const vector<vector<long> >& diagonalMatrix)
{
    vector<vector<long> > rowMatrix = vector<vector<long> >(diagonalMatrix[0].size(), vector<long>(diagonalMatrix.size(),0));
    for(unsigned int i=0, sz1 = rowMatrix.size(); i < sz1; i++)
        for(unsigned int j=0, sz2 = rowMatrix[i].size(); j < sz2; j++)
            rowMatrix[i][j] = diagonalMatrix[myModulu(j-i,sz1)][i];
    return rowMatrix;
}

PTMatrix::PTMatrix(vector<vector<long> > matrix, bool isDiagonal, bool saveRowMatrix)
{   
    _d = matrix.size();
    if(isDiagonal)
    {
        _diagonalMatrix = matrix;
        if(saveRowMatrix)
            _rowMatrix = PTMatrix::Diagonal2Row(matrix);
    }
    else
    {
        _diagonalMatrix = PTMatrix::Row2Diagonal(matrix);
        if(saveRowMatrix)
            _rowMatrix = matrix;
    }
        
}

PTMatrix::PTMatrix(unsigned int dimension, int numbersLimit, bool saveRowMatrix)
/*this constructor create a random matrix
params:
 sizes: - the size of the matrix
 numbersLimit = the values limit, that means that all the values in the matrix would be between 0 to numbersLimit (not included)
*/
{
    _d = dimension;
    _diagonalMatrix = vector<vector<long> >(dimension, vector<long>(dimension));
    for(unsigned int i=0; i < dimension; i++)
        for(unsigned int j=0; j< dimension; j++)
            _diagonalMatrix[i][j] = rand() % numbersLimit;

    if(saveRowMatrix)
        _rowMatrix = PTMatrix::Diagonal2Row(_diagonalMatrix);
}


vector<vector<long> > PTMatrix::getRowMatrix() const
{
    return _rowMatrix;
}

vector<vector<long> > PTMatrix::getDiagonalMatrix() const
{
    return _diagonalMatrix;
}

void PTMatrix::print(string label) const{
    if(label.compare("")!=0)
        cout << label << endl;
    unsigned int cellSize = largestNumInMatrixDigits(_diagonalMatrix)+1; //+1 for space
    
    cout << " ";
    
    for(unsigned int i=0; i< _d*cellSize; i++)
        cout << "-";
    cout << endl;
    for(unsigned int i=0; i< _d; i++)
    {
        cout << "|";
        for(unsigned int j=0; j< _d; j++)
                printNum((*this)(i,j), cellSize);
        cout << "|" << endl;
    }
    cout << " ";
    for(unsigned int i=0; i< _d*cellSize; i++)
        cout << "-";
    cout << endl;
}

void PTMatrix::printDiagonal(string label) const
{
    if(label.compare("")!=0)
        cout << label << endl;

    unsigned int cellSize = largestNumInMatrixDigits(_diagonalMatrix)+1; //+1 for space
    
    cout << " ";
    
    for(unsigned int i=0; i< _d*cellSize; i++)
        cout << "-";
    cout << endl;
    for(unsigned int i=0; i< _d; i++)
    {
        cout << "|";
        for(unsigned int j=0; j< _d; j++)
                printNum(_diagonalMatrix[i][j], cellSize);
        cout << "|" << endl;
    }
    cout << " ";
    for(unsigned int i=0; i< _d*cellSize; i++)
        cout << "-";
    cout << endl;
}

int PTMatrix::getDimension() const
{
    return _d;
}

//operators

long& PTMatrix::operator()(unsigned int row, unsigned int column){
    if(row >= _d || column >= _d)
    {
        throw out_of_range("Error, indices out of bound!");
    }
    int i = row, j = column; //casting to int so the subtraction be ok
    return _diagonalMatrix[myModulu(j-i,_d)][row];
}

const long& PTMatrix::operator()(unsigned int row, unsigned int column) const{
    if(row >= _d || column >= _d)
    {
        throw out_of_range("Error, indices out of bound!");
    }
    int i = row, j = column; //casting to int so the subtraction be ok
    return _diagonalMatrix[myModulu(j-i,_d)][row];
}

PTMatrix PTMatrix::operator*(const PTMatrix& other) const{
    //check sizes
    if(getDimension() != other.getDimension())
        throw MatricesSizesNotMatch(getDimension(), other.getDimension());
    
    vector<vector<long> > res(getDimension(), vector<long>(getDimension(),0));
    for(unsigned int i=0; i < res.size(); i++)
        for(unsigned int j=0; j < res[i].size(); j++)
            for(unsigned int k = 0; k < other.getDimension(); k++)
                res[i][j] += (*this)(i,k)*other(k,j);
    
    return PTMatrix(res, false);
}

PTMatrix PTMatrix::operator*=(const PTMatrix& other){ return (*this) = (*this)*other; }

bool PTMatrix::isEqual(const vector<long>& other) const
{
    vector<vector<long> > rowMatrix = PTMatrix::Diagonal2Row(_diagonalMatrix);

    for(int i = 0; i < _d; i++)
        for(int j = 0; j < _d; j++)
        {
            if (rowMatrix[i][j] != other[i * _d + j])
                return false;
        }

    return true; 
}

vector<long> PTMatrix::MatrixEncoding(const vector<vector<long> >& rowMatrix)
{
    vector<long> v;
    for(int i=0; i < rowMatrix.size(); i++)
    {
        v.insert(v.end(), rowMatrix[i].begin(), rowMatrix[i].end());
    }
    return v;
}
EncryptedMatrix PTMatrix::encrypt(const EncryptedArray& ea, const FHEPubKey& publicKey, bool saveRow) const{
    vector<Ctxt> encDiagonalMatrix(_d, Ctxt(publicKey));
    Ctxt encRowMatrix(publicKey);

    unsigned int nslots = ea.size();
    for(unsigned int i=0; i< _d; i++)
        {
            vector<long> temp = _diagonalMatrix[i];
            temp.resize(nslots,0);
            ea.encrypt(encDiagonalMatrix[i], publicKey, temp);
        }
    
    if (saveRow)
    {
        vector<long> flatMatrix = PTMatrix::MatrixEncoding(_rowMatrix);
        flatMatrix.resize(nslots, 0);
        ea.encrypt(encRowMatrix, publicKey, flatMatrix);
    }

    return EncryptedMatrix(encDiagonalMatrix, encRowMatrix, _d, saveRow);
}

EncryptedMatrix PTMatrix::encrypt(const FHEPubKey& publicKey, bool saveRow) const{
    EncryptedArray ea(publicKey.getContext());
    return encrypt(ea, publicKey, saveRow);
}

PTMatrix PTMatrix::sigmaPermutation(unsigned int d){
    long dSquare = d * d;
    vector<vector<long> > rowMatrix(dSquare, vector<long>(dSquare, 0));
    for(int l=0; l < dSquare; l++)
        for(int i=0; i < d; i++)
            for(int j=0; j < d; j++)
            {
                if(l == d*i + myModulu(i+j,d)) rowMatrix[d*i+j][l] = 1;
            }
    return PTMatrix(rowMatrix, false, false);
}
PTMatrix PTMatrix::tauPermutation(unsigned int d){
    long dSquare = d * d;
    vector<vector<long> > rowMatrix(dSquare, vector<long>(dSquare, 0));
    for(int l=0; l < dSquare; l++)
        for(int i=0; i < d; i++)
            for(int j=0; j < d; j++)
            {
                if(l == d*myModulu(i+j,d) + j) rowMatrix[d*i+j][l] = 1;
            }
    return PTMatrix(rowMatrix, false, false);
}
PTMatrix PTMatrix::phiPermutation(unsigned int d, int k){
    long dSquare = d * d;
    vector<vector<long> > rowMatrix(dSquare, vector<long>(dSquare, 0));
    for(int l=0; l < dSquare; l++)
        for(int i=0; i < d; i++)
            for(int j=0; j < d; j++)
            {
                if(l == d*i + myModulu(k+j,d)) rowMatrix[d*i+j][l] = 1;
            }
    return PTMatrix(rowMatrix, false, false);
}
PTMatrix PTMatrix::psiPermutation(unsigned int d, int k){
    long dSquare = d * d;
    vector<vector<long> > rowMatrix(dSquare, vector<long>(dSquare, 0));
    for(int l=0; l < dSquare; l++)
        for(int i=0; i < d; i++)
            for(int j=0; j < d; j++)
            {
                if(l == d*myModulu(i+k,d) + j) rowMatrix[d*i+j][l] = 1;
            }
    return PTMatrix(rowMatrix, false, false);
}

vector<ZZX> PTMatrix::DiagonalEncoding(const EncryptedArray& ea)
{
    vector<ZZX> zzxVec(_d);
    for(int i = 0; i < _d; i++)
    {
        vector<long> temp = _diagonalMatrix[i];
        temp.resize(ea.size());
        ea.encode(zzxVec[i], temp);
    }
    return zzxVec;
}

ZZX VectorEncoding(const vector<long>& v, const EncryptedArray& ea)
{
    ZZX zzxVec;
    vector<long> temp = v;
    temp.resize(ea.size());
    ea.encode(zzxVec, temp);
    return zzxVec;
}


/* --------------------- EncryptedMatrix class -------------*/
EncryptedMatrix::EncryptedMatrix(
                 const vector<Ctxt>& encDiagonalMatrix,
                 const Ctxt& encRowMatrix,
                 int dimension, bool haveRowMatrix): 
                 _rowMatrix(encRowMatrix),
                 _diagonalMatrix(encDiagonalMatrix),
                 _d(dimension),
                 _haveRowMatrix(haveRowMatrix){}

PTMatrix EncryptedMatrix::decrypt(const EncryptedArray& ea, const FHESecKey& secretKey) const {
    vector<vector<long> > diagonalMatrix(_d);
    for(unsigned int i=0; i < _d; i++){
        ea.decrypt(_diagonalMatrix[i], secretKey, diagonalMatrix[i]);
        diagonalMatrix[i].resize(_d, 0);
    }
    return PTMatrix(diagonalMatrix, true, false);
}

PTMatrix EncryptedMatrix::decrypt(const FHESecKey& secretKey) const {
    EncryptedArray ea(secretKey.getContext());
    return decrypt(ea, secretKey);
}

//matrix multyplication by vector. NOTE: this return a column vector! so don't use it to create a matrix (unless you want it to be column vectors matrix)
Ctxt EncryptedMatrix::operator*(const Ctxt& vec) const{
    EncryptedArray ea(vec.getContext());
    Ctxt result(vec.getPubKey());
    int len = _diagonalMatrix.size();
    
    //TODO: Still not perfectlly working
    
    Ctxt fixedVec = vec;
    cout << "ea.size " << ea.size() << " matrix row size " << _d << endl;
    if(ea.size() != _d) //Fix the problem that if the size of the vector is not nslots, the zero padding make the rotation push zeros to the begining of the vector
    {
        //replicate the vector to fill instead of zero padding
        for(unsigned int length = _d; length < ea.size(); length*=2){
            Ctxt copyVec = fixedVec;
            ea.shift(copyVec, length);  //shift length to right
            fixedVec+=copyVec;
        }
    }
    
    for(int i=0; i < len; i++)
    {
        Ctxt rotatedVec(fixedVec);   //copy vec
        ea.rotate(rotatedVec, -i);   //rotate it i right (-i left)
        rotatedVec *= _diagonalMatrix[i];
        result += rotatedVec;
    }
    return result;
}

Ctxt EncryptedMatrix::LinTrans1(const vector<vector<long> >& matrix) const{
    EncryptedArray ea(_diagonalMatrix[0].getContext());
    Ctxt result(_diagonalMatrix[0].getPubKey());
    int len = matrix.size();
    
    //TODO: Still not perfectlly working
    Ctxt fixedVec = _rowMatrix;
    if(ea.size() != len) //Fix the problem that if the size of the vector is not nslots, the zero padding make the rotation push zeros to the begining of the vector
    {
        //replicate the vector to fill instead of zero padding
        for(unsigned int length =len; length < ea.size(); length*=2){
            Ctxt copyVec = fixedVec;
            ea.shift(copyVec, length);  //shift length to right
            fixedVec+=copyVec;
        }
    }

    int sqrtD = (int)sqrt((double)_d);
    for(int i = -sqrtD + 1; i < sqrtD; i++)
    {
        Ctxt temp(_diagonalMatrix[0].getPubKey());
        for(int j = 0; j < sqrtD; j++)
        {
            Ctxt rotatedVec(fixedVec);   //copy vec
            ea.rotate(rotatedVec, -j);   //rotate it i right (-i left)
            vector<long> rotatedVec2 = matrix[myModulu(sqrtD*i+j, len)];
            int step = sqrtD * i;
            if (step > 0)
                // right rotation with abs(step) steps
                std::rotate(rotatedVec2.begin(), rotatedVec2.begin()+rotatedVec2.size()-step, rotatedVec2.end());
            else
            {
                // left rotation with abs(step) steps
                std::rotate(rotatedVec2.begin(), rotatedVec2.begin()-step, rotatedVec2.end());
            }
            
            ZZX U;
            if(ea.size() != rotatedVec2.size())
                rotatedVec2.resize(ea.size());
                
            ea.encode(U, rotatedVec2);
            rotatedVec.multByConstant(U);
            temp += rotatedVec;
        }
        ea.rotate(temp, -(sqrtD * i));
        result += temp;
    }
    return result;
}

Ctxt EncryptedMatrix::LinTrans2(const vector<ZZX>& matrix) const{
    EncryptedArray ea(_diagonalMatrix[0].getContext());
    Ctxt result(_diagonalMatrix[0].getPubKey());
    int len = matrix.size();
    
    //TODO: Still not perfectlly working
    
    Ctxt fixedVec = _rowMatrix;
    if(ea.size() != len) //Fix the problem that if the size of the vector is not nslots, the zero padding make the rotation push zeros to the begining of the vector
    {
        //replicate the vector to fill instead of zero padding
        for(unsigned int length =len; length < ea.size(); length*=2){
            Ctxt copyVec = fixedVec;
            ea.shift(copyVec, length);  //shift length to right
            fixedVec+=copyVec;
        }
    }
    
    for(int i=0; i < _d; i++)
    {
        Ctxt rotatedVec(fixedVec);   //copy vec
        ea.rotate(rotatedVec, -(_d*i));   //rotate it i right (-i left)
        rotatedVec.multByConstant(matrix[_d*i]);
        result += rotatedVec;
    }

    return result;
}

Ctxt EncryptedMatrix::LinTrans3(const Ctxt& vec, const vector<ZZX>& matrix, int d, int k){
    EncryptedArray ea(vec.getContext());
    Ctxt result(vec.getPubKey());
    int len = matrix.size();
    
    //TODO: Still not perfectlly working
    
    Ctxt fixedVec = vec;
    if(ea.size() != len) //Fix the problem that if the size of the vector is not nslots, the zero padding make the rotation push zeros to the begining of the vector
    {
        //replicate the vector to fill instead of zero padding
        for(unsigned int length =len; length < ea.size(); length*=2){
            Ctxt copyVec = fixedVec;
            ea.shift(copyVec, length);  //shift length to right
            fixedVec+=copyVec;
        }
    }
    

    Ctxt rotatedVec(fixedVec);   //copy vec
    ea.rotate(rotatedVec, -k);   //rotate it i right (-i left)
    rotatedVec.multByConstant(matrix[k]);
    result = rotatedVec;


    rotatedVec =fixedVec;   //copy vec
    ea.rotate(rotatedVec, d-k);   //rotate it i right (-i left)
    rotatedVec.multByConstant(matrix[myModulu(k-d, len)]);

    result += rotatedVec;
    return result;
}

Ctxt EncryptedMatrix::LinTrans4(const Ctxt& vec, int d, int k){
    EncryptedArray ea(vec.getContext());
    Ctxt result(vec.getPubKey());
    int len = d * d;
    
    //TODO: Still not perfectlly working
    
    Ctxt fixedVec = vec;
    if(ea.size() != len) //Fix the problem that if the size of the vector is not nslots, the zero padding make the rotation push zeros to the begining of the vector
    {
        //replicate the vector to fill instead of zero padding
        for(unsigned int length =len; length < ea.size(); length*=2){
            Ctxt copyVec = fixedVec;
            ea.shift(copyVec, length);  //shift length to right
            fixedVec+=copyVec;
        }
    }
    

    Ctxt rotatedVec(fixedVec);   //copy vec
    ea.rotate(rotatedVec, -(d*k));   //rotate it i right (-i left)
    result = rotatedVec;

    return result;
}

int EncryptedMatrix::getD()
{
    return _d;
}

Ctxt EncryptedMatrix::getRowMatrix()
{
    return _rowMatrix;
}

Ctxt EncryptedMatrix::operator*( EncryptedMatrix& other) 
{
    //check sizes
    if(_d != other.getD())
        throw MatricesSizesNotMatch(_d, other.getD());

    Ctxt vec = _diagonalMatrix[0]; //save it for faster vec.getPubKey() in the loop
    EncryptedArray ea(vec.getContext());

    vector<Ctxt> A(_d, Ctxt(vec.getPubKey())), B(_d, Ctxt(vec.getPubKey()));

    Timer a0Init;
    a0Init.start();
    A[0] = LinTrans1(PTMatrix::sigmaPermutation(_d).getDiagonalMatrix());
    a0Init.stop();
    std::cout << "Time taken for the a[0]: " << a0Init.elapsed_time() << std::endl;   

    Timer b0Init;
    b0Init.start();
    B[0] = other.LinTrans2(PTMatrix::tauPermutation(_d).DiagonalEncoding(ea));
    b0Init.stop();
    std::cout << "Time taken for the b[0]: " << b0Init.elapsed_time() << std::endl;

    for(int k = 1; k < _d; k++)
    {
        Timer akInit;
        akInit.start();
        A[k] = LinTrans3(A[0], PTMatrix::phiPermutation(_d, k).DiagonalEncoding(ea), _d, k);
        akInit.stop();
        std::cout << "Time taken for the a[" << k <<"]: " << akInit.elapsed_time() << std::endl;

        Timer bkInit;
        bkInit.start();
        B[k] = LinTrans4(B[0], _d, k);
        bkInit.stop();
        std::cout << "Time taken for the b[" << k <<"]: " << bkInit.elapsed_time() << std::endl;
    }

    Timer multiplySum;
    multiplySum.start();
    Ctxt result = A[0];
    result.multiplyBy(B[0]) ;
    for(int k = 1; k < _d; k++)
    {
        Ctxt temp = A[k];
        temp *= B[k];
        result += temp;
    }
    multiplySum.stop();
    std::cout << "Time taken for multiplySum " << multiplySum.elapsed_time() << std::endl;
    return result;
}

void test(int dimension, long p)
{
    long m = 0;                   // Specific modulus
	//long p = 113;//16487    32003             // Plaintext base [default=2], should be a prime number
	long r = 1;                   // Lifting [default=1]
	long L = 450;                  // Number of levels in the modulus chain [default=heuristic]
	long c = 3;                   // Number of columns in key-switching matrix [default=2]
	long w = 64;                  // Hamming weight of secret key
	long d = 1;                   // Degree of the field extension [default=1]
	long k = 80;                  // Security parameter [default=80] 
    long s = 0;                   // Minimum number of slots [default=0]

    //int dimension = 8; 

    Timer tInit;
	tInit.start();
    cout << "current p is: " << p << endl;
    m = FindM1(k, L, c, p, d, s, 0, dimension, false);           // Find a value for m given the specified values
    if (m == -1) return;

    //m=32003-1;//1907 is a safe prime
    // m=16487-1; 
    std::cout << "Initializing context... " << std::flush;
	FHEcontext context(m, p, r); 	          // Initialize context
	buildModChain(context, L, c);             // Modify the context, adding primes to the modulus chain
	std::cout << "OK!" << std::endl;
    cout<<"securitylevel="<<context.securityLevel()<<endl;

    //return 0;
	std::cout << "Generating keys... " << std::flush;
	FHESecKey secretKey(context);                    // Construct a secret key structure
	const FHEPubKey& publicKey = secretKey;                 // An "upcast": FHESecKey is a subclass of FHEPubKey
	secretKey.GenSecKey(w);                          // Actually generate a secret key with Hamming weight
	addSome1DMatrices(secretKey);                    // Extra information for relinearization

    ZZX G =  context.alMod.getFactorsOverZZ()[0]; 
    EncryptedArray ea(context, G);

    tInit.stop();
	std::cout << "Time taken for the initialization: " << tInit.elapsed_time() << std::endl;
    cout << "m: " << m << endl;
    cout << "nslots: " << ea.size() << endl;   

    Timer ptMatrixInit;
	ptMatrixInit.start();
    PTMatrix ptMatrix1(dimension,3, true);
    EncryptedMatrix encMatrix1 = ptMatrix1.encrypt(publicKey, true);

    
    PTMatrix ptMatrix2(dimension, 3, true);
    EncryptedMatrix encMatrix2 = ptMatrix2.encrypt(publicKey, true);
    ptMatrixInit.stop();
    std::cout << "Time taken for the ptMatrix initialization: " << ptMatrixInit.elapsed_time() << std::endl;
    ptMatrix1.print();
    ptMatrix2.print();

    
    PTMatrix ptResult = ptMatrix1 * ptMatrix2;
    cout << "matrix multiplication result:" << endl;
    ptResult.print();

    Timer totalTime;
	totalTime.start();
    Ctxt result = encMatrix1 * encMatrix2;
    totalTime.stop();
	std::cout << "matrix mul total time " << totalTime.elapsed_time() << std::endl;

    vector<long> temp(ea.size(), 0);
    ea.decrypt(result, secretKey, temp);
    for(int i = 0; i < dimension * dimension; i++)
        cout << temp[i] << ' ';
    cout << endl << endl << endl << endl;

    if (ptResult.isEqual(temp))
        cout << "Correct decryption!" << endl;

}
int main()
{
    ifstream file("prime.txt");
    string line;
    long prime;
    int dimension;
    cout << "dimension: ";
    cin >> dimension;
    if (file.is_open())
    {
        while(!file.eof())
        {
            getline(file, line); //get number of rows 
            prime = stoi(line);
            test(dimension, prime);

        }
        file.close();
    }
    else
    {
        cout << "file not exist" << endl;
        return 1;
    }
    
    return 0;
}

