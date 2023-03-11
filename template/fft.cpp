class FFT {
public:
    vector<complex<double>> v;
    int n;
    vector<complex<double> > fft() {
        if (n == 1) {
            return v;
        }
        vector<complex<double>> v1(n / 2), v2(n / 2);
        for (int i = 0; i < n; i += 2) {
            v1[i / 2] = v[i];
            v2[i / 2] = v[i + 1];
        }
        FFT f1(v1, n / 2), f2(v2, n / 2);
        vector<complex<double> > r1 = f1.fft();
        vector<complex<double> > r2 = f2.fft();
        vector<complex<double> > r(n);
        for (int i = 0; i < n; i++) {
            complex<double> w(cos(-2 * M_PI * i / n), sin(-2 * M_PI * i / n));
            r[i] = r1[i % (n / 2)] + w * r2[i % (n / 2)];
        }
        return r;
    }   
    vector<complex<double> > nfft() { 
        vector<complex<double> > r = _nfft();
        for (int i = 0; i < n; i++) {
            r[i] /= n;
        }
        return r;
    }
    vector<complex<double> > _nfft() {
        if (n == 1) {
            return v;
        }
        vector<complex<double>> v1(n / 2), v2(n / 2);
        for (int i = 0; i < n; i += 2) {
            v1[i / 2] = v[i];
            v2[i / 2] = v[i + 1];
        }
        FFT f1(v1, n / 2), f2(v2, n / 2);
        vector<complex<double> > r1 = f1._nfft();
        vector<complex<double> > r2 = f2._nfft();
        vector<complex<double> > r(n);
        for (int i = 0; i < n; i++) {
            complex<double> w(cos(2 * M_PI * i / n), sin(2 * M_PI * i / n));
            r[i] = r1[i % (n / 2)] + w * r2[i % (n / 2)];
            // r[i] /= n;
        }
        return r;
    }

    FFT(vector<complex<double> > t, int n) {
        this->n = n;
        v.resize(n);
        for (int i = 0; i < t.size(); i++) {
            v[i] = t[i];
        }
    }    
    FFT(vector<int > t, int n) {
        this->n = n;
        v.resize(n);
        for (int i = 0; i < t.size(); i++) {
            v[i] = t[i];
        }
    }
};



class NTT {
public:
    vector<mint> v;
    int n;
    // 998244353
    mint w, w_inv;
    NTT(const vector<mint> &t, int n) {
        this->n = n;
        assert((998244353 - 1) % n == 0);
        w = mint(3).power((998244353 - 1) / n);
        w_inv = w.inv();
        v.resize(n);
        for (int i = 0; i < t.size(); i++) {
            v[i] = t[i];
        }
    }
    void ntt_inplace() {
        int len = 0;
        while ((1 << len) < n) {
            len++;
        }
        vector<int> rev(n);
        for(int i=0;i<n;i++) {
            rev[i] = rev[i>>1] >> 1 | (i&1) << (len-1);
        }
        for(int i = 0;i < n; i++) {
            if(i<rev[i]) {
                swap(v[i],v[rev[i]]);
            }
        }
        for(int i = 1; i < n; i<<=1){
            mint wn = mint(3).inv().power((mod - 1) / (i << 1));
            for(int p=i<<1,j=0;j<n;j+=p){
                mint x=1;
                for(int k=0;k<i;k++,x = x * wn){
                    mint u=v[j+k], t = x*v[i+j+k];
                    v[j+k] = u + t;
                    v[i+j+k] = u - t;
                }
            }
        }
    }
    void intt_inplace() {
        int len = 0;
        while ((1 << len) < n) {
            len++;
        }
        vector<int> rev(n);
        for(int i=0;i<n;i++) {
            rev[i] = rev[i>>1] >> 1 | (i&1) << (len-1);
        }
        for(int i = 0;i < n; i++) {
            if(i<rev[i]) {
                swap(v[i],v[rev[i]]);
            }
        }
        for(int i = 1; i < n; i<<=1){
            mint wn = mint(3).power((mod - 1) / (i << 1));
            for(int p=i<<1,j=0;j<n;j+=p){
                mint x=1;
                for(int k=0;k<i;k++,x = x * wn){
                    mint u=v[j+k], t = x*v[i+j+k];
                    v[j+k] = u + t;
                    v[i+j+k] = u - t;
                }
            }
        }
        mint mul = mint(n).inv();
        for (int i = 0; i < n; i++) {
            v[i] *= mul;
        }
    }
    
    vector<mint> ntt() {
        if (n == 1) {
            return v;
        }
        vector<mint> v1(n / 2), v2(n / 2);
        for (int i = 0; i < n; i += 2) {
            v1[i / 2] = v[i];
            v2[i / 2] = v[i + 1];
        }
        NTT f1(v1, n / 2), f2(v2, n / 2);
        vector<mint > r1 = f1.ntt();
        vector<mint > r2 = f2.ntt();
        vector<mint > r(n);
        mint ww = 1;
        for (int i = 0; i < n; i++) {
            r[i] = r1[i % (n / 2)] + ww * r2[i % (n / 2)];
            ww *= w_inv;
        }
        return r;
    }
    vector<mint > intt() { 
        vector<mint > r = _intt();
        for (int i = 0; i < n; i++) {
            r[i] /= n;
        }
        return r;
    }
    vector<mint> _intt() {
        if (n == 1) {
            return v;
        }
        vector<mint> v1(n / 2), v2(n / 2);
        for (int i = 0; i < n; i += 2) {
            v1[i / 2] = v[i];
            v2[i / 2] = v[i + 1];
        }
        NTT f1(v1, n / 2), f2(v2, n / 2);
        vector<mint > r1 = f1._intt();
        vector<mint > r2 = f2._intt();
        vector<mint > r(n);
        mint ww = 1;
        for (int i = 0; i < n; i++) {
            r[i] = r1[i % (n / 2)] + ww * r2[i % (n / 2)];
            ww *= w;
        }
        return r;
    }
};

class Poly {
public:
    vector<mint> v;
    int deg = 0;
    Poly(vector<mint> v, int deg=-1) {
        this->v = v;
        if (deg == -1) {
            this->deg = v.size() - 1;
        } else {
            this->deg = deg;
        }
    }
    void show() {
        for (int i = 0; i <= deg; i++) {
            cerr << v[i] << " \n"[i == deg];
        }
    }
    // bool operator<(const Poly& b) {
    //     return this->deg >= b.deg;
    // } 
};

Poly mul(const Poly &x, const Poly &y) {
    int d = x.deg + y.deg;
    int n = 1;
    while (n < d + 1) {
        n *= 2;
    }
    NTT nx = NTT(x.v, n);
    NTT ny = NTT(y.v, n);
    nx.ntt_inplace();
    ny.ntt_inplace();
    // vector<mint> xx = nx.ntt();
    // vector<mint> yy = ny.ntt();
    for (int i = 0; i < nx.v.size(); i++) {
        nx.v[i] *= ny.v[i];
    }
    nx.intt_inplace();
    return Poly(nx.v, d);
}

int my_cmp(const Poly &x, const Poly &y) {return x.deg > y.deg;}
