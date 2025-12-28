#include <algorithm>
#include <array>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstring>
#include <filesystem>
#include <numeric>
#include <cassert>
#include <iomanip>

using namespace std;

template<class T>
static void write_le(ostream &os, T v){
    for(size_t i = 0; i<sizeof(T); i++) os.put((char)((v >> (8*i)) & 0xFF));
}

template<class T>
static T read_le(istream &is) {
    T v = 0; 
    for(size_t i = 0; i<sizeof(T); i++) {
        int c = is.get();
        if(c == EOF) throw runtime_error("EOF");
        v |= ((uint64_t)(unsigned char)c)<<(8*i);
    }
    return v;
}

// Fenwick tre for kumulative frekvenser 
struct Fenwick {
    int n;
    vector<uint32_t> bit;
    uint64_t total;
    Fenwick(int n_ = 257): n(n_), bit(n_+1 ,0), total(0) {}
        
    void add(int idx, uint32_t delta) {
        for(int i = idx + 1; i <= n; i += i & -i) {
            bit[i] += delta;    
        }
        total += delta;
    }

    uint64_t sumPrefix(int idx) const {
        uint64_t s = 0; 
        for(int i = idx + 1; i > 0; i -= i& - i) s += bit[i];
        return s;
    }

    int findByCumulative(uint32_t off) const {
        int lo = -1, hi = n - 1;
        while (hi - lo > 1) {
            int mid = lo + ((hi - lo) >> 1);
            if (sumPrefix(mid) <= off) lo = mid;
            else hi = mid;
        }
        return hi;
    }
};

// Range Encoder/Decoder (lo, hi)
static constexpr uint64_t TOP = (1ull << 56);

struct RangeEncoder {
    uint64_t low   = 0;
    uint32_t range = 0xFFFFFFFFu;
    uint8_t  cache = 0;
    uint32_t ff_count = 0;           
    bool     started = false;
    ostream &out;

    explicit RangeEncoder(ostream &os) : out(os) {}

    inline void shift_low() {
        uint32_t top = (uint32_t)(low >> 24);
        if (!started) {
            cache = (uint8_t)top;
            started = true;
        } else if ( ((uint32_t)low < 0xFF000000u) || (low >> 32) ) {
            uint8_t carry = (uint8_t)(low >> 32);
            out.put((char)(cache + carry));
            uint8_t fill = carry ? 0x00 : 0xFF;
            while (ff_count--) out.put((char)fill);
            cache = (uint8_t)top;
            ff_count = 0;
        } else {
            ++ff_count;
        }
        low = ((uint32_t)low) << 8;
    }

    inline void normalize() {
        while (range < 0x01000000u) {
            range <<= 8;
            shift_low();
        }
    }

    inline void encode(uint32_t lo, uint32_t hi, uint32_t total) {
        uint32_t freq = hi - lo;
        if (!freq) { ++hi; ++freq; }

        // quotient-based split
        uint32_t q = range / total;
        low   += (uint64_t)q * lo;
        range  = (uint32_t)((uint64_t)q * freq);

        normalize();
    }


    inline void finish() {
        low += 0x01000000u;
        for (int i = 0; i < 4; ++i) shift_low();
    }
};



struct RangeDecoder {
    uint32_t code  = 0;
    uint32_t range = 0xFFFFFFFFu;
    istream &in;

    size_t  budget     = 0;
    size_t  bytes_read = 0;

    RangeDecoder(istream &is, size_t budget_bytes) : in(is), budget(budget_bytes) { init(); }
    explicit RangeDecoder(istream &is) : in(is), budget(SIZE_MAX) { init(); }

    inline unsigned char get_byte_or_zero() {
        if (bytes_read >= budget) return 0;
        int c = in.get();
        if (c == EOF) return 0;
        ++bytes_read;
        return (unsigned char)c;
    }

    inline void init() {
        code = 0;
        for (int i = 0; i < 4; ++i) code = (code << 8) | get_byte_or_zero();
    }

    inline void normalize() {
        while (range < 0x01000000u) {
            code  = (code << 8) | get_byte_or_zero();
            range <<= 8;
        }
    }

    inline uint32_t get_count(uint32_t total) const {
        uint32_t q = range / total;
        uint32_t t = code / q;
        if (t >= total) t = total - 1;
        return t;
    }

    inline void remove_symbol(uint32_t lo, uint32_t hi, uint32_t total) {
        uint32_t q = range / total;
        code  -= (uint32_t)((uint64_t)q * lo);
        range  = (uint32_t)((uint64_t)q * (hi - lo));
        normalize();
    }

};

struct RangeModel {
    Fenwick fw;
    std::vector<uint32_t> freq;
    static constexpr int ALPH=257;
    static constexpr int MATCH_SYM = 256;
    static constexpr uint32_t MAX_TOTAL = 1u << 20;
    uint64_t init_total;
    RangeModel(): fw(ALPH), freq(ALPH, 1), init_total(0) {
        for(int s = 0; s < ALPH; s++) fw.add(s,1);
        init_total = fw.total;
    }

    void rescale_if_needed() {
        if (fw.total <= MAX_TOTAL) return;
        for (int s=0; s<ALPH; ++s) {
            uint32_t f = (freq[s] + 1) >> 1;
            if (!f) f = 1;
            freq[s] = f;
        }
        fw = Fenwick(ALPH);
        for (int s=0; s<ALPH; ++s) fw.add(s, freq[s]);
    }

    void encodeSymbol(RangeEncoder& enc, int sym) {
        uint32_t cumLow  = (sym == 0) ? 0u : (uint32_t)fw.sumPrefix(sym - 1);
        uint32_t cumHigh = (uint32_t)fw.sumPrefix(sym);
        uint32_t total   = (uint32_t)fw.total;
        assert(cumLow < cumHigh && cumHigh <= total);

        enc.encode(cumLow, cumHigh, total);

        fw.add(sym, 1);
        ++freq[sym];
        rescale_if_needed();
    }

    int decodeSymbol(RangeDecoder& dec) {
        uint32_t total = (uint32_t)fw.total;
        uint32_t off   = dec.get_count(total);

        int sym = fw.findByCumulative(off);

        uint32_t cumLow  = (sym == 0) ? 0u : (uint32_t)fw.sumPrefix(sym - 1);
        uint32_t cumHigh = (uint32_t)fw.sumPrefix(sym);
        assert(cumLow <= off && off < cumHigh);

        dec.remove_symbol(cumLow, cumHigh, total);

        fw.add(sym, 1);
        ++freq[sym];
        rescale_if_needed();
        return sym;
    }

};

static void rc_pack_file(const string& inPath, const string& outPath) {
    ifstream is(inPath, ios::binary);
    if (!is) throw runtime_error("can't open input");
    vector<uint8_t> data{istreambuf_iterator<char>(is), {}};

    ofstream os(outPath, ios::binary);
    if (!os) throw runtime_error("can't open output");

    os.write("RC01", 4);
    write_le<uint64_t>(os, (uint64_t)data.size());

    ostringstream payload(ios::out|ios::binary);
    RangeEncoder enc(payload);
    const uint32_t TOTAL = 257;
    for (uint8_t b : data) {
        uint32_t lo = (uint32_t)b;
        uint32_t hi = lo + 1;
        enc.encode(lo, hi, TOTAL);
    }
    enc.finish();
    const string pay = payload.str();
    write_le<uint32_t>(os, (uint32_t)pay.size());
    os.write(pay.data(), pay.size());
}

static void rc_unpack_file(const string& inPath, const string& outPath) {
    ifstream is(inPath, ios::binary);
    if (!is) throw runtime_error("can't open input");

    char magic[4]; is.read(magic, 4);
    if (is.gcount()!=4 || strncmp(magic, "RC01", 4)!=0) throw runtime_error("bad magic");
    uint64_t n = read_le<uint64_t>(is);
    uint32_t cbytes = read_le<uint32_t>(is);
    string payload(cbytes, '\0');
    is.read(payload.data(), cbytes);
    if ((uint32_t)is.gcount()!=cbytes) throw runtime_error("truncated payload");

    istringstream pis(payload, ios::in|ios::binary);
    RangeDecoder dec(pis, cbytes);
    const uint32_t TOTAL = 257;

    vector<uint8_t> out; out.reserve(n);
    for (uint64_t i=0;i<n;++i) {
        uint32_t off = dec.get_count(TOTAL);
        dec.remove_symbol(off, off+1, TOTAL);
        out.push_back((uint8_t)off);
    }
    if (dec.bytes_read != cbytes) throw runtime_error("range payload not fully consumed");

    ofstream os(outPath, ios::binary);
    if (!os) throw runtime_error("can't open output");
    os.write((char*)out.data(), out.size());
}

// BWT & bruk av sirkular suffiks array
static vector<int> build_sa(const vector<uint8_t>& s) {
    const int n = (int)s.size();
    if (n == 0) return {};
    vector<int> sa(n), r(n), nr(n), key(n);

    for (int i = 0; i < n; ++i) { sa[i] = i; r[i] = (int)s[i]; }

    auto counting_sort = [&](const vector<int>& key, int max_key, vector<int>& idx) {
        vector<int> cnt(max(256, max_key + 2), 0);
        for (int i = 0; i < n; ++i) ++cnt[key[idx[i]] + 1];
        for (int i = 1; i < (int)cnt.size(); ++i) cnt[i] += cnt[i - 1];
        vector<int> out(n);
        for (int i = 0; i < n; ++i) out[cnt[key[idx[i]]]++] = idx[i];
        idx.swap(out);
    };

    for (int k = 1; k < n; k <<= 1) {
        for (int pos = 0; pos < n; ++pos) key[pos] = r[(pos + k) % n];
        counting_sort(key, *max_element(key.begin(), key.end()), sa);

        for (int pos = 0; pos < n; ++pos) key[pos] = r[pos];
        counting_sort(key, *max_element(key.begin(), key.end()), sa);

        // re-rank
        nr[sa[0]] = 0;
        int distinct = 1;
        for (int i = 1; i < n; ++i) {
            int a = sa[i - 1], b = sa[i];
            if (r[a] != r[b] || r[(a + k) % n] != r[(b + k) % n]) ++distinct;
            nr[b] = distinct - 1;
        }
        r.swap(nr);
        if (distinct == n) break; // all unique
    }
    return sa;
}

static void bwt_transform(const vector<uint8_t>& src, vector<uint8_t>& dst, uint32_t& primary) {
    const int n = (int)src.size();
    dst.resize(n);
    auto sa = build_sa(src); // cyclic order
    primary = 0;
    for (int i = 0; i < n; ++i) {
        int j = sa[i];
        dst[i] = src[(j + n - 1) % n];
        if (j == 0) primary = i;
    }
}

static void bwt_inverse(const vector<uint8_t>& bwt, uint32_t primary, vector<uint8_t>& out) {
    const int n = (int)bwt.size();
    out.resize(n);
    if (primary >= (uint32_t)n) throw runtime_error("BWT primary out of range");

    array<uint32_t, 256> cnt{}; 
    for (uint8_t c : bwt) cnt[c]++;

    array<uint32_t, 256> first{};
    int sum = 0;
    for (int c = 0; c < 256; ++c) { first[c] = sum; sum += cnt[c]; }

    // LF mapping
    vector<uint32_t> tally(256, 0), next(n);
    for (int i = 0; i < n; ++i) {
        uint8_t c = bwt[i];
        next[i] = first[c] + tally[c]++;
    }

    int j = (int)primary;
    for (int i = n - 1; i >= 0; --i) {
        out[i] = bwt[j];
        j = next[j];
    }
}

// LZP
struct LZP {
    static constexpr uint32_t HASH_SIZE = 1<<20;
    static constexpr uint32_t HASH_MASK = HASH_SIZE - 1;
    vector<uint8_t> table;
    LZP(): table(HASH_SIZE, 0) {}
    static inline uint32_t hash4(uint8_t a, uint8_t b, uint8_t c, uint8_t d) {
        uint32_t x = ((uint32_t)a<<24)|((uint32_t)b<<16)|((uint32_t)c<<8)|d;
        // MurmurHash3 konstanter
        x ^= x>>17;
        x *= 0x85ebca6b;
        x ^= x>>13;
        x *= 0xc2b2ae35;
        x ^= x>>16;
        return x & HASH_MASK;
    }

    //encode
    template<class Emit>
    void encode(const vector<uint8_t>& data, Emit emit) {
        uint8_t c1 = 0, c2 = 0, c3 = 0, c4 = 0;
        for(size_t i = 0; i < data.size(); i++) {
            uint8_t cur = data[i];
            uint32_t h = hash4(c1,c2,c3,c4);
            uint8_t pred = table[h];
            if (pred==cur) { emit(RangeModel::MATCH_SYM); }
            else { emit((int)cur); table[h] = cur; }
            c1 = c2; c2 = c3; c3 = c4; c4 = cur;
        }
    }

    //decode
    template<class NextSym>
    void decode(size_t n, NextSym nextSym, vector<uint8_t>& out) {
        out.clear();
        out.reserve(n);
        uint8_t c1 = 0, c2 = 0, c3 = 0, c4 = 0;
        for(size_t i = 0; i<n; i++) {
            int sym = nextSym();
            if(sym == RangeModel::MATCH_SYM) {
                uint32_t h = hash4(c1,c2,c3,c4);
                uint8_t b = table[h]; 
                out.push_back(b); 
                c1=c2; c2=c3; c3=c4; c4=b;
            } else {
                uint8_t b = (uint8_t)sym; 
                uint32_t h = hash4(c1,c2,c3,c4);
                table[h] = b;
                out.push_back(b);
                c1=c2; c2=c3; c3=c4; c4=b;
            }
        }
    }
};

// Kompressor/Dekompressor
static constexpr uint32_t DEFAULT_BLOCK = 1u<<20;
void compress_file(const string& inPath, const string& outPath, uint32_t blockSize=DEFAULT_BLOCK) {
    ifstream is(inPath, ios::binary);
    if (!is) throw runtime_error("Kan ikke åpne input");

    ofstream os(outPath, ios::binary);
    if (!os) throw runtime_error("Kan ikke åpne output");

    // Header
    os.write("BSC1", 4);
    write_le<uint32_t>(os, blockSize);

    vector<uint8_t> all{ istreambuf_iterator<char>(is), istreambuf_iterator<char>() };
    cerr << "[INFO] input bytes = " << all.size() << "\n";
    const uint64_t n = all.size();
    write_le<uint64_t>(os, n);

    size_t off = 0;
    vector<uint8_t> block, bwt;

    while (off < n) {
        const uint32_t blen = (uint32_t)min<uint64_t>(blockSize, n - off);
        block.assign(all.begin() + off, all.begin() + off + blen);

    #ifdef TEST_ONLY_RANGECODER_IO
        // RANGE CODER: encode raw bytes som 0..256-1 (TOTAL=257)
        ostringstream oss(ios::out | ios::binary);
        RangeEncoder enc(oss);
        const uint32_t TOTAL = 257;
        for (uint8_t b : block) {
            uint32_t lo = (uint32_t)b, hi = lo + 1;
            enc.encode(lo, hi, TOTAL);
        }
        enc.finish();
        const string payload = oss.str();

        uint32_t primary = 0; // ikke brukt
        cerr << "[RC-ONLY] blen=" << blen << " cbytes=" << payload.size() << "\n";

        write_le<uint32_t>(os, blen);
        write_le<uint32_t>(os, primary);
        write_le<uint32_t>(os, (uint32_t)payload.size());
        os.write(payload.data(), payload.size());
    #else
        //BWT+LZP+RC PIPELINE
        vector<uint8_t> bwt;
        uint32_t primary = 0;
        bwt_transform(block, bwt, primary);

        ostringstream oss(ios::out | ios::binary);
        RangeEncoder enc(oss);
        RangeModel m;
        LZP lzp;
        lzp.encode(bwt, [&](int sym){ m.encodeSymbol(enc, sym); });
        enc.finish();
        const string payload = oss.str();

        cerr << "[DBG] payload bytes = " << payload.size() << "\n";

        write_le<uint32_t>(os, blen);
        write_le<uint32_t>(os, primary);
        write_le<uint32_t>(os, (uint32_t)payload.size());
        os.write(payload.data(), payload.size());
    #endif

        if (!os) throw runtime_error("I/O error while writing block");
        off += blen;
    }


    os.flush();        
}

void decompress_file(const string& inPath, const string& outPath) {
    ifstream is(inPath, ios::binary);
    if(!is) throw runtime_error("Kan ikke åpne input");
    char magic[4]; is.read(magic, 4);
    if(is.gcount()!=4 || strncmp(magic, "BSC1", 4)!= 0) throw runtime_error("Bad magic :(");
    uint32_t blockSize = read_le<uint32_t>(is);
    (void)blockSize;
    uint64_t orig = read_le<uint64_t>(is);
    if (orig == 0) { ofstream os(outPath, ios::binary); return; }

    vector<uint8_t> out;
    out.reserve(orig);

    while (out.size() < orig) {
        if (!is) throw runtime_error("EOF blokk");
        uint32_t blen    = read_le<uint32_t>(is);
        uint32_t primary = read_le<uint32_t>(is);
        uint32_t cbytes  = read_le<uint32_t>(is);

        string payload(cbytes, '\0');
        is.read(payload.data(), cbytes);
        if ((uint32_t)is.gcount() != cbytes) throw runtime_error("EOF payload");
        if (cbytes < 4) throw runtime_error("Compressed block too small (need >= 4 bytes)");

        istringstream pis(payload, ios::in | ios::binary);
        RangeDecoder dec(pis, cbytes);

#ifdef TEST_ONLY_RANGECODER_IO
        // RANGE CODER: decode raw bytes (TOTAL=257)
        const uint32_t TOTAL = 257;
        vector<uint8_t> raw; raw.reserve(blen);
        for (uint32_t i = 0; i < blen; ++i) {
            uint32_t offv = dec.get_count(TOTAL);
            dec.remove_symbol(offv, offv + 1, TOTAL);
            raw.push_back((uint8_t)offv);
        }
        if (dec.bytes_read != cbytes) {
            ostringstream msg; msg << "Range decoder consumed " << dec.bytes_read
                                   << " bytes, but payload is " << cbytes;
            throw runtime_error(msg.str());
        }
        out.insert(out.end(), raw.begin(), raw.end());
#else
        // Full pipeline path: LZP + inverse BWT
        RangeModel m; vector<uint8_t> lzp_out; LZP lzp;
        auto nextSym = [&](){ return m.decodeSymbol(dec); };
        lzp.decode(blen, nextSym, lzp_out);

        if (dec.bytes_read != cbytes) {
            ostringstream msg; msg << "Range decoder consumed " << dec.bytes_read
                                   << " bytes, but payload is " << cbytes;
            throw runtime_error(msg.str());
        }

        vector<uint8_t> raw;
        bwt_inverse(lzp_out, primary, raw);
        out.insert(out.end(), raw.begin(), raw.end());
#endif
    }

    ofstream os(outPath, ios::binary);
    if(!os) throw runtime_error("Kan ikke åpne output");
    os.write((char*)out.data(), out.size());
}

static void test_rangecoder_model_free() {
    auto encode_seq = [](const vector<int>& seq) -> string {
        ostringstream oss(ios::out | ios::binary);
        RangeEncoder enc(oss);
        const uint32_t TOTAL = 257;
        for (int s : seq) {
            uint32_t lo = (uint32_t)s;
            uint32_t hi = lo + 1;
            enc.encode(lo, hi, TOTAL);
        }
        enc.finish();
        return oss.str();
    };

    auto decode_seq = [](const string& payload, size_t n) -> vector<int> {
        istringstream is(payload, ios::in | ios::binary);
        RangeDecoder dec(is, payload.size());
        const uint32_t TOTAL = 257;
        vector<int> out;
        out.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            uint32_t off = dec.get_count(TOTAL); // 0..256
            uint32_t lo  = off, hi = off + 1;
            dec.remove_symbol(lo, hi, TOTAL);
            out.push_back((int)off);
        }
        if (dec.bytes_read != payload.size())
            throw runtime_error("rangecoder: exact-consumption failed");
        return out;
    };

    auto round = [&](const vector<int>& seq) {
        string pay = encode_seq(seq);
        auto got = decode_seq(pay, seq.size());
        if (got != seq) {
            for (size_t i=0;i<seq.size();++i) if (seq[i]!=got[i]) {
                cerr << "[RC] first mismatch at i="<<i<<" exp="<<seq[i]<<" got="<<got[i]<<"\n";
                break;
            }
            throw runtime_error("rangecoder model-free mismatch");
        }
    };

    // Deterministic cases
    round({}); // empty
    round({0}); round({256});
    round({1,2,3,4,5});
    round(vector<int>(1000, 0)); // long run
    { vector<int> z(10000); iota(z.begin(), z.end(), 0); for (auto& v: z) v%=257; round(z); }

    // Pseudo-random (fixed seed)
    {
        vector<int> r(20000);
        uint32_t x = 123456789u;
        for (auto& v : r) { x ^= x<<13; x ^= x>>17; x ^= x<<5; v = (int)(x % 257u); }
        round(r);
    }
    cerr << "[RC] model-free tests passed\n";
}

static void test_rangecoder_with_model() {
    auto encode_seq = [](const vector<int>& seq) -> string {
        ostringstream oss(ios::binary);
        RangeEncoder enc(oss);
        RangeModel m;
        for (int s : seq) m.encodeSymbol(enc, s);
        enc.finish();
        return oss.str();
    };
    auto decode_seq = [](const string& pay, size_t n) -> vector<int> {
        istringstream is(pay, ios::binary);
        RangeDecoder dec(is, pay.size());
        RangeModel m;
        vector<int> out; out.reserve(n);
        for (size_t i=0;i<n;++i) out.push_back(m.decodeSymbol(dec));
        if (dec.bytes_read != pay.size())
            throw runtime_error("rangemodel: exact-consumption failed");
        return out;
    };
    auto round = [&](const vector<int>& seq){
        string p = encode_seq(seq);
        auto got = decode_seq(p, seq.size());
        if (got != seq) throw runtime_error("rangemodel mismatch");
    };

    round({});
    round({0,1,2,3,4,5,256,0,0,0,5,5,5});
    { vector<int> many(20000); uint32_t x=1; for (auto& v: many){ x=1664525*x+1013904223u; v=int(x%257u);} round(many); }
    cerr << "[RC] adaptive-model tests passed\n";
}

static void test_bwt_only() {
    auto round = [](const vector<uint8_t>& s){
        vector<uint8_t> t,r; uint32_t p=~0u;
        bwt_transform(s, t, p);
        bwt_inverse(t, p, r);
        if (r != s) throw runtime_error("BWT roundtrip failed");
    };
    round({});
    round({'b','a','n','a','n','a'});
    {
        vector<uint8_t> a(1<<15); uint32_t x=7;
        for (auto& b: a){ x=1103515245*x+12345; b=(uint8_t)(x>>24); }
        round(a);
    }
    cerr << "[BWT] tests passed\n";
}

static void test_lzp_only() {
    auto round = [](const vector<uint8_t>& data){
        LZP lzpe, lzpd;
        vector<int> syms;
        syms.reserve(data.size());
        // Encode
        {
            uint8_t c1=0,c2=0,c3=0,c4=0;
            for (uint8_t cur : data) {
                uint32_t h = LZP::hash4(c1,c2,c3,c4);
                uint8_t pred = lzpe.table[h];
                if (pred == cur) syms.push_back(256);
                else { syms.push_back(cur); lzpe.table[h] = cur; }
                c1=c2; c2=c3; c3=c4; c4=cur;
            }
        }
        // Decode
        vector<uint8_t> out;
        size_t idx = 0;
        auto nextSym = [&]() { return syms[idx++]; };
        lzpd.decode(data.size(), nextSym, out);
        if (out != data) throw runtime_error("LZP roundtrip failed");
    };
    round({});
    round({0,0,0,0,0});
    round({1,2,3,4,5,6,7});
    {
        vector<uint8_t> a(200000); uint32_t x=9;
        for (auto& b: a){ x=1664525*x+1013904223u; b=(uint8_t)(x>>24); }
        round(a);
    }
    cerr << "[LZP] tests passed\n";
}

static void debug_single_symbol() {
    ostringstream oss(ios::binary);
    RangeEncoder enc(oss);
    enc.encode(1, 2, 257);
    enc.finish();
    string pay = oss.str();

    // Print payload bytes as hex
    std::cerr << "[DBG] payload.size=" << pay.size() << " bytes: ";
    for (unsigned char c : pay) std::cerr << std::hex << std::uppercase << std::setw(2) << std::setfill('0') << (int)c << ' ';
    cerr << std::dec << "\n";

    istringstream is(pay, ios::binary);
    RangeDecoder dec(is, pay.size());
    uint32_t q = dec.range / 257u;
    uint32_t t = dec.get_count(257u);
    cerr << "[DBG] range=0x" << std::hex << dec.range
              << " q=0x" << q
              << " code=0x" << dec.code
              << " t=" << std::dec << t << "\n";
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // unit tests
    if (argc == 2 && string(argv[1]) == "rc-tests") {
        try {
            test_rangecoder_model_free();
            test_rangecoder_with_model();
            test_bwt_only();
            test_lzp_only();
        } catch (const exception& e) {
            cerr << "Range coder tests FAILED: " << e.what() << "\n";
            return 1;
        }
        cerr << "Range coder tests OK\n";
        return 0;
    }

    if (argc != 4) {
        cerr << "Usage:\n"
             << "  " << argv[0] << " c <input> <output.bsc>      # BWT+LZP+RC\n"
             << "  " << argv[0] << " d <input.bsc> <output>      # inverse\n"
             << "  " << argv[0] << " rc-c <input> <output.rc>    # raw RC only\n"
             << "  " << argv[0] << " rc-d <input.rc> <output>    # raw RC only\n";
        return 2;
    }

    try {
        string mode = argv[1], inPath = argv[2], outPath = argv[3];
        if      (mode == "c")     compress_file(inPath, outPath);   // full pipeline
        else if (mode == "d")     decompress_file(inPath, outPath); // full inverse
        else if (mode == "rc-c")  rc_pack_file(inPath, outPath);    // RC only
        else if (mode == "rc-d")  rc_unpack_file(inPath, outPath);  // RC only
        else { cerr << "Unknown mode\n"; return 2; }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}