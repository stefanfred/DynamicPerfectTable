#ifndef DYNAMICPERFECTHASHING_DYNPERFECTHASH_HPP
#define DYNAMICPERFECTHASHING_DYNPERFECTHASH_HPP


#include <vector>
#include <cstdint>
#include <tuple>
#include <algorithm>
#include <limits>
#include <iostream>


inline uint64_t computeM_u32(uint32_t d) {
    return UINT64_C(0xFFFFFFFFFFFFFFFF) / d + 1;
}
inline uint64_t mul128_u32(uint64_t lowbits, uint32_t d) {
    return ((__uint128_t) lowbits * d) >> 64;
}
inline uint32_t fastdiv_u32(uint32_t a, uint64_t M) {
    return (uint32_t) (mul128_u32(M, a));
}
inline uint32_t fastmod_u32(uint32_t a, uint64_t M, uint32_t d) {
    uint64_t lowbits = M * a;
    return (uint32_t) (mul128_u32(lowbits, d));
}

template<typename H, typename V>
class DynPerfectHashTable {
private:
    typedef uint8_t S;
    size_t n;
    size_t bucketCount;
    size_t capacity;
    std::vector<S> seeds;
    std::vector<H> keys;
    std::vector<V> values;
    std::vector<bool> occupied;
    static constexpr uint64_t lambda_log = 2;
    static constexpr double maxLoad = 0.50;
    uint64_t bucketMask;
    uint64_t slotShift;

    uint64_t slotM;
    uint64_t bucketM;

    inline size_t getSlot(H h, S s) {
        auto i = s ^ h;
        i = (i*i)>>32;
        return (i*capacity)>>32;
        //return fastmod_u32(uint32_t((i*i)>>32), slotM, capacity);
        //return (i * i) >> slotShift;
    }

    inline size_t getBucket(H h) {
        return (uint32_t(h) * bucketCount)>>32;
        //return fastmod_u32(uint32_t(h), bucketM, bucketCount);
        //return h & bucketMask;
    }

    // build and swap
    DynPerfectHashTable(const DynPerfectHashTable<H, V> &other) {
    }

    void rebuildAll() {
        DynPerfectHashTable(this);
    }

    void setCapacity(size_t keyCount) {
        n = keyCount;
        size_t power = 1;
        slotShift = 64;
        auto cap = static_cast<size_t>(double(keyCount) / maxLoad);
        /*while (power < cap) {
            power *= 2;
            slotShift--;
        }
        capacity = power;
        bucketCount = capacity >> lambda_log;*/

        capacity = cap;
        slotM = computeM_u32(cap);
        bucketCount = cap / 4;
        bucketM = computeM_u32(bucketCount);

        bucketMask = bucketCount - 1;
        seeds.resize(bucketCount);
        keys.resize(capacity);
        values.resize(capacity);
        occupied.resize(capacity);
    }

    bool searchBucket(const std::vector<std::tuple<H, V>> &pairs, size_t bucketIndex) {
        std::vector<size_t> slots;
        slots.resize(pairs.size());
        for (S seed = 0; seed < std::numeric_limits<S>::max(); seed++) {
            bool succ = true;
            for (size_t i = 0; i < pairs.size(); i++) {
                const std::tuple<H, V> &p = pairs[i];
                size_t slot = getSlot(std::get<0>(p), seed);
                if (occupied[slot]) {
                    succ = false;
                    break;
                }
                slots[i] = slot;
            }
            if (succ) {
                for (size_t i = 0; i < pairs.size(); i++) {
                    const std::tuple<H, V> &p = pairs[i];
                    const size_t slot = slots[i];
                    keys[slot] = std::get<0>(p);
                    values[slot] = std::get<1>(p);
                }
                std::sort(slots.begin(), slots.end());
                if (std::unique(slots.begin(), slots.end()) == slots.end()) {
                    seeds[bucketIndex] = seed;
                    for (const auto slot: slots) {
                        occupied[slot] = true;
                    }
                    return true;
                }
            }
        }
        exit(111);
        return false;
    }

public:
    DynPerfectHashTable(size_t size) {

    }

    DynPerfectHashTable(std::vector<std::tuple<H, V>> pairs) {
        std::vector<std::tuple<size_t, std::vector<std::tuple<H, V>>>> buckets;
        setCapacity(pairs.size());
        buckets.resize(bucketCount);
        for (size_t i = 0; i < bucketCount; i++) {
            std::get<0>(buckets[i]) = i;
        }
        for (const std::tuple<H, V> &pair: pairs) {
            std::get<1>(buckets[getBucket(std::get<0>(pair))]).push_back(pair);
        }
        std::cout << "Keys reorganized" << std::endl;
        std::sort(buckets.begin(), buckets.end(), [](const auto &a, const auto &b) {
            return std::get<1>(a).size() > std::get<1>(b).size();
        });
        std::cout << "Keys sorted" << std::endl;
        size_t cnt = 0;
        for (size_t i = 0; i < buckets.size(); i++) {
            const std::tuple<size_t, std::vector<std::tuple<H, V>>> &b = buckets[i];
            searchBucket(std::get<1>(b), std::get<0>(b));
            cnt += std::get<1>(b).size();
            if (i % 100000 == 0) {
                std::cout << "Found " << (double(cnt) / double(pairs.size())) << " " << std::get<1>(b).size()
                          << std::endl;
            }
        }

        std::cout << "Bits " << (sizeof(S) * 8 >> lambda_log) << std::endl;
        std::cout << "Load " << (double(n)/double(capacity)) << std::endl;
    }

    inline size_t query(H h) {
        size_t bucket = getBucket(h);
        S s = seeds[bucket];
        return getSlot(h, s);
    }

    DynPerfectHashTable() : DynPerfectHashTable(16) {

    }

    inline void update(H h, V v) {

    }

    inline void add(H h, V v) {

    }

    inline bool contains(H h) {

    }

    void iterate() {

    }

    inline V get(H h) {
        return values[query(h)];
    }
};


#endif //DYNAMICPERFECTHASHING_DYNPERFECTHASH_HPP
