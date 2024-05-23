#ifndef DYNAMICPERFECTHASHING_DYNPERFECTHASH_HPP
#define DYNAMICPERFECTHASHING_DYNPERFECTHASH_HPP


#include <vector>
#include <cstdint>
#include <tuple>
#include <algorithm>
#include <limits>
#include <iostream>
#include <random>


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

//#define POWER_OF_TWO_CAPACITY

template<typename H, typename V>
class DynPerfectHashTable {
private:
    typedef uint8_t S;
    typedef uint32_t P;
    size_t n;
    size_t bucketCount;
    size_t capacity;
    std::vector<S> seeds;
    std::vector<H> keys;
    std::vector<bool> occupiedKeys;
    std::vector<V> values;
    std::vector<bool> occupiedValues;
    std::vector<P> bucketPointer;
    static constexpr double maxLoad = 0.85;

#ifdef POWER_OF_TWO_CAPACITY
    static constexpr uint64_t lambda_log = 2;
    uint64_t bucketMask;
    uint64_t slotShift;
#else
    static constexpr double lambda = 4;
#endif


    std::random_device rd;     // Only used once to initialise (seed) engine
    std::mt19937 rng;    // Random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<S> uni; // Guaranteed unbiased


    void addToKeyStorage(H key, size_t bucketIndex) {
        size_t pos = bucketIndex * capacity / bucketCount; // ToDo
        while (occupiedKeys[pos]) pos = (pos + 1 == capacity) ? 0 : (pos + 1);
        occupiedKeys[pos] = true;
        keys[pos] = key;
    }

    void removeFromKeyStorage(H key, size_t bucketIndex) {
        //ToDo


    }

    void getAllKeysOfBucket(size_t bucketIndex, std::vector<H> &keysOfBucket) {
        size_t pos = bucketIndex * capacity / bucketCount;
        H key;
        while (occupiedKeys[pos] && getBucket((key = keys[pos])) != bucketIndex)
            pos = (pos + 1 == capacity) ? 0 : (pos + 1);
        if (!occupiedKeys[pos]) {
            return;
        }
        do {
            if ((getBucket((key = keys[pos])) == bucketIndex)) {
                keysOfBucket.push_back(key);
            }
            pos = (pos + 1 == capacity) ? 0 : (pos + 1);
        } while (occupiedKeys[pos]);
    }

    void getAllKeyValuesOfBucketAndFree(size_t bucketIndex, std::vector<std::tuple<H, V>> &keyValuesOfBucket) {
        std::vector<H> keysOfBucket;
        getAllKeysOfBucket(bucketIndex, keysOfBucket);
        for (const H &key: keysOfBucket) {
            size_t slot = query(key);
            occupiedValues[slot] = false;
            keyValuesOfBucket.push_back({key, values[slot]});
        }
    }

    size_t getBucketSize(size_t bucketIndex) {
        std::vector<H> keysOfBucket;
        getAllKeysOfBucket(bucketIndex, keysOfBucket);
        return keysOfBucket.size();
    }

    uint64_t getCostOfBucket(size_t size) {
        return 1 << (size << 1);
    }

    inline size_t getSlot(H h, S s) {
        auto i = s ^ h;
#ifdef POWER_OF_TWO_CAPACITY
        return (i * i) >> slotShift;
#else
        return ((__uint128_t) (i*i)*capacity) >> 64;
        //i = (i * i) >> 32;
        //return (i * capacity) >> 32;
#endif
    }

    inline size_t getBucket(H h) {
#ifdef POWER_OF_TWO_CAPACITY
        return h & bucketMask;
#else
        //return ((__uint128_t) h*bucketCount) >> 64;
        return (uint32_t(h) * bucketCount) >> 32;
#endif
    }

    void setCapacity(size_t keyCount) {
        n = keyCount;
        auto cap = static_cast<size_t>(double(keyCount) / maxLoad);

#ifdef POWER_OF_TWO_CAPACITY
        size_t power = 1;
        slotShift = 64;
        while (power < cap) {
            power *= 2;
            slotShift--;
        }
        capacity = power;
        bucketCount = capacity >> lambda_log;
        bucketMask = bucketCount - 1;
#else
        capacity = cap;
        bucketCount = int(double(n + lambda) / lambda);
#endif
        seeds.resize(bucketCount);
        keys.resize(capacity);
        values.resize(capacity);
        bucketPointer.resize(capacity);
        occupiedValues.resize(capacity);
        occupiedKeys.resize(capacity);
    }

    bool searchBucketHeuristically(const std::vector<std::tuple<H, V>> &pairs, size_t bucketIndex, uint64_t recDepth) {
        //std::cout << "HEURISTICALLY " << pairs.size() << " " << bucketIndex << " " << recDepth << std::endl;
        S seed = recDepth == 0 ? 0 : uni(rng);
        uint64_t selfCost = getCostOfBucket(pairs.size());
        S prevS = seeds[bucketIndex];
        S seedOfLowestCost;
        std::vector<size_t> slots;
        slots.resize(pairs.size());
        uint64_t lowestCost = std::numeric_limits<uint64_t>::max();
        for (uint64_t cnt = 0;
             cnt < std::numeric_limits<S>::max(); cnt++, seed = (seed + 1) % std::numeric_limits<S>::max()) {
            if (recDepth > 0 && seed == prevS) {
                continue;
            }
            uint64_t currentSeedCost = 0;
            for (size_t i = 0; i < pairs.size(); i++) {
                const std::tuple<H, V> &p = pairs[i];
                size_t slot = getSlot(std::get<0>(p), seed);
                slots[i] = slot;
                if (occupiedValues[slot]) {
                    size_t otherBucketIndex = bucketPointer[slot];
                    if (otherBucketIndex != bucketIndex) {
                        currentSeedCost += getCostOfBucket(getBucketSize(otherBucketIndex));
                        if (currentSeedCost > lowestCost) {
                            break;
                        }
                    }
                }
            }
            if (currentSeedCost < lowestCost) {
                std::sort(slots.begin(), slots.end());
                if (std::unique(slots.begin(), slots.end()) == slots.end()) {
                    lowestCost = currentSeedCost;
                    seedOfLowestCost = seed;

                    if (currentSeedCost < 2 * selfCost) {
                        break;
                    }
                }
            }
        }
        if (lowestCost == std::numeric_limits<uint64_t>::max()) {
            // all local collisions
            return false;
        }

        // collect key value of other buckets
        std::unordered_map<size_t, std::vector<std::tuple<H, V>>> savedBuckets;
        for (size_t i = 0; i < pairs.size(); i++) {
            const std::tuple<H, V> &p = pairs[i];
            size_t slot = getSlot(std::get<0>(p), seedOfLowestCost);
            if (occupiedValues[slot]) {
                size_t otherBucketIndex = bucketPointer[slot];
                if (otherBucketIndex != bucketIndex && !savedBuckets.contains(otherBucketIndex)) {
                    std::vector<std::tuple<H, V>> keyValuesBucket;
                    getAllKeyValuesOfBucketAndFree(otherBucketIndex, keyValuesBucket);
                    savedBuckets[otherBucketIndex] = keyValuesBucket;
                }
            }
        }

        //insert self
        for (size_t i = 0; i < pairs.size(); i++) {
            const std::tuple<H, V> &p = pairs[i];
            size_t slot = getSlot(std::get<0>(p), seedOfLowestCost);
            occupiedValues[slot] = true;
            values[slot] = std::get<1>(p);
            bucketPointer[slot] = bucketIndex;
        }
        seeds[bucketIndex] = seedOfLowestCost;

        //insert other
        for (auto &it: savedBuckets) {
            if (!searchBucketHeuristically(it.second, it.first, recDepth + 1)) {
                return false;
            }
        }
        return true;
    }

    bool searchBucket(const std::vector<std::tuple<H, V>> &pairs, size_t bucketIndex) {
        for (const auto &p: pairs) {
            addToKeyStorage(std::get<0>(p), bucketIndex);
        }
        std::vector<size_t> slots;
        slots.resize(pairs.size());
        for (S seed = 0; seed < std::numeric_limits<S>::max(); seed++) {
            bool succ = true;
            for (size_t i = 0; i < pairs.size(); i++) {
                const std::tuple<H, V> &p = pairs[i];
                size_t slot = getSlot(std::get<0>(p), seed);
                if (occupiedValues[slot]) {
                    succ = false;
                    break;
                }
                slots[i] = slot;
            }
            if (succ) {
                for (size_t i = 0; i < pairs.size(); i++) {
                    const std::tuple<H, V> &p = pairs[i];
                    const size_t slot = slots[i];
                    values[slot] = std::get<1>(p);
                }
                std::sort(slots.begin(), slots.end());
                if (std::unique(slots.begin(), slots.end()) == slots.end()) {
                    seeds[bucketIndex] = seed;
                    for (const auto slot: slots) {
                        occupiedValues[slot] = true;
                        bucketPointer[slot] = bucketIndex;
                    }
                    return true;
                }
            }
        }
        return searchBucketHeuristically(pairs, bucketIndex, 0);
    }

public:
    DynPerfectHashTable(size_t size) {

    }

    DynPerfectHashTable(std::vector<std::tuple<H, V>> pairs) : rng(rd()) {
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

        std::cout << "Bits " << (sizeof(S) *8.0*double (bucketCount)/double(n) ) << std::endl;
        std::cout << "Load " << (double(n) / double(capacity)) << std::endl;
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
