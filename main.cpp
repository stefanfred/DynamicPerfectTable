#include <iostream>
#include <random>
#include <chrono>
#include "DynPerfectHash.hpp"

#define DO_NOT_OPTIMIZE(value) asm volatile("" : : "r,m"(value) : "memory")


int main() {


    std::random_device rd;
    std::mt19937_64 gen(0);
    std::uniform_int_distribution<uint64_t> dis;

    size_t size = 1<<22;
    size_t queries = 1e8;
    std::vector<std::tuple<uint64_t , uint64_t>> pairs;
    pairs.reserve(size);
    for (size_t i = 0; i < size; ++i) {
        pairs.emplace_back(dis(gen), dis(gen));
    }

    std::cout<< "create dyn"<<std::endl;
    DynPerfectHashTable<uint64_t,uint64_t> table(pairs);
    std::cout<< "create std"<<std::endl;
    std::unordered_map<uint64_t, uint64_t > map;
    for(auto& p :pairs) {
        map[std::get<0>(p)]=std::get<1>(p);
    }

    std::cout<< "check out"<<std::endl;
    for(auto& p :pairs) {
        if(table.get(std::get<0>(p)) != std::get<1>(p)) {
            exit(1);
        }
    }

    std::cout<< "check perfect"<<std::endl;
    std::vector<bool> taken;
    taken.resize(size*4);
    for(auto& p :pairs) {
        uint64_t pos=table.query(std::get<0>(p));
        if(taken[pos]) {
            exit(1);
        }
        taken[pos]=true;
    }

    std::cout<< "prepare query"<<std::endl;
    std::vector<uint64_t> queryInputs;
    queryInputs.reserve(queries);
    for (int i = 0; i < queries; ++i) {
        uint64_t pos = dis(gen) % size;
        queryInputs.push_back(std::get<0>(pairs[pos]));
    }
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    std::cout<< "measure std"<<std::endl;
    begin = std::chrono::steady_clock::now();
    uint64_t sum=0;
    for (int i = 0; i < queries; ++i) { DO_NOT_OPTIMIZE(map[queryInputs[i]]); }
    end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / double(queries)) << "[ns]" << std::endl;

    std::cout<< "measure dyn"<<std::endl;
    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < queries; ++i) { DO_NOT_OPTIMIZE(table.get(queryInputs[i])); }
    end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / double(queries)) << "[ns]" << std::endl;
    DO_NOT_OPTIMIZE(sum);


    return 0;
}
