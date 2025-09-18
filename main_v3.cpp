#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <chrono>
#include <algorithm>
#include <random>
#include "louds-trie.hpp"

using namespace std;
using namespace std::chrono;
using namespace louds;

void verify_keys(Trie* trie, const vector<string>& expected_keys, 
                const vector<string>& unexpected_keys) {
    for (const string& key : expected_keys) {
        if (trie->lookup(key) == -1) {
            cerr << "ERROR: Expected key not found: " << key << endl;
            assert(false);
        }
    }
    for (const string& key : unexpected_keys) {
        if (trie->lookup(key) != -1) {
            cerr << "ERROR: Unexpected key found: " << key << endl;
            assert(false);
        }
    }
}

void test_basic_merge() {
    cout << "Testing Basic Merge" << endl;
    
    Trie trie1;
    trie1.add("apple");
    trie1.add("banana");
    trie1.add("cherry");
    trie1.build();
    
    cout << "Trie 1: " << trie1.n_keys() << " keys, " 
         << trie1.n_nodes() << " nodes, "
         << trie1.size() << " bytes" << endl;
    
    Trie trie2;
    trie2.add("apricot");
    trie2.add("blueberry");
    trie2.add("date");
    trie2.build();
    
    cout << "Trie 2: " << trie2.n_keys() << " keys, "
         << trie2.n_nodes() << " nodes, "
         << trie2.size() << " bytes" << endl;
    
    auto start = high_resolution_clock::now();
    Trie* merged1 = Trie::merge_trie(trie1, trie2);
    auto end = high_resolution_clock::now();
    double time1 = duration_cast<microseconds>(end - start).count();
    
    cout << "Approach 1 (Extract-Merge-Rebuild):" << endl;
    cout << "  Time: " << time1 << " μs" << endl;
    cout << "  Merged: " << merged1->n_keys() << " keys, "
         << merged1->n_nodes() << " nodes, "
         << merged1->size() << " bytes" << endl;
    
    start = high_resolution_clock::now();
    Trie* merged2 = Trie::merge_trie_efficient(trie1, trie2);
    end = high_resolution_clock::now();
    double time2 = duration_cast<microseconds>(end - start).count();
    
    cout << "Approach 2 (Efficient Traversal):" << endl;
    cout << "  Time: " << time2 << " μs" << endl;
    cout << "  Merged: " << merged2->n_keys() << " keys, "
         << merged2->n_nodes() << " nodes, "
         << merged2->size() << " bytes" << endl;
    
    vector<string> expected = {"apple", "apricot", "banana", "blueberry", "cherry", "date"};
    vector<string> unexpected = {"grape", "kiwi", "mango"};
    
    verify_keys(merged1, expected, unexpected);
    verify_keys(merged2, expected, unexpected);
    
    assert(merged1->n_keys() == merged2->n_keys());
    assert(merged1->n_nodes() == merged2->n_nodes());
    
    delete merged1;
    delete merged2;
    cout << "Basic merge test passed" << endl;
}

void test_overlapping_keys() {
    cout << "Testing Overlapping Keys" << endl;
    
    Trie trie1;
    trie1.add("apple");
    trie1.add("banana");
    trie1.add("cherry");
    trie1.add("date");
    trie1.build();
    
    Trie trie2;
    trie2.add("banana");     // duplicate
    trie2.add("cherry");     // duplicate  
    trie2.add("elderberry");
    trie2.add("fig");
    trie2.build();
    
    Trie* merged1 = Trie::merge_trie(trie1, trie2);
    Trie* merged2 = Trie::merge_trie_efficient(trie1, trie2);
    
    cout << "Approach 1: " << merged1->n_keys() << " keys (should be 6)" << endl;
    cout << "Approach 2: " << merged2->n_keys() << " keys (should be 6)" << endl;
    
    vector<string> expected = {"apple", "banana", "cherry", "date", "elderberry", "fig"};
    verify_keys(merged1, expected, {});
    verify_keys(merged2, expected, {});
    
    assert(merged1->n_keys() == 6);
    assert(merged2->n_keys() == 6);
    
    delete merged1;
    delete merged2;
    cout << "Overlapping keys test passed" << endl;
}

void test_empty_merge() {
    cout << "Testing Empty Trie Merge" << endl;
    
    Trie trie1;
    trie1.add("alpha");
    trie1.add("beta");
    trie1.add("gamma");
    trie1.build();
    
    Trie trie2;
    trie2.build();
    
    Trie* merged1 = Trie::merge_trie(trie1, trie2);
    Trie* merged2 = Trie::merge_trie_efficient(trie1, trie2);
    
    assert(merged1->n_keys() == 3);
    assert(merged2->n_keys() == 3);
    
    Trie* merged3 = Trie::merge_trie(trie2, trie1);
    Trie* merged4 = Trie::merge_trie_efficient(trie2, trie1);
    
    assert(merged3->n_keys() == 3);
    assert(merged4->n_keys() == 3);
    
    vector<string> expected = {"alpha", "beta", "gamma"};
    verify_keys(merged1, expected, {});
    verify_keys(merged2, expected, {});
    verify_keys(merged3, expected, {});
    verify_keys(merged4, expected, {});
    
    delete merged1;
    delete merged2;
    delete merged3;
    delete merged4;
    cout << "Empty merge test passed" << endl;
}

void test_common_prefixes() {
    cout << "Testing Common Prefixes" << endl;
    
    Trie trie1;
    trie1.add("car");
    trie1.add("card");
    trie1.add("care");
    trie1.add("careful");
    trie1.add("carefully");
    trie1.build();
    
    Trie trie2;
    trie2.add("car");
    trie2.add("cargo");
    trie2.add("career");
    trie2.add("carrot");
    trie2.add("carry");
    trie2.build();
    
    Trie* merged1 = Trie::merge_trie(trie1, trie2);
    Trie* merged2 = Trie::merge_trie_efficient(trie1, trie2);
    
    cout << "Approach 1: " << merged1->n_keys() << " keys, "
         << merged1->n_nodes() << " nodes" << endl;
    cout << "Approach 2: " << merged2->n_keys() << " keys, "
         << merged2->n_nodes() << " nodes" << endl;
    
    vector<string> expected = {
        "car", "card", "care", "career", "careful", 
        "carefully", "cargo", "carrot", "carry"
    };
    
    verify_keys(merged1, expected, {});
    verify_keys(merged2, expected, {});
    
    assert(merged1->n_keys() == 9);
    assert(merged2->n_keys() == 9);
    
    delete merged1;
    delete merged2;
    cout << "Common prefixes test passed" << endl;
}

void test_performance_comparison() {
    cout << "Performance Comparison" << endl;
    
    vector<string> words1, words2;
    
    vector<string> prefixes = {"app", "ban", "car", "dat", "ele", "fig", "gra", "hon"};
    vector<string> suffixes = {"", "le", "ana", "rot", "eful", "efully", "ing", "ed", "er"};
    
    for (const auto& prefix : prefixes) {
        for (const auto& suffix : suffixes) {
            if ((prefix[0] - 'a') % 2 == 0) {
                words1.push_back(prefix + suffix);
            } else {
                words2.push_back(prefix + suffix);
            }
        }
    }
    
    for (int i = 0; i < 10; i++) {
        string word = "common" + to_string(i);
        words1.push_back(word);
        words2.push_back(word);
    }
    
    sort(words1.begin(), words1.end());
    sort(words2.begin(), words2.end());
    
    Trie trie1, trie2;
    for (const auto& word : words1) {
        trie1.add(word);
    }
    trie1.build();
    
    for (const auto& word : words2) {
        trie2.add(word);
    }
    trie2.build();
    
    cout << "Trie 1: " << trie1.n_keys() << " keys" << endl;
    cout << "Trie 2: " << trie2.n_keys() << " keys" << endl;
    
    auto start = high_resolution_clock::now();
    Trie* merged1 = Trie::merge_trie(trie1, trie2);
    auto end = high_resolution_clock::now();
    double time1 = duration_cast<microseconds>(end - start).count();
    
    start = high_resolution_clock::now();
    start = high_resolution_clock::now();
    Trie* merged2 = Trie::merge_trie_efficient(trie1, trie2);
    end = high_resolution_clock::now();
    double time2 = duration_cast<microseconds>(end - start).count();
    
    cout << "\nResults:" << endl;
    cout << "  Approach 1 (Extract-Merge-Rebuild): " << time1 << " μs" << endl;
    cout << "  Approach 2 (Efficient Traversal): " << time2 << " μs" << endl;
    cout << "  Speedup: " << (time1 / time2) << "x" << endl;
    cout << "  Merged trie: " << merged1->n_keys() << " keys, "
         << merged1->n_nodes() << " nodes" << endl;
    
    assert(merged1->n_keys() == merged2->n_keys());
    assert(merged1->n_nodes() == merged2->n_nodes());
    
    delete merged1;
    delete merged2;
    cout << "Performance comparison complete" << endl;
}

void test_large_scale() {
    cout << "Large Scale Test" << endl;
    
    vector<string> keys1, keys2;
    
    for (int i = 0; i < 1000; i++) {
        keys1.push_back("a" + to_string(i * 2));
        keys2.push_back("a" + to_string(i * 2 + 1));
    }
    
    for (int i = 0; i < 100; i++) {
        string common = "common" + to_string(i);
        keys1.push_back(common);
        keys2.push_back(common);
    }
    
    sort(keys1.begin(), keys1.end());
    sort(keys2.begin(), keys2.end());
    
    Trie trie1, trie2;
    for (const auto& key : keys1) {
        trie1.add(key);
    }
    trie1.build();
    
    for (const auto& key : keys2) {
        trie2.add(key);
    }
    trie2.build();
    
    cout << "Trie 1: " << trie1.n_keys() << " keys, "
         << trie1.size() << " bytes" << endl;
    cout << "Trie 2: " << trie2.n_keys() << " keys, "
         << trie2.size() << " bytes" << endl;
    
    auto start = high_resolution_clock::now();
    Trie* merged1 = Trie::merge_trie(trie1, trie2);
    auto end = high_resolution_clock::now();
    double time1 = duration_cast<milliseconds>(end - start).count();
    
    start = high_resolution_clock::now();
    Trie* merged2 = Trie::merge_trie_efficient(trie1, trie2);
    end = high_resolution_clock::now();
    double time2 = duration_cast<milliseconds>(end - start).count();
    
    cout << "\nMerge results:" << endl;
    cout << "  Approach 1: " << time1 << " ms, "
         << merged1->n_keys() << " keys, "
         << merged1->size() << " bytes" << endl;
    cout << "  Approach 2: " << time2 << " ms, "
         << merged2->n_keys() << " keys, "
         << merged2->size() << " bytes" << endl;
    
    for (int i = 0; i < 100; i += 10) {
        string key = "a" + to_string(i);
        assert(merged1->lookup(key) != -1);
        assert(merged2->lookup(key) != -1);
    }
    
    delete merged1;
    delete merged2;
    cout << "Large scale test passed" << endl;
}

int main() {
    cout << "Testing LOUDS Trie Merge" << endl;
    
    try {
        test_basic_merge();
        test_overlapping_keys();
        test_empty_merge();
        test_common_prefixes();
        test_performance_comparison();
        test_large_scale();
        
        cout << "All tests passed successfully" << endl;
    } catch (const exception& e) {
        cerr << "Test failed with exception: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}