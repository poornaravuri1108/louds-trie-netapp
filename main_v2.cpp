#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>
#include <sstream>

#include "louds-trie.hpp"

using namespace std;
using namespace std::chrono;

void test_merge_functionality() {
  cout << "Testing Merge Functionality" << endl;
  
  {
    cout << "Test 1: Basic merge with disjoint sets" << endl;
    louds::Trie trie1, trie2;
    
    trie1.add("apple");
    trie1.add("banana");
    trie1.add("cherry");
    trie1.build();
    cout << "Trie1: " << trie1.n_keys() << " keys, " << trie1.n_nodes() << " nodes" << endl;
    
    trie2.add("date");
    trie2.add("elderberry");
    trie2.add("fig");
    trie2.build();
    cout << "Trie2: " << trie2.n_keys() << " keys, " << trie2.n_nodes() << " nodes" << endl;
    
    louds::Trie* merged = louds::Trie::merge_trie(trie1, trie2);
    cout << "Merged: " << merged->n_keys() << " keys, " << merged->n_nodes() << " nodes" << endl;

    vector<string> merged_keys = merged->get_all_keys();
    cout << "Keys in merged trie: ";
    for (size_t i = 0; i < merged_keys.size(); ++i) {
      if (i > 0) cout << ", ";
      cout << "\"" << merged_keys[i] << "\"";
    }
    cout << endl;
    
    assert(merged->lookup("apple") != -1);
    assert(merged->lookup("banana") != -1);
    assert(merged->lookup("cherry") != -1);
    assert(merged->lookup("date") != -1);
    assert(merged->lookup("elderberry") != -1);
    assert(merged->lookup("fig") != -1);
    assert(merged->lookup("grape") == -1);
    
    cout << "Test 1 passed!" << endl;
    delete merged;
  }
  
  {
    cout << "Test 2: Merge with overlapping keys" << endl;
    louds::Trie trie1, trie2;
    
    trie1.add("apple");
    trie1.add("banana");
    trie1.add("cherry");
    trie1.build();
    
    trie2.add("banana");
    trie2.add("cherry");
    trie2.add("date");
    trie2.build();
    
    louds::Trie* merged = louds::Trie::merge_trie(trie1, trie2);
    cout << "Merged: " << merged->n_keys() << " keys (should be 4)" << endl;
    
    vector<string> merged_keys = merged->get_all_keys();
    cout << "Keys in merged trie: ";
    for (size_t i = 0; i < merged_keys.size(); ++i) {
      if (i > 0) cout << ", ";
      cout << "\"" << merged_keys[i] << "\"";
    }
    cout << endl;

    assert(merged->n_keys() == 4);
    
    assert(merged->lookup("apple") != -1);
    assert(merged->lookup("banana") != -1);
    assert(merged->lookup("cherry") != -1);
    assert(merged->lookup("date") != -1);
    
    cout << "Test 2 passed!" << endl;
    delete merged;
  }
  
  {
    cout << "Test 3: Merge with empty trie" << endl;
    louds::Trie trie1, trie2;
    
    trie1.add("alpha");
    trie1.add("beta");
    trie1.build();
    
    trie2.build();
    trie2.build();
    
    louds::Trie* merged = louds::Trie::merge_trie(trie1, trie2);

    cout << "Merged: " << merged->n_keys() << " keys (should be 2)" << endl;
    vector<string> merged_keys = merged->get_all_keys();
    cout << "Keys in merged trie: ";
    for (size_t i = 0; i < merged_keys.size(); ++i) {
      if (i > 0) cout << ", ";
      cout << "\"" << merged_keys[i] << "\"";
    }
    cout << endl;
    
    assert(merged->n_keys() == 2);
    assert(merged->lookup("alpha") != -1);
    assert(merged->lookup("beta") != -1);
    
    cout << "Test 3 passed!" << endl;
    delete merged;
  }
  
  {
    cout << "Test 4: Common prefixes" << endl;
    louds::Trie trie1, trie2;
    
    trie1.add("car");
    trie1.add("card");
    trie1.add("care");
    trie1.add("careful");
    trie1.build();
    
    trie2.add("car");
    trie2.add("career");
    trie2.add("cargo");
    trie2.build();
    
    louds::Trie* merged = louds::Trie::merge_trie(trie1, trie2);
    cout << "Merged: " << merged->n_keys() << " keys, " << merged->n_nodes() << " nodes" << endl;
    vector<string> merged_keys = merged->get_all_keys();
    cout << "Keys in merged trie: ";
    for (size_t i = 0; i < merged_keys.size(); ++i) {
      if (i > 0) cout << ", ";
      cout << "\"" << merged_keys[i] << "\"";
    }
    cout << endl;
    assert(merged->n_keys() == 6);  // car, card, care, careful, cargo, career
    
    assert(merged->lookup("car") != -1);
    assert(merged->lookup("card") != -1);
    assert(merged->lookup("care") != -1);
    assert(merged->lookup("careful") != -1);
    assert(merged->lookup("cargo") != -1);
    assert(merged->lookup("career") != -1);
    assert(merged->lookup("cart") == -1);  
    
    cout << "Test 4 passed!" << endl;
    delete merged;
  }
  
  {
    cout << "Test 5: Performance test" << endl;
    louds::Trie trie1, trie2;
    
    vector<string> keys1;
    for (int i = 0; i < 1000; i += 2) {
      stringstream ss;
      ss << "key" << i;
      keys1.push_back(ss.str());
    }
    sort(keys1.begin(), keys1.end());
    
    vector<string> keys2;
    for (int i = 1; i < 1000; i += 2) {
      stringstream ss;
      ss << "key" << i;
      keys2.push_back(ss.str());
    }
    for (int i = 100; i < 200; i += 2) {
      stringstream ss;
      ss << "key" << i;
      keys2.push_back(ss.str());
    }
    sort(keys2.begin(), keys2.end());
    
    for (const string& key : keys1) {
      trie1.add(key);
    }
    trie1.build();
    cout << "Trie1: " << trie1.n_keys() << " keys" << endl;
    
    for (const string& key : keys2) {
      trie2.add(key);
    }
    trie2.build();
    cout << "Trie2: " << trie2.n_keys() << " keys" << endl;
    
    auto start = high_resolution_clock::now();
    louds::Trie* merged = louds::Trie::merge_trie(trie1, trie2);
    auto end = high_resolution_clock::now();
    double elapsed = duration_cast<microseconds>(end - start).count();
    
    cout << "Merged: " << merged->n_keys() << " keys in " << elapsed << " μs" << endl;
    cout << "Size: " << merged->size() << " bytes" << endl;
    
    assert(merged->lookup("key100") != -1);
    assert(merged->lookup("key500") != -1);
    assert(merged->lookup("key999") != -1);
    assert(merged->lookup("key1000") == -1);
    
    cout << "Test 5 passed!" << endl;
    delete merged;
  }
  
  cout << "\n=== All Merge Tests Passed! ===" << endl;
}

int main(int argc, char** argv) {
  ios_base::sync_with_stdio(false);
  
  if (argc > 1 && string(argv[1]) == "--test") {
    test_merge_functionality();
    return 0;
  }
  
  cout << "Original LOUDS Trie Test (reading from stdin)" << endl;
  cout << "To test merge functionality, run with --test flag" << endl;
  
  vector<string> keys;
  string line;
  while (getline(cin, line)) {
    keys.push_back(line);
  }
  
  if (keys.empty()) {
    cout << "No input provided. Run with --test to test merge functionality." << endl;
    return 0;
  }

  high_resolution_clock::time_point begin = high_resolution_clock::now();
  louds::Trie trie;
  for (uint64_t i = 0; i < keys.size(); ++i) {
    trie.add(keys[i]);
  }
  trie.build();
  high_resolution_clock::time_point end = high_resolution_clock::now();
  double elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  cout << "build = " << (elapsed / keys.size()) << " ns/key" << endl;

  cout << "#keys = " << trie.n_keys() << endl;
  cout << "#nodes = " << trie.n_nodes() << endl;
  cout << "size = " << trie.size() << " bytes" << endl;

  vector<uint64_t> ids(keys.size());

  begin = high_resolution_clock::now();
  for (uint64_t i = 0; i < keys.size(); ++i) {
    ids[i] = trie.lookup(keys[i]);
  }
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  cout << "seq. lookup = " << (elapsed / keys.size()) << " ns/key" << endl;

  sort(ids.begin(), ids.end());
  for (uint64_t i = 0; i < ids.size(); ++i) {
    assert(ids[i] == (int64_t)i);
  }

  random_shuffle(keys.begin(), keys.end());

  begin = high_resolution_clock::now();
  for (uint64_t i = 0; i < keys.size(); ++i) {
    assert(trie.lookup(keys[i]) != -1);
  }
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  cout << "rnd. lookup = " << (elapsed / keys.size()) << " ns/key" << endl;

  return 0;
}