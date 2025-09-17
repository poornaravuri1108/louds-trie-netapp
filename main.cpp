#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>

#include "louds-trie.hpp"

using namespace std;
using namespace std::chrono;

int main() {
  ios_base::sync_with_stdio(false);
  vector<string> keys;
  string line;
  while (getline(cin, line)) {
    keys.push_back(line);
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
    assert(ids[i] == i);
  }

  random_shuffle(keys.begin(), keys.end());

  begin = high_resolution_clock::now();
  for (uint64_t i = 0; i < keys.size(); ++i) {
    assert(trie.lookup(keys[i]) != -1);
  }
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  cout << "rnd. lookup = " << (elapsed / keys.size()) << " ns/key" << endl;

  // Test merge_trie
  cout << "\n=== Testing merge_trie ===" << endl;
  louds::Trie trie1;
  trie1.add("apple");
  trie1.add("application");
  trie1.add("banana");
  trie1.build();
  cout << "Trie1: " << trie1.n_keys() << " keys" << endl;

  louds::Trie trie2;
  trie2.add("apple");  // duplicate
  trie2.add("apricot");
  trie2.add("cherry");
  trie2.build();
  cout << "Trie2: " << trie2.n_keys() << " keys" << endl;

  louds::Trie* merged = louds::Trie::merge_trie(trie1, trie2);
  cout << "Merged: " << merged->n_keys() << " keys (expected 5)" << endl;

  // Verify lookups
  assert(merged->lookup("apple") >= 0);
  assert(merged->lookup("application") >= 0);
  assert(merged->lookup("apricot") >= 0);
  assert(merged->lookup("banana") >= 0);
  assert(merged->lookup("cherry") >= 0);
  assert(merged->lookup("notfound") == -1);
  cout << "All merge tests passed!" << endl;

  delete merged;

  return 0;
}