#include "louds-trie.hpp"

#ifdef _MSC_VER
 #include <intrin.h>
 #include <immintrin.h>
#else  // _MSC_VER
 #include <x86intrin.h>
#endif  // _MSC_VER

#include <cassert>
#include <vector>
#include <algorithm>


namespace louds {
namespace {

uint64_t Popcnt(uint64_t x) {
#ifdef _MSC_VER
  return __popcnt64(x);
#else  // _MSC_VER
  return __builtin_popcountll(x);
#endif  // _MSC_VER
}

uint64_t Ctz(uint64_t x) {
#ifdef _MSC_VER
  return _tzcnt_u64(x);
#else  // _MSC_VER
  return __builtin_ctzll(x);
#endif  // _MSC_VER
}

struct BitVector {
  struct Rank {
    uint32_t abs_hi;
    uint8_t abs_lo;
    uint8_t rels[3];

    uint64_t abs() const {
      return ((uint64_t)abs_hi << 8) | abs_lo;
    }
    void set_abs(uint64_t abs) {
      abs_hi = (uint32_t)(abs >> 8);
      abs_lo = (uint8_t)abs;
    }
  };

  vector<uint64_t> words;
  vector<Rank> ranks;
  vector<uint32_t> selects;
  uint64_t n_bits;

  BitVector() : words(), ranks(), selects(), n_bits(0) {}

  uint64_t get(uint64_t i) const {
    return (words[i / 64] >> (i % 64)) & 1UL;
  }
  void set(uint64_t i, uint64_t bit) {
    if (bit) {
      words[i / 64] |= (1UL << (i % 64));
    } else {
      words[i / 64] &= ~(1UL << (i % 64));
    }
  }

  void add(uint64_t bit) {
    if (n_bits % 256 == 0) {
      words.resize((n_bits + 256) / 64);
    }
    set(n_bits, bit);
    ++n_bits;
  }
  // build builds indexes for rank and select.
  void build() {
    uint64_t n_blocks = words.size() / 4;
    uint64_t n_ones = 0;
    ranks.resize(n_blocks + 1);
    for (uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
      ranks[block_id].set_abs(n_ones);
      for (uint64_t j = 0; j < 4; ++j) {
        if (j != 0) {
          uint64_t rel = n_ones - ranks[block_id].abs();
          ranks[block_id].rels[j - 1] = rel;
        }

        uint64_t word_id = (block_id * 4) + j;
        uint64_t word = words[word_id];
        uint64_t n_pops = Popcnt(word);
        uint64_t new_n_ones = n_ones + n_pops;
        if (((n_ones + 255) / 256) != ((new_n_ones + 255) / 256)) {
          uint64_t count = n_ones;
          while (word != 0) {
            uint64_t pos = Ctz(word);
            if (count % 256 == 0) {
              selects.push_back(((word_id * 64) + pos) / 256);
              break;
            }
            word ^= 1UL << pos;
            ++count;
          }
        }
        n_ones = new_n_ones;
      }
    }
    ranks.back().set_abs(n_ones);
    selects.push_back(words.size() * 64 / 256);
  }

  // rank returns the number of 1-bits in the range [0, i).
  uint64_t rank(uint64_t i) const {
    uint64_t word_id = i / 64;
    uint64_t bit_id = i % 64;
    uint64_t rank_id = word_id / 4;
    uint64_t rel_id = word_id % 4;
    uint64_t n = ranks[rank_id].abs();
    if (rel_id != 0) {
      n += ranks[rank_id].rels[rel_id - 1];
    }
    n += Popcnt(words[word_id] & ((1UL << bit_id) - 1));
    return n;
  }
  // select returns the position of the (i+1)-th 1-bit.
  uint64_t select(uint64_t i) const {
    const uint64_t block_id = i / 256;
    uint64_t begin = selects[block_id];
    uint64_t end = selects[block_id + 1] + 1UL;
    if (begin + 10 >= end) {
      while (i >= ranks[begin + 1].abs()) {
        ++begin;
      }
    } else {
      while (begin + 1 < end) {
        const uint64_t middle = (begin + end) / 2;
        if (i < ranks[middle].abs()) {
          end = middle;
        } else {
          begin = middle;
        }
      }
    }
    const uint64_t rank_id = begin;
    i -= ranks[rank_id].abs();

    uint64_t word_id = rank_id * 4;
    if (i < ranks[rank_id].rels[1]) {
      if (i >= ranks[rank_id].rels[0]) {
        word_id += 1;
        i -= ranks[rank_id].rels[0];
      }
    } else if (i < ranks[rank_id].rels[2]) {
      word_id += 2;
      i -= ranks[rank_id].rels[1];
    } else {
      word_id += 3;
      i -= ranks[rank_id].rels[2];
    }
    return (word_id * 64) + Ctz(_pdep_u64(1UL << i, words[word_id]));
  }

  uint64_t size() const {
    return sizeof(uint64_t) * words.size()
      + sizeof(Rank) * ranks.size()
      + sizeof(uint32_t) * selects.size();
  }
};

struct Level {
  BitVector louds;
  BitVector outs;
  vector<uint8_t> labels;
  uint64_t offset;

  Level() : louds(), outs(), labels(), offset(0) {}

  uint64_t size() const;
};

uint64_t Level::size() const {
  return louds.size() + outs.size() + labels.size();
}

inline void child_range(const std::vector<Level>& Lv,
  uint64_t lev, uint64_t node_id,
  uint64_t& b, uint64_t& e) {
  b = e = 0;
  if (lev + 1 >= Lv.size()) return;
  const Level& ch = Lv[lev + 1];
  if (ch.louds.n_bits == 0) return;

  uint64_t start_pos = (node_id != 0) ? (ch.louds.select(node_id - 1) + 1) : 0;
  uint64_t pos = start_pos;
  while (pos < ch.louds.n_bits && !ch.louds.get(pos)) ++pos;  
  uint64_t k = pos - start_pos;                                
  b = start_pos - node_id;                                     
  e = b + k;                                                   
}

inline bool is_terminal_at(const std::vector<Level>& Lv, uint64_t lev_plus_1, uint64_t child_id) {
  if (lev_plus_1 >= Lv.size()) return false;
  const Level& L = Lv[lev_plus_1];
  return (child_id < L.outs.n_bits) && (L.outs.get(child_id) != 0);
}

inline void append_parent_one(std::vector<Level>& out_levels, uint64_t lev, bool is_first_parent_at_level) {
  const uint64_t target = lev + 1;
  if (out_levels.size() <= target) out_levels.resize(target + 1);
  if (lev == 0 && is_first_parent_at_level) return;
  out_levels[target].louds.add(1);
}

inline void emit_child(std::vector<Level>& out_levels, uint64_t lev, uint8_t label) {
  Level& nextL = out_levels[lev + 1];
  if (nextL.louds.n_bits == 0) {
    nextL.louds.add(0);
    nextL.louds.add(1);
  } else {
    nextL.louds.set(nextL.louds.n_bits - 1, 0);
    nextL.louds.add(1);
  }
  nextL.labels.push_back(label);
  nextL.outs.add(0);
}

}  // namespace

class TrieImpl {
 public:
  TrieImpl();
  ~TrieImpl() {}

  void add(const string &key);
  void build();
  int64_t lookup(const string &query) const;

  uint64_t n_keys() const {
    return n_keys_;
  }
  uint64_t n_nodes() const {
    return n_nodes_;
  }
  uint64_t size() const {
    return size_;
  }

  void collect_all_keys(vector<string>& keys) const;
  const vector<Level>& get_levels() const { return levels_; }

 private:
  vector<Level> levels_;
  uint64_t n_keys_;
  uint64_t n_nodes_;
  uint64_t size_;
  string last_key_;

  void collect_keys_recursive(vector<string>& keys, uint64_t node_id, 
    uint64_t level, string& prefix) const;

  // Friend declaration for merge functions
  friend Trie* Trie::merge_trie(const Trie& trie1, const Trie& trie2);
  friend Trie* Trie::merge_trie_direct_quadratic(const Trie& trie1, const Trie& trie2);
  friend Trie* Trie::merge_trie_direct_linear(const Trie& trie1, const Trie& trie2);
};

TrieImpl::TrieImpl()
  : levels_(2), n_keys_(0), n_nodes_(1), size_(0), last_key_() {
  levels_[0].louds.add(0);
  levels_[0].louds.add(1);
  levels_[1].louds.add(1);
  levels_[0].outs.add(0);
  levels_[0].labels.push_back(' ');
}

void TrieImpl::add(const string &key) {
  assert(key > last_key_);
  if (key.empty()) {
    levels_[0].outs.set(0, 1);
    ++levels_[1].offset;
    ++n_keys_;
    last_key_ = key;
    return;
  }
  if (key.length() + 1 >= levels_.size()) {
    levels_.resize(key.length() + 2);
  }
  uint64_t i = 0;
  for ( ; i < key.length(); ++i) {
    Level &level = levels_[i + 1];
    uint8_t byte = key[i];
    if ((i == last_key_.length()) || (byte != level.labels.back())) {
      level.louds.set(levels_[i + 1].louds.n_bits - 1, 0);
      level.louds.add(1);
      level.outs.add(0);
      level.labels.push_back(key[i]);
      ++n_nodes_;
      break;
    }
  }
  for (++i; i < key.length(); ++i) {
    Level &level = levels_[i + 1];
    level.louds.add(0);
    level.louds.add(1);
    level.outs.add(0);
    level.labels.push_back(key[i]);
    ++n_nodes_;
  }
  levels_[key.length() + 1].louds.add(1);
  ++levels_[key.length() + 1].offset;
  levels_[key.length()].outs.set(levels_[key.length()].outs.n_bits - 1, 1);
  ++n_keys_;
  last_key_ = key;
}

void TrieImpl::build() {
  uint64_t offset = 0;
  for (uint64_t i = 0; i < levels_.size(); ++i) {
    Level &level = levels_[i];
    level.louds.build();
    level.outs.build();
    offset += levels_[i].offset;
    level.offset = offset;
    size_ += level.size();
  }
}

int64_t TrieImpl::lookup(const string &query) const {
  if (query.length() >= levels_.size()) {
    return -1;
  }
  uint64_t node_id = 0;
  for (uint64_t i = 0; i < query.length(); ++i) {
    const Level &level = levels_[i + 1];
    uint64_t node_pos;
    if (node_id != 0) {
      node_pos = level.louds.select(node_id - 1) + 1;
      node_id = node_pos - node_id;
    } else {
      node_pos = 0;
    }

    // Linear search implementation
    // for (uint8_t byte = query[i]; ; ++node_pos, ++node_id) {
    //   if (level.louds.get(node_pos) || level.labels[node_id] > byte) {
    //     return -1;
    //   }
    //   if (level.labels[node_id] == byte) {
    //     break;
    //   }
    // }

    // Binary search implementation
    uint64_t end = node_pos;
    uint64_t word = level.louds.words[end / 64] >> (end % 64);
    if (word == 0) {
      end += 64 - (end % 64);
      word = level.louds.words[end / 64];
      while (word == 0) {
        end += 64;
        word = level.louds.words[end / 64];
      }
    }
    end += Ctz(word);
    uint64_t begin = node_id;
    end = begin + end - node_pos;

    uint8_t byte = query[i];
    while (begin < end) {
      node_id = (begin + end) / 2;
      if (byte < level.labels[node_id]) {
        end = node_id;
      } else if (byte > level.labels[node_id]) {
        begin = node_id + 1;
      } else {
        break;
      }
    }
    if (begin >= end) {
      return -1;
    }
  }
  const Level &level = levels_[query.length()];
  if (!level.outs.get(node_id)) {
    return -1;
  }
  return level.offset + level.outs.rank(node_id);
}

void TrieImpl::collect_keys_recursive(vector<string>& keys, uint64_t node_id, 
  uint64_t level_idx, string& prefix) const {
  
  if (level_idx == 0 && levels_[0].outs.get(0)) {
    keys.push_back("");
  }
  
  if (level_idx >= levels_.size()) {
    return;
  }
  
  const Level& level = levels_[level_idx];
  
  if (level_idx > 0 && node_id < level.outs.n_bits && level.outs.get(node_id)) {
    keys.push_back(prefix);
  }
  
  if (level_idx + 1 >= levels_.size()) {
    return;
  }
  
  const Level& child_level = levels_[level_idx + 1];
  
  if (child_level.louds.n_bits == 0) {
    return;
  }
  
  uint64_t child_pos;
  
  if (node_id != 0) {
    if (node_id - 1 >= child_level.louds.rank(child_level.louds.n_bits)) {
      return;
    }
    child_pos = child_level.louds.select(node_id - 1) + 1;
  } else {
    child_pos = 0;
  }
  
  uint64_t child_id = child_pos - node_id;
  
  while (child_pos < child_level.louds.n_bits && !child_level.louds.get(child_pos)) {
    if (child_id < child_level.labels.size()) {
      prefix.push_back(child_level.labels[child_id]);
      collect_keys_recursive(keys, child_id, level_idx + 1, prefix);
      prefix.pop_back();
    }
    child_pos++;
    child_id++;
  }
}

void TrieImpl::collect_all_keys(vector<string>& keys) const {
  string prefix;
  collect_keys_recursive(keys, 0, 0, prefix);
}

Trie::Trie() : impl_(new TrieImpl) {}

Trie::~Trie() {
  delete impl_;
}

void Trie::add(const string &key) {
  impl_->add(key);
}

void Trie::build() {
  impl_->build();
}

int64_t Trie::lookup(const string &query) const {
  return impl_->lookup(query);
}

uint64_t Trie::n_keys() const {
  return impl_->n_keys();
}

uint64_t Trie::n_nodes() const {
  return impl_->n_nodes();
}

uint64_t Trie::size() const {
  return impl_->size();
}

std::vector<std::string> Trie::get_all_keys() const {
  std::vector<std::string> keys;
  impl_->collect_all_keys(keys);
  std::sort(keys.begin(), keys.end());
  return keys;
}

// APPROACH 1: Simple Extract-Merge-Rebuild
Trie* Trie::merge_trie(const Trie& trie1, const Trie& trie2) {
  Trie* merged = new Trie();
  
  vector<string> keys1, keys2;
  trie1.impl_->collect_all_keys(keys1);
  trie2.impl_->collect_all_keys(keys2);

  assert(std::is_sorted(keys1.begin(), keys1.end()));
  assert(std::is_sorted(keys2.begin(), keys2.end()));
  
  vector<string> merged_keys;
  merged_keys.reserve(keys1.size() + keys2.size());
  
  size_t i = 0, j = 0;
  while (i < keys1.size() && j < keys2.size()) {
    if (keys1[i] < keys2[j]) {
      merged_keys.push_back(keys1[i++]);
    } else if (keys1[i] > keys2[j]) {
      merged_keys.push_back(keys2[j++]);
    } else {
      merged_keys.push_back(keys1[i]);
      i++;
      j++;
    }
  }
  
  while (i < keys1.size()) {
    merged_keys.push_back(keys1[i++]);
  }
  while (j < keys2.size()) {
    merged_keys.push_back(keys2[j++]);
  }
  
  for (const string& key : merged_keys) {
    merged->add(key);
  }
  merged->build();
  
  return merged;
}


Trie* Trie::merge_trie_direct_quadratic(const Trie& t1, const Trie& t2) {
  Trie* out = new Trie();
  auto& out_impl   = *out->impl_;
  auto& out_levels = out_impl.levels_;
  const auto& L1   = t1.impl_->get_levels();
  const auto& L2   = t2.impl_->get_levels();

  out_impl.n_keys_  = 0;
  out_impl.n_nodes_ = 1; // root
  out_impl.size_    = 0;
  out_levels.resize(2);

  if (!L1.empty() && L1[0].outs.get(0)) { out_levels[0].outs.set(0,1); ++out_levels[1].offset; ++out_impl.n_keys_; }
  if (!L2.empty() && L2[0].outs.get(0) && !out_levels[0].outs.get(0)) { out_levels[0].outs.set(0,1); ++out_levels[1].offset; ++out_impl.n_keys_; }

  struct Pair { bool h1, h2; uint64_t id1, id2; };
  std::vector<Pair> curr(1, Pair{ !L1.empty(), !L2.empty(), 0, 0 }), next;

  for (uint64_t lev = 0; !curr.empty(); ++lev) {
    if (out_levels.size() <= lev + 1) out_levels.resize(lev + 2);
    next.clear();

    for (size_t pidx = 0; pidx < curr.size(); ++pidx) {
      const auto& p = curr[pidx];

      struct Child { uint8_t lab; bool h1, h2; uint64_t c1, c2; };
      std::vector<Child> children;

      if (p.h1) {
        uint64_t b=0,e=0; child_range(L1, lev, p.id1, b, e);
        for (uint64_t i = b; i < e; ++i) {
          uint8_t lab = L1[lev + 1].labels[i];
          auto it = std::find_if(children.begin(), children.end(),
                                 [lab](const Child& c){ return c.lab == lab; });
          if (it == children.end()) children.push_back({lab, true, false, i, 0});
          else                      { it->h1 = true; it->c1 = i; }
        }
      }
      if (p.h2) {
        uint64_t b=0,e=0; child_range(L2, lev, p.id2, b, e);
        for (uint64_t i = b; i < e; ++i) {
          uint8_t lab = L2[lev + 1].labels[i];
          auto it = std::find_if(children.begin(), children.end(),
                                 [lab](const Child& c){ return c.lab == lab; });
          if (it == children.end()) children.push_back({lab, false, true, 0, i});
          else                      { it->h2 = true; it->c2 = i; }
        }
      }

      append_parent_one(out_levels, lev, (pidx == 0));

      if (children.empty()) {
        continue;
      }

      std::sort(children.begin(), children.end(),
                [](const Child& a, const Child& b){ return a.lab < b.lab; });

      for (const auto& ch : children) {
        emit_child(out_levels, lev, ch.lab);
        bool term = (ch.h1 && is_terminal_at(L1, lev + 1, ch.c1))
                 || (ch.h2 && is_terminal_at(L2, lev + 1, ch.c2));
        if (term) {
          Level& here = out_levels[lev + 1];
          here.outs.set(here.outs.n_bits - 1, 1);
          if (out_levels.size() <= lev + 2) out_levels.resize(lev + 3);
          ++out_levels[lev + 2].offset;      
          ++out_impl.n_keys_;                 
        }
        next.push_back(Pair{ ch.h1, ch.h2, ch.c1, ch.c2 });
        ++out_impl.n_nodes_;                  
      }
    }
    curr.swap(next); 
  }

  out_impl.build();
  return out;
}

Trie* Trie::merge_trie_direct_linear(const Trie& t1, const Trie& t2) {
  Trie* out = new Trie();
  auto& out_impl   = *out->impl_;
  auto& out_levels = out_impl.levels_;
  const auto& L1   = t1.impl_->get_levels();
  const auto& L2   = t2.impl_->get_levels();

  out_impl.n_keys_  = 0;
  out_impl.n_nodes_ = 1; 
  out_impl.size_    = 0;
  out_levels.resize(2);

  if (!L1.empty() && L1[0].outs.get(0)) { out_levels[0].outs.set(0,1); ++out_levels[1].offset; ++out_impl.n_keys_; }
  if (!L2.empty() && L2[0].outs.get(0) && !out_levels[0].outs.get(0)) { out_levels[0].outs.set(0,1); ++out_levels[1].offset; ++out_impl.n_keys_; }

  struct Pair { bool h1, h2; uint64_t id1, id2; };
  std::vector<Pair> curr(1, Pair{ !L1.empty(), !L2.empty(), 0, 0 }), next;

  for (uint64_t lev = 0; !curr.empty(); ++lev) {
    if (out_levels.size() <= lev + 1) out_levels.resize(lev + 2);
    next.clear();

    for (size_t pidx = 0; pidx < curr.size(); ++pidx) {
      const auto& p = curr[pidx];

      uint64_t b1=0,e1=0, b2=0,e2=0;
      if (p.h1) child_range(L1, lev, p.id1, b1, e1);
      if (p.h2) child_range(L2, lev, p.id2, b2, e2);

      append_parent_one(out_levels, lev, (pidx == 0));
      append_parent_one(out_levels, lev, (pidx == 0));

      if (b1 == e1 && b2 == e2) {
        continue; 
      }

      while (b1 < e1 || b2 < e2) {
        bool take1=false, take2=false;
        uint8_t lab=0;

        if (b1 < e1 && b2 < e2) {
          uint8_t l1 = L1[lev + 1].labels[b1];
          uint8_t l2 = L2[lev + 1].labels[b2];
          if (l1 == l2) { lab = l1; take1 = take2 = true; }
          else if (l1 < l2) { lab = l1; take1 = true; }
          else               { lab = l2; take2 = true; }
        } else if (b1 < e1) { lab = L1[lev + 1].labels[b1]; take1 = true; }
        else                { lab = L2[lev + 1].labels[b2]; take2 = true; }

        emit_child(out_levels, lev, lab);

        bool term = (take1 && is_terminal_at(L1, lev + 1, b1))
                 || (take2 && is_terminal_at(L2, lev + 1, b2));
        if (term) {
          Level& here = out_levels[lev + 1];
          here.outs.set(here.outs.n_bits - 1, 1);
          if (out_levels.size() <= lev + 2) out_levels.resize(lev + 3);
          ++out_levels[lev + 2].offset;
          ++out_impl.n_keys_;
        }

        next.push_back(Pair{ take1, take2, take1 ? b1 : 0, take2 ? b2 : 0 });
        ++out_impl.n_nodes_;

        if (take1) ++b1;
        if (take2) ++b2;
      }
    }
    curr.swap(next);  
  }

  out_impl.build();
  return out;
}



}  // namespace louds