#ifndef _QUICKLAYERORDEREDHEAP_HPP
#define _QUICKLAYERORDEREDHEAP_HPP

#include "lohify.hpp"
#include "PrimitiveVector.hpp"
#include <functional>

template <typename T>
class QuickLayerOrderedHeap {
protected:
  PrimitiveVector<unsigned long> _partition_indices;

  T*_data;
  unsigned long _n;

  void _put_min_and_max_as_first_and_last_elements_in_layer() {
    for (unsigned long layer_i=0; layer_i< n_layers(); ++layer_i) {
      // Put min and max elements at start and end of layer:
      T*layer_s = layer_begin(layer_i);
      T*layer_e = layer_end(layer_i);

      T*min_ptr = std::min_element(layer_s, layer_e);
      std::swap(*(layer_s), *min_ptr);

      T*max_ptr = std::max_element(layer_s, layer_e);
      std::swap(*(layer_e-1), *max_ptr);
    }
  }

public:
  QuickLayerOrderedHeap(T* data_param, unsigned long n, std::function<bool(const T&,const T&)> compare = [](const auto & lhs, const auto & rhs){return lhs<rhs;}):
    _data(data_param),
    _n(n)
  {
    max_quick_lohify(data_param, data_param+n, _partition_indices, compare);
    _put_min_and_max_as_first_and_last_elements_in_layer();
  }

  T* layer_begin(unsigned long layer_i) const {
    if (layer_i == 0)
      return _data;
    return _data+_partition_indices[layer_i-1];
  }

  T* layer_end(unsigned long layer_i) const {
    if (layer_i+1 == n_layers())
      return _data+_n;
    return _data+_partition_indices[layer_i];
  }

  unsigned long layer_size(unsigned long layer_i) const {
    return layer_end(layer_i) - layer_begin(layer_i);
  }

  const T & min_in_layer(unsigned long layer_i) const {
    return *layer_begin(layer_i);
  }

  const T & max_in_layer(unsigned long layer_i) const {
    return *(layer_end(layer_i)-1);
  }

  void verify() const {
    // note: only used for debugging -- verifies array is LOHified but is slow
    for (unsigned long layer_i=0; layer_i+1<n_layers(); ++layer_i)
      assert(max_in_layer(layer_i) <= min_in_layer(layer_i+1));
  }

  unsigned long n_layers() const {
    return _partition_indices.size()+1;
  }

  unsigned long n() const {
    return _n;
  }

  friend std::string str(const QuickLayerOrderedHeap & loh) {
    std::string res = "";
    for (unsigned long layer_i=0; layer_i<loh.n_layers(); ++layer_i) {
      res += str(loh.layer_begin(layer_i), loh.layer_end(layer_i));
      if (layer_i+1 != loh.n_layers())
	res += "<= ";
    }
    return res;
  }

};

#endif
