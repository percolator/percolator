#ifndef _LAYERORDEREDHEAP_HPP
#define _LAYERORDEREDHEAP_HPP

#include "LayerArithmetic.hpp"
#include "lohify.hpp"

template <typename T>
class LayerOrderedHeap {
protected:
  LayerArithmetic*_la;

  T*_data;
  unsigned long _n;
  unsigned long _n_layers;

  void _lohify(T*x, T*x_end, unsigned long*partition_ranks, unsigned long*partition_ranks_end) {
    max_random_lohify(x, x_end, partition_ranks, partition_ranks_end);
  }

  void _put_min_and_max_as_first_and_last_elements_in_layer() {
    for (unsigned long layer_i=0; layer_i< n_layers(); ++layer_i) {
      // Put min and max elements at start and end of layer:
      T*layer_s = layer_begin(layer_i);
      T*layer_e = layer_end(layer_i);

      T*min_ptr = std::min_element(layer_s, layer_e);
      std::swap(*layer_s, *min_ptr);

      T*max_ptr = std::max_element(layer_s, layer_e);
      std::swap(*(layer_e-1), *max_ptr);
    }
  }

public:
  LayerOrderedHeap(T* data_param, unsigned long n, LayerArithmetic*la):
    _la(la),
    _data(data_param),
    _n(n)
  {
    _la->guarantee_n_elements_supported(n);
    
    // Note: log search (e.g., via upper bound) is not used because
    // _la may be shared and have become huge from some other
    // client. The cost of linear search is already amortized out by
    // the cost of _lohify, which is already in Omega(number of pivots).
    unsigned long*pivot_ptr;
    for (pivot_ptr=_la->begin(); pivot_ptr!=_la->end(); ++pivot_ptr) {
      if (*pivot_ptr >= n)
	break;
    }

    _n_layers = pivot_ptr + 1 - _la->begin();

    // pivot_ptr is one past the final relevant pivot:
    _lohify(_data, _data + _n, _la->begin(), pivot_ptr);
    _put_min_and_max_as_first_and_last_elements_in_layer();

    // fixme: keep until release / runtime plots; however, comment out then (it slows things down)
    //verify();
  }

  T* layer_begin(unsigned long layer_i) const {
    if (layer_i == 0)
      return _data;
    return _data+_la->get_partition_index(layer_i-1);
  }

  T* layer_end(unsigned long layer_i) const {
    if (layer_i+1 == n_layers())
      return _data+_n;
    return _data+_la->get_partition_index(layer_i);
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
    for (unsigned long layer_i=0; layer_i+1<n_layers(); ++layer_i)
      assert(max_in_layer(layer_i) <= min_in_layer(layer_i+1));
  }

  unsigned long n_layers() const {
    return _n_layers;
  }

  unsigned long n() const {
    return _n;
  }

  friend std::string str(const LayerOrderedHeap & loh){
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
