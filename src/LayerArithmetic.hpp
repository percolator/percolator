#ifndef _LAYERARITHMETIC_HPP
#define _LAYERARITHMETIC_HPP

#include "PrimitiveVector.hpp"

#include <math.h>

class LayerArithmetic {
private:
  double _alpha;
  PrimitiveVector<unsigned long> _partition_indices;
  unsigned long _n_supported;

  double _power_of_alpha;
  double _partition_index;

  void _add_partition_index_for_next_layer() {
    _partition_indices.push_back( (unsigned long)(ceil(_partition_index)) );

    _power_of_alpha *= _alpha;
    _partition_index += _power_of_alpha;

    _n_supported = _partition_indices.back()+1;
  }

public:
  LayerArithmetic(double alpha):
    _alpha(alpha),
    _n_supported(0),
    _power_of_alpha(1.),
    _partition_index(1.)
  {
    assert(alpha>=1.);
  }

  LayerArithmetic(double alpha, unsigned long n_elements_to_support):
    LayerArithmetic(alpha)
  {
    guarantee_n_elements_supported(n_elements_to_support);
    // Supports n elements, but no more:
    if (n_elements_to_support <= _partition_indices.back()+1)
      _partition_indices.pop_back();
  }

  unsigned long get_partition_index(unsigned long partition_i) const {
    assert(partition_i < _partition_indices.size());
    return _partition_indices[partition_i];
  }

  unsigned long get_partition_index(unsigned long partition_i) {
    while (partition_i >= _partition_indices.size())
      _add_partition_index_for_next_layer();
    return _partition_indices[partition_i];
  }

  unsigned long get_flat_layer_start_index(unsigned long partition_i) {
    if (partition_i == 0)
      return 0;
    return get_partition_index(partition_i-1);
  }

  unsigned long get_flat_layer_end_index(unsigned long partition_i) {
    return get_partition_index(partition_i);
  }

  void guarantee_n_elements_supported(unsigned long n) {
    while (_n_supported < n)
      _add_partition_index_for_next_layer();
  }

  double alpha() const {
    return _alpha;
  }

  unsigned long* begin() {
    return _partition_indices.begin();
  }
  unsigned long* end() {
    return _partition_indices.end();
  }
  unsigned long n_partitions() const {
    return _partition_indices.size();
  }
};

#endif
