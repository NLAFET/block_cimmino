#ifndef PAIRCOMPARISON_HPP
#define PAIRCOMPARISON_H

enum position_in_pair{first, second};
enum comparison_direction{ascending, descending};

template<class T, class K, position_in_pair POSITION, comparison_direction DIRECTION> bool pair_comparison(const std::pair<T, K>& pair1, const std::pair<T, K>& pair2){
  switch(POSITION){
    case first :
      switch(DIRECTION){
        case ascending: return pair1.first < pair2.first;
        case descending: return pair1.first > pair2.first;
      }
    case second :
      switch(DIRECTION){
        case ascending: return pair1.second < pair2.second;
        case descending: return pair1.second > pair2.second;
      }
  }
}

template<class T, class K, comparison_direction DIRECTION> bool elt_comparison(const T& elt1, const K& elt2){
    switch(DIRECTION){
        case ascending: return elt1 < elt2;
        case descending: return elt1 > elt2;
    }
}

template <class ForwardIterator>
  std::size_t max_element_index (ForwardIterator first, ForwardIterator last )
{
  ForwardIterator highest = first;
  std::size_t index = 0;
  std::size_t i = 0;
  if (first==last) return index;
  while (++first!=last) {
    ++i;
    if (*first>*highest) {
      highest=first;
      index = i;
    }
  }
  return index;
}

template <class ForwardIterator>
  std::size_t min_element_index ( ForwardIterator first, ForwardIterator last )
{
  ForwardIterator lowest = first;
  std::size_t index = 0;
  std::size_t i = 0;
  if (first==last) return index;
  while (++first!=last) {
    ++i;
    if (*first<*lowest) {
      lowest=first;
      index = i;
    }
  }
  return index;
}

#endif
