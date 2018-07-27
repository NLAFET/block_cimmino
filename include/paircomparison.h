#ifndef PAIRCOMPARISON_HPP
#define PAIRCOMPARISON_H

enum position_in_pair{first, second};
enum comparison_direction{inferior, superior};

template<class PAIR_TYPE, position_in_pair POSITION, comparison_direction DIRECTION> bool pair_comparison(const PAIR_TYPE& element1, const PAIR_TYPE& element2){
  switch(POSITION){
    case first :
      switch(DIRECTION){
        case inferior: return element1.first < element2.first;
        case superior: return element1.first > element2.first;
      }
    case second :
      switch(DIRECTION){
        case inferior: return element1.second < element2.second;
        case superior: return element1.second > element2.second;
      }
  }
}

#endif
