#ifndef _Arithmetic_H
#define _Arithmetic_H

#include <iostream>

using namespace std;


class Arithmetic
{
public:
  template <typename T>
    static T add(T lhs, T rhs);
  template <typename T>
    static void addEq(T & lhs, T rhs);

  template <typename T>
    static T sub(T lhs, T rhs);
  template <typename T>
    static void subEq(T & lhs, T rhs);

  template <typename T>
    static T mult(T lhs, T rhs);
  template <typename T>
    static void multEq(T & lhs, T rhs);

  template <typename T>
    static T div(T lhs, T rhs);
  template <typename T>
    static void divEq(T & lhs, T rhs);

  template <typename T>
    static T min(T lhs, T rhs);
  template <typename T>
    static void minEq(T & lhs, T rhs);

  template <typename T>
    static T max(T lhs, T rhs);
  template <typename T>
    static void maxEq(T & lhs, T rhs);

private:

};

#endif

