// Written by Oliver Serang 2009
// see license for more information

#ifndef _CACHE_H
#define _CACHE_H

#include <iostream>
using namespace std;

template <typename C, typename R, typename A>
  class MemberFunction
{
 protected:
  R (C::*func)(const A &) const;
 public:
  virtual ~MemberFunction() {}
  MemberFunction(R (C::*f)(const A &) const)
    {
      func = f;
    }
  virtual R operator ()(const A & arg, const C * obj) const
  {
    return (obj->*func)(arg);
  }
};

template <typename C, typename R, typename A>
  class LastCachedMemberFunction
{
  // note: this class assumes that the object doesn't change between
  // the two calls-- it would be good to verify that somehow. For now,
  // I'm not sure how that's possible. It could only be defined to
  // work for subclasses of some defined class that clears the cache when
  // the object is modified. Unfortunately, I don't know how to do that,
  // since it would necessitate the object calling something every
  // time it uses a non const function. This is an interesting puzzle.
  // 
  // alternatively, it could require the object to be "locked" or
  // something like that. Maybe...

 protected:
  mutable A arg;
  mutable const C * myObj;
  mutable R returnValueGivenArg;
  MemberFunction<C,R,A> mf;

  // this is necessary since the default value of arg may not be tolerated by the function
  // on the other hand, if the default value is tolerated, then the
  // return value won't be set if the function is called the first time
  // with default argument
  mutable bool set;
 public:
  string name;
 virtual ~LastCachedMemberFunction() {}
 LastCachedMemberFunction(R (C::*f)(const A &) const, const string & n) :
  mf(f)
    {
      name = n;
      set = false;
      myObj = NULL;
    }

  virtual const R & operator ()(const A & a, const C * obj) const
  {
    #ifndef NOCACHE
    if ( myObj == obj && arg == a && set )
      return returnValueGivenArg;
    #endif

    returnValueGivenArg = mf.operator()(a, obj);

    myObj = obj;
    arg = a;
    set = true;
    return returnValueGivenArg;
  }
};


#endif