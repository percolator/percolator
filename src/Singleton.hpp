#ifndef _SINGLETON_H
#define _SINGLETON_H

template <typename T, typename ...ARGS>
class Singleton {
public:
  static T & get_instance(ARGS...args) {
    static T instance(args...);
    return instance;
  }
};

#endif
