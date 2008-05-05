template <typename T>
T Arithmetic::add(T lhs, T rhs)
{
  return lhs + rhs;
}

template <typename T>
void Arithmetic::addEq(T & lhs, T rhs)
{
  lhs += rhs;
}

template <typename T>
T Arithmetic::sub(T lhs, T rhs)
{
  return lhs - rhs;
}

template <typename T>
void Arithmetic::subEq(T & lhs, T rhs)
{
  lhs -= rhs;
}

template <typename T>
T Arithmetic::mult(T lhs, T rhs)
{
  return lhs * rhs;
}

template <typename T>
void Arithmetic::multEq(T & lhs, T rhs)
{
  lhs *= rhs;
}

template <typename T>
T Arithmetic::div(T lhs, T rhs)
{
  return lhs / rhs;
}

template <typename T>
void Arithmetic::divEq(T & lhs, T rhs)
{
  lhs /= rhs;
}

template <typename T>
T Arithmetic::min(T lhs, T rhs)
{
  return lhs < rhs ? lhs : rhs;
}

template <typename T>
void Arithmetic::minEq(T & lhs, T rhs)
{
  lhs = lhs < rhs ? lhs : rhs;
}

template <typename T>
T Arithmetic::max(T lhs, T rhs)
{
  return lhs > rhs ? lhs : rhs;
}

template <typename T>
void Arithmetic::maxEq(T & lhs, T rhs)
{
  lhs = lhs > rhs ? lhs : rhs;
}
