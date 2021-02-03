// file      : examples/cxx/tree/binary/boost/boost-archive-insertion.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : not copyrighted - public domain

#ifndef BOOST_ARCHIVE_INSERTION_HXX
#define BOOST_ARCHIVE_INSERTION_HXX

#include <cstddef> // std::size_t
#include <string>

#include <xsd/cxx/tree/buffer.hxx>
#include <xsd/cxx/tree/ostream.hxx>

#include <boost/cstdint.hpp>

namespace xsd
{
  namespace cxx
  {
    namespace tree
    {
      // as_size
      //
      template <typename Archive, typename T>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, ostream_common::as_size<T> x)
      {
        std::size_t v (static_cast<std::size_t> (x.x_));
        s.impl () << v;
        return s;
      }

      // 8-bit
      //
      template <typename Archive, typename T>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, ostream_common::as_int8<T> x)
      {
        boost::int8_t v (static_cast<boost::int8_t> (x.x_));
        s.impl () << v;
        return s;
      }

      template <typename Archive, typename T>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, ostream_common::as_uint8<T> x)
      {
        boost::uint8_t v (static_cast<boost::uint8_t> (x.x_));
        s.impl () << v;
        return s;
      }


      // 16-bit
      //
      template <typename Archive, typename T>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, ostream_common::as_int16<T> x)
      {
        boost::int16_t v (static_cast<boost::int16_t> (x.x_));
        s.impl () << v;
        return s;
      }

      template <typename Archive, typename T>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, ostream_common::as_uint16<T> x)
      {
        boost::uint16_t v (static_cast<boost::uint16_t> (x.x_));
        s.impl () << v;
        return s;
      }


      // 32-bit
      //
      template <typename Archive, typename T>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, ostream_common::as_int32<T> x)
      {
        boost::int32_t v (static_cast<boost::int32_t> (x.x_));
        s.impl () << v;
        return s;
      }

      template <typename Archive, typename T>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, ostream_common::as_uint32<T> x)
      {
        boost::uint32_t v (static_cast<boost::uint32_t> (x.x_));
        s.impl () << v;
        return s;
      }


      // 64-bit
      //
      template <typename Archive, typename T>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, ostream_common::as_int64<T> x)
      {
        boost::int64_t v (static_cast<boost::int64_t> (x.x_));
        s.impl () << v;
        return s;
      }

      template <typename Archive, typename T>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, ostream_common::as_uint64<T> x)
      {
        boost::uint64_t v (static_cast<boost::uint64_t> (x.x_));
        s.impl () << v;
        return s;
      }


      // Boolean
      //
      template <typename Archive, typename T>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, ostream_common::as_bool<T> x)
      {
        bool v (static_cast<bool> (x.x_));
        s.impl () << v;
        return s;
      }


      // Floating-point
      //
      template <typename Archive, typename T>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, ostream_common::as_float32<T> x)
      {
        float v (static_cast<float> (x.x_));
        s.impl () << v;
        return s;
      }

      template <typename Archive, typename T>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, ostream_common::as_float64<T> x)
      {
        double v (static_cast<double> (x.x_));
        s.impl () << v;
        return s;
      }


      // Insertion of std::basic_string.
      //
      template <typename Archive, typename C>
      inline ostream<Archive>&
      operator<< (ostream<Archive>& s, const std::basic_string<C>& x)
      {
        s.impl () << x;
        return s;
      }


      // Insertion of a binary buffer.
      //
      template <typename Archive, typename C>
      ostream<Archive>&
      operator<< (ostream<Archive>& s, const buffer<C>& x)
      {
        // Boost.Serialization needs an lvalue.
        //
        std::size_t size (x.size());
        s.impl () << size;
        s.impl ().save_binary (x.data (), x.size ());
        return s;
      }
    }
  }
}

#endif // BOOST_ARCHIVE_INSERTION_HXX
