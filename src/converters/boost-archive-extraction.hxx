// file      : examples/cxx/tree/binary/boost/boost-archive-insertion.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : not copyrighted - public domain

#ifndef BOOST_ARCHIVE_EXTRACTION_HXX
#define BOOST_ARCHIVE_EXTRACTION_HXX

#include <cstddef> // std::size_t
#include <string>

#include <xsd/cxx/tree/buffer.hxx>
#include <xsd/cxx/tree/istream.hxx>

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
      inline istream<Archive>&
      operator>> (istream<Archive>& s, istream_common::as_size<T>& x)
      {
        std::size_t r;
        s.impl () >> r;
        x.x_ = static_cast<T> (r);
        return s;
      }

      // 8-bit
      //
      template <typename Archive, typename T>
      inline istream<Archive>&
      operator>> (istream<Archive>& s, istream_common::as_int8<T>& x)
      {
        boost::int8_t r;
        s.impl () >> r;
        x.x_ = static_cast<T> (r);
        return s;
      }

      template <typename Archive, typename T>
      inline istream<Archive>&
      operator>> (istream<Archive>& s, istream_common::as_uint8<T>& x)
      {
        boost::uint8_t r;
        s.impl () >> r;
        x.x_ = static_cast<T> (r);
        return s;
      }


      // 16-bit
      //
      template <typename Archive, typename T>
      inline istream<Archive>&
      operator>> (istream<Archive>& s, istream_common::as_int16<T>& x)
      {
        boost::int16_t r;
        s.impl () >> r;
        x.x_ = static_cast<T> (r);
        return s;
      }

      template <typename Archive, typename T>
      inline istream<Archive>&
      operator>> (istream<Archive>& s, istream_common::as_uint16<T>& x)
      {
        boost::uint16_t r;
        s.impl () >> r;
        x.x_ = static_cast<T> (r);
        return s;
      }


      // 32-bit
      //
      template <typename Archive, typename T>
      inline istream<Archive>&
      operator>> (istream<Archive>& s, istream_common::as_int32<T>& x)
      {
        boost::int32_t r;
        s.impl () >> r;
        x.x_ = static_cast<T> (r);
        return s;
      }

      template <typename Archive, typename T>
      inline istream<Archive>&
      operator>> (istream<Archive>& s, istream_common::as_uint32<T>& x)
      {
        boost::uint32_t r;
        s.impl () >> r;
        x.x_ = static_cast<T> (r);
        return s;
      }


      // 64-bit
      //
      template <typename Archive, typename T>
      inline istream<Archive>&
      operator>> (istream<Archive>& s, istream_common::as_int64<T>& x)
      {
        boost::int64_t r;
        s.impl () >> r;
        x.x_ = static_cast<T> (r);
        return s;
      }

      template <typename Archive, typename T>
      inline istream<Archive>&
      operator>> (istream<Archive>& s, istream_common::as_uint64<T>& x)
      {
        boost::uint64_t r;
        s.impl () >> r;
        x.x_ = static_cast<T> (r);
        return s;
      }


      // Boolean
      //
      template <typename Archive, typename T>
      inline istream<Archive>&
      operator>> (istream<Archive>& s, istream_common::as_bool<T>& x)
      {
        bool r;
        s.impl () >> r;
        x.x_ = static_cast<T> (r);
        return s;
      }


      // Floating-point
      //
      template <typename Archive, typename T>
      inline istream<Archive>&
      operator>> (istream<Archive>& s, istream_common::as_float32<T>& x)
      {
        float r;
        s.impl () >> r;
        x.x_ = static_cast<T> (r);
        return s;
      }

      template <typename Archive, typename T>
      inline istream<Archive>&
      operator>> (istream<Archive>& s, istream_common::as_float64<T>& x)
      {
        double r;
        s.impl () >> r;
        x.x_ = static_cast<T> (r);
        return s;
      }

      // Extraction of std::basic_string.
      //

      template <typename Archive, typename C>
      inline istream<Archive>&
      operator>> (istream<Archive>& s, std::basic_string<C>& x)
      {
        s.impl () >> x;
        return s;
      }


      // Extraction of a binary buffer.
      //
      template <typename Archive, typename C>
      istream<Archive>&
      operator>> (istream<Archive>& s, buffer<C>& x)
      {
        std::size_t size;
        s.impl () >> size;
        x.size (size);
        s.impl ().load_binary (x.data (), size);
        return s;
      }
    }
  }
}

#endif // BOOST_ARCHIVE_EXTRACTION_HXX
