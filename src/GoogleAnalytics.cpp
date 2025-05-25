#include "GoogleAnalytics.h"

#include <sstream>
#include <iostream>
#include "Globals.h"
#include "Version.h"

#define  NO_BOOST_DATE_TIME_INLINE
#include <boost/asio.hpp>
#include <boost/functional/hash.hpp>
#include <boost/functional/hash_fwd.hpp>

// adapted from https://github.com/crux-toolkit/crux-toolkit/blob/master/src/util/crux-utils.cpp
bool GoogleAnalytics::parseUrl(std::string url, std::string* host, std::string* path) {
  if (!host || !path) {
    return false;
  }
  // find protocol
  size_t protocolSuffix = url.find("://");
  if (protocolSuffix != std::string::npos) {
    url = url.substr(protocolSuffix + 3);
  }
  size_t pathBegin = url.find('/');
  if (pathBegin == std::string::npos) {
    *host = url;
    *path = "/";
  } else {
    *host = url.substr(0, pathBegin);
    *path = url.substr(pathBegin);
  }
  if (host->empty()) {
    *host = *path = "";
    return false;
  }
  return true;
}

void GoogleAnalytics::httpRequest(const std::string& url, const std::string& data) {
  // Parse URL into host and path components
  std::string host, path;
  if (!parseUrl(url, &host, &path)) {
    if (VERB > 2) {
      std::cerr << "Warning: Failed parsing URL " << url << std::endl;
    }
    return;
  }

  using namespace boost::asio;

  // Establish TCP connection to host on port 80
  boost::asio::io_context service;
  boost::asio::ip::tcp::resolver resolver(service);
  auto endpoints = resolver.resolve(host, "80");
  auto endpoint = endpoints.begin();
  boost::asio::ip::tcp::socket sock(service);
  sock.connect(endpoint);
  
  std::size_t seed = 0;
  boost::hash_combine(seed, ip::host_name());
  boost::hash_combine(seed, sock.local_endpoint().address().to_string());
  std::stringstream stream;
  stream << std::hex << seed;
  
  std::string placeholder = "CID_PLACEHOLDER";
  std::string cid = stream.str();
  
  std::string newData(data);
  
  if (VERB > 3) {
    std::cerr << "Analytics data string: " << newData << std::endl;
  }
  
  newData.replace(newData.find(placeholder), placeholder.length(), cid);
  
  // Determine method (GET if no data; otherwise POST)
  std::string method = newData.empty() ? "GET" : "POST";
  std::ostringstream lengthString;
  lengthString << newData.length();
  
  std::string contentLengthHeader = newData.empty()
    ? ""
    : "Content-Length: " + lengthString.str() + "\r\n";
  // Send the HTTP request
  std::string request =
    method + " " + path + " HTTP/1.1\r\n"
    "Host: " + host + "\r\n" +
    contentLengthHeader +
    "Connection: close\r\n"
    "\r\n" + newData;
  sock.send(buffer(request));
}

void GoogleAnalytics::postToAnalytics(const std::string& appName) {
  // Post data to Google Analytics
  // For more information, see: https://developers.google.com/analytics/devguides/collection/protocol/v1/devguide
  try {
    std::stringstream paramBuilder;
    paramBuilder << "v=1"                // Protocol verison
                 << "&tid=UA-71657517-2" // Tracking ID
                 << "&cid=CID_PLACEHOLDER" // Unique device ID
                 << "&t=event"           // Hit type
                 << "&ec=percolator"     // Event category
                 << "&ea=" << appName    // Event action
                 << "&el="               // Event label
#ifdef _MSC_VER
                      "win"
#elif __APPLE__
                      "mac"
#else
                      "linux"
#endif
                   << '-' << VERSION;
    httpRequest(
      "http://www.google-analytics.com/collect",
      paramBuilder.str());
  } catch (...) {
  }
}
