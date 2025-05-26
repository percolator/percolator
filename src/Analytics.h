#ifndef ANALYTICS_H
#define ANALYTICS_H

#include <boost/asio.hpp>
#include <boost/functional/hash.hpp>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>
#include "Version.h"

using boost::asio::ip::tcp;

// Define your PostHog API key here
static const std::string POSTHOG_API_KEY = "phc_VdyCaLRunO8zJrq6SY31nezkC368ORUOdFtOj0FPeyQ";

std::string getPlatform() {
#ifdef _WIN32
    return "Windows";
#elif __APPLE__
    return "macOS";
#elif __linux__
    return "Linux";
#else
    return "Unknown";
#endif
}

// Generate a hashed host ID using hostname and platform
std::string generateHostId() {
    std::size_t seed = 0;
    boost::hash_combine(seed, boost::asio::ip::host_name());
    boost::hash_combine(seed, getPlatform());
    std::ostringstream oss;
    oss << std::hex << seed;
    return oss.str();
}

// Helper to escape string values for JSON
std::string jsonEscape(const std::string& input) {
    std::ostringstream ss;
    for (char c : input) {
        switch (c) {
            case '"': ss << "\\\""; break;
            case '\\': ss << "\\\\"; break;
            case '\b': ss << "\\b"; break;
            case '\f': ss << "\\f"; break;
            case '\n': ss << "\\n"; break;
            case '\r': ss << "\\r"; break;
            case '\t': ss << "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20)
                    ss << "\\u" << std::hex << std::setw(4) << std::setfill('0') << (int)c;
                else
                    ss << c;
        }
    }
    return ss.str();
}

void postToPostHog(const std::string& eventName) {
    try {
        boost::asio::io_context io_context;

        tcp::resolver resolver(io_context);
        auto endpoints = resolver.resolve("app.posthog.com", "80");
        tcp::socket socket(io_context);
        boost::asio::connect(socket, endpoints);

        std::string platform = getPlatform();
        std::string hostId = generateHostId();

        std::ostringstream json;
        json << "{"
             << "\"event\":\"" << jsonEscape(eventName) << "\"," 
             << "\"distinct_id\":\"" << jsonEscape(hostId) << "\"," 
             << "\"properties\":{"
             <<   "\"platform\":\"" << jsonEscape(platform) << "\"," 
             <<   "\"source\":\"percolator\"," 
             <<   "\"library\":\"boost.asio\""
             << "}}";

        std::string body = json.str();

        std::ostringstream request;
        request << "POST /capture/ HTTP/1.1\r\n"
                << "Host: app.posthog.com\r\n"
                << "Content-Type: application/json\r\n"
                << "User-Agent: Percolator/" << VERSION << "\r\n"
                << "Content-Length: " << body.length() << "\r\n"
                << "Authorization: Bearer " << POSTHOG_API_KEY << "\r\n"
                << "Connection: close\r\n\r\n"
                << body;

        boost::asio::write(socket, boost::asio::buffer(request.str()));

        boost::asio::streambuf response;
        boost::asio::read_until(socket, response, "\r\n");

        std::istream response_stream(&response);
        std::string http_version;
        unsigned int status_code;
        response_stream >> http_version >> status_code;


        if (status_code != 200) {
            std::cerr << "PostHog request failed, HTTP status: " << status_code << "\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "PostHog error: " << e.what() << "\n";
    }
}

#endif // ANALYTICS_H
