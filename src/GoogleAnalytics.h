
#include <string>

#ifndef GOOGLE_ANALYTICS_H_
#define GOOGLE_ANALYTICS_H_

class GoogleAnalytics
{

public:
    static bool parseUrl(std::string url, std::string *host, std::string *path);
    static void httpRequest(const std::string &url, const std::string &data);
    static void postToAnalytics(const std::string &appName);
};


#endif