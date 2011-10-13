#ifndef CHECK_VERSION_H
#define CHECK_VERSION_H

#include <stdio.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <netdb.h>
#include <string.h>


class Check_version
{
    int create_tcp_socket();
    char *get_ip(const char *host);

public:
    Check_version(int version);
};

#endif // CHECK_VERSION_H
