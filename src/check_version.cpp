#include "check_version.h"

#include <stdio.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <netdb.h>
#include <string>
#include <sstream>
#include <iostream>

using namespace std;

#define PORT 80

Check_version::Check_version(int version)
{

    cout<<"\nThis is PRANK v."<<version<<".\nChecking if updates are available at http://http://code.google.com/p/prank-msa.\n";

    struct sockaddr_in *remote;
    char buf[BUFSIZ+1];

    int sock = create_tcp_socket();
    char *ip = get_ip("prank-msa.googlecode.com");

    remote = (struct sockaddr_in *)malloc(sizeof(struct sockaddr_in *));
    remote->sin_family = AF_INET;
    int tmpres = inet_pton(AF_INET, ip, (void *)(&(remote->sin_addr.s_addr)));
    if( tmpres < 0)
    {
        perror("Can't set remote->sin_addr.s_addr");
        exit(1);
    }
    else if(tmpres == 0)
    {
        fprintf(stderr, "%s is not a valid IP address\n", ip);
        exit(1);
    }
    remote->sin_port = htons(PORT);

    if(connect(sock, (struct sockaddr *)remote, sizeof(struct sockaddr)) < 0){
        perror("Could not connect");
        exit(1);
    }

    char get[] = "GET /git/VERSION_HISTORY HTTP/1.0\r\nHost: prank-msa.googlecode.com\r\nUser-Agent: HTMLGET 1.0\r\n\r\n";

    //Send the query to the server
    int sent = 0;
    while(sent < (int)strlen(get))
    {
        tmpres = send(sock, get+sent, strlen(get)-sent, 0);
        if(tmpres == -1)
        {
            perror("Can't send query");
            exit(1);
        }
        sent += tmpres;
    }


    stringstream output;

    //now it is time to receive the page
    memset(buf, 0, sizeof(buf));
    int htmlstart = 0;
    char * htmlcontent;
    while((tmpres = recv(sock, buf, BUFSIZ, 0)) > 0){
        if(htmlstart == 0)
        {
            /* Under certain conditions this will not work.
            * If the \r\n\r\n part is splitted into two messages
            * it will fail to detect the beginning of HTML content
            */
            htmlcontent = strstr(buf, "\r\n\r\n");
            if(htmlcontent != NULL)
            {
                htmlstart = 1;
                htmlcontent += 4;
            }
        }
        else
        {
            htmlcontent = buf;
        }
        if(htmlstart)
        {
            output<<htmlcontent;
        }

        memset(buf, 0, tmpres);
    }
    if(tmpres < 0)
    {
      perror("Error receiving data");
    }

    bool print_this = true;
    bool has_printed = false;
    string s;

    while( getline(output,s) )
    {
        istringstream ss(s);
        int d;
        char v,p;
        while( ss >> v >> p >> d )
        {
            if(v=='v' && p=='.' && int(d*10000) <= int(version*10000)+10)
            {
               print_this = false;
            }
        }

        if(print_this)
        {
            if(!has_printed)
                cout<<"\nFound updates. Changes in the more recent versions:\n\n";

            has_printed = true;
            cout<<s<<endl;
        }
        else
        {
            break;
        }
    }

    if(!has_printed)
        cout<<"\nNo updates are available.\n\n";

    free(remote);
    free(ip);
    close(sock);
}

int Check_version::create_tcp_socket()
{
  int sock;
  if((sock = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0){
    perror("Can't create TCP socket");
    exit(1);
  }
  return sock;
}


char *Check_version::get_ip(const char *host)
{
  struct hostent *hent;
  int iplen = 15; //XXX.XXX.XXX.XXX
  char *ip = (char *)malloc(iplen+1);
  memset(ip, 0, iplen+1);
  if((hent = gethostbyname(host)) == NULL)
  {
    herror("Can't get IP");
    exit(1);
  }
  if(inet_ntop(AF_INET, (void *)hent->h_addr_list[0], ip, iplen) == NULL)
  {
    perror("Can't resolve host");
    exit(1);
  }
  return ip;
}

