/***************************************************************************
 *   Copyright (C) 2005 by Ari Loytynoja   *
 *   ari@ebi.ac.uk   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef DBMATRIX_H
#define DBMATRIX_H

#ifndef RFOR
#define RFOR(i,n) for(i=n; i>=0; i--)
#endif
#ifndef FOR
#define FOR(i,n) for(i=0; i<n; i++)
#endif

#include <string>
#include <iostream>
#include <cassert>
#include <cmath>

class DbMatrix
{
private:
    int x;
    int y;
    int z;
    int w;
    bool xar;
    bool yar;
    bool zar;
    bool war;

    std::string name;
    double* data;
    int i,j,k,l;
public:
    DbMatrix(int x, std::string name="");
    DbMatrix(int x, int y, std::string name="");
    DbMatrix(int x, int y, int z, std::string name="");
    DbMatrix(int x, int y, int z, int w, std::string name="");

    ~DbMatrix();

    void allocate();
    void initialise(double v = 0);

    double g(int xa, int ya=0, int za = 0, int wa = 0)
    { /*
        if (!(xa>=0&&ya>=0&&za>=0&&wa>=0&&xa<x&&ya<y&&za<z&&wa<w))std::cout<<name<<" "<<xa<<" "<<ya<<" "<<za<<" "<<wa<<std::endl;/**/
        assert(xa>=0);
        assert(xa<x);
        assert(ya>=0);
        assert(ya<y);
        assert(za>=0);
        assert(za<z);
        assert(wa>=0);
        assert(wa<w);
        return data[xa + ya*x + za*x*y + wa*x*y*z];
    }
    /*    void s(double v, int xa, int ya=0, int za = 0, int wa = 0) {  if(!(xa>=0&&ya>=0&&za>=0&&wa>=0&&xa<x&&ya<y&&za<z&&wa<w))std::cout<<"s "<<name<<" "<<v<<" "<<xa<<" "<<ya<<" "<<za<<" "<<wa<<std::endl; assert(xa>=0); assert(xa<x); assert(ya>=0); assert(ya<y); assert(za>=0); assert(za<z);  assert(wa>=0); assert(wa<w); data[xa + ya*x + za*x*y + wa*x*y*z] = v; }*/
    void s(double v, int xa, int ya=0, int za = 0, int wa = 0);
    void a(double v, int xa, int ya=0, int za = 0, int wa = 0)
    {
        assert(xa>=0);
        assert(xa<x);
        assert(ya>=0);
        assert(ya<y);
        assert(za>=0);
        assert(za<z);
        assert(wa>=0);
        assert(wa<w);
        data[xa + ya*x + za*x*y + wa*x*y*z] += v;
    }
    void d(double v, int xa, int ya=0, int za = 0, int wa = 0)
    {
        assert(xa>=0);
        assert(xa<x);
        assert(ya>=0);
        assert(ya<y);
        assert(za>=0);
        assert(za<z);
        assert(wa>=0);
        assert(wa<w);
        data[xa + ya*x + za*x*y + wa*x*y*z] /= v;
    }
    void m(double v, int xa, int ya=0, int za = 0, int wa = 0)
    {
        assert(xa>=0);
        assert(xa<x);
        assert(ya>=0);
        assert(ya<y);
        assert(za>=0);
        assert(za<z);
        assert(wa>=0);
        assert(wa<w);
        data[xa + ya*x + za*x*y + wa*x*y*z] *= v;
    }
    void alog(double v, int xa, int ya=0, int za = 0, int wa = 0)
    {
        assert(xa>=0);
        assert(xa<x);
        assert(ya>=0);
        assert(ya<y);
        assert(za>=0);
        assert(za<z);
        assert(wa>=0);
        assert(wa<w);
        data[xa + ya*x + za*x*y + wa*x*y*z] = sumLogs(data[xa + ya*x + za*x*y + wa*x*y*z], v);
    }
    void clog(int xa, int ya=0, int za = 0, int wa = 0)
    {
        assert(xa>=0);
        assert(xa<x);
        assert(ya>=0);
        assert(ya<y);
        assert(za>=0);
        assert(za<z);
        assert(wa>=0);
        assert(wa<w);
        data[xa + ya*x + za*x*y + wa*x*y*z] = log(data[xa + ya*x + za*x*y + wa*x*y*z]);
    }
    void printName()
    {
        std::cout<<"Name "<<name<<": x = "<<x<<", y = "<<y<<", z = "<<z<<", w = "<<w<<std::endl;
    }
    void print();
    void print(int m);
    double sumLogs(double a, double b);

    void resize(int i);
    void copyData(double *tmp,int new_x,int new_y,int new_z,int new_w);

    void allowResize(bool xr, bool yr=false, bool zr=false, bool wr=false);
};

#endif
