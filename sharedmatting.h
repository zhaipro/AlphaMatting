#ifndef SHAREDMSTTING_H
#define SHAREDMSTTING_H

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class Point {
public:
    int x, y;
    Point(int _x=0, int _y=0);
};

class Scalar {
public:
    double val[3];
    Scalar(double b=0, double g=0, double r=0);
    Scalar(uint8_t *bgr);
    void operator += (const Scalar &a);
    void operator /= (double a);
};

struct Tuple
{
    Scalar f;
    Scalar b;
    double   sigmaf;
    double   sigmab;

    int flag;
};

struct Ftuple
{
    Scalar f;
    Scalar b;
    double   alphar;
    double   confidence;
};

/*程序中认定Point中 x为行，y为列，可能错误，但对程序结果没有影响*/
class SharedMatting
{
public:
    SharedMatting();
    ~SharedMatting();

    void loadImage(unsigned char *im, unsigned char *trimap, int64_t w, int64_t h);
    void expandKnown(unsigned char *alpha);
    void sample(Point p, vector<Point>& f, vector<Point>& b);
    void gathering();
    void refineSample(unsigned char *alpha);
    void localSmooth(unsigned char *alpha);
    void solveAlpha(unsigned char *alpha);
    void release();

    double mP(int i, int j, Scalar f, Scalar b);
    double nP(int i, int j, Scalar f, Scalar b);
    double eP(int i1, int j1, int i2, int j2);
    double pfP(Point p, vector<Point>& f, vector<Point>& b);
    double aP(int i, int j, double pf, Scalar f, Scalar b);
    double gP(Point p, Point fp, Point bp, double pf);
    double dP(Point s, Point d);
    double sigma2(Point p);
    double distanceColor2(Scalar cs1, Scalar cs2);
    double comalpha(Scalar c, Scalar f, Scalar b);

private:
    unsigned char *m_trimap;

    vector<struct Tuple> tuples;
    vector<struct Ftuple> ftuples;

    int height;
    int width;
    int kI;
    int kG;
    int * unknownIndex;//Unknown的索引信息；
    double kC;

    int step;
    int channels;
    uint8_t* data;
};

#endif
