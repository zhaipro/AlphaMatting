#ifndef SHAREDMSTTING_H
#define SHAREDMSTTING_H

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class Point {
public:
    int x, y;
    Point();
    Point(int _x, int _y);
};

class Scalar {
public:
    double val[3];
    Scalar();
    Scalar(double b, double g, double r);
};

struct labelPoint
{
    int x;
    int y;
    int label;
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

    void loadImage(unsigned char *_data, int64_t w, int64_t h);
    void loadTrimap(unsigned char *data);
    void expandKnown();
    void sample(Point p, vector<Point>& f, vector<Point>& b);
    void gathering();
    void refineSample(unsigned char *alpha);
    void localSmooth(unsigned char *alpha);
    void solveAlpha(unsigned char *alpha);
    void Sample(vector<vector<Point> > &F, vector<vector<Point> > &B);
    void release();

    double mP(int i, int j, Scalar f, Scalar b);
    double nP(int i, int j, Scalar f, Scalar b);
    double eP(int i1, int j1, int i2, int j2);
    double pfP(Point p, vector<Point>& f, vector<Point>& b);
    double aP(int i, int j, double pf, Scalar f, Scalar b);
    double gP(Point p, Point fp, Point bp, double pf);
    double gP(Point p, Point fp, Point bp, double dpf, double pf);
    double dP(Point s, Point d);
    double sigma2(Point p);
    double distanceColor2(Scalar cs1, Scalar cs2);
    double comalpha(Scalar c, Scalar f, Scalar b);

private:
    unsigned char *pImg;
    unsigned char *trimap;

    vector<Point> uT;
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
