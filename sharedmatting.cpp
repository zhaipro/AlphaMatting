// https://www.inf.ufrgs.br/~eslgastal/SharedMatting/
// https://blog.csdn.net/gzj2013/article/details/82685389
#include "sharedmatting.h"
#include <time.h>

// ?????????
// #define min(a, b) ((a) < (b))?(a):(b)
// #define max(a, b) ((a) > (b))?(a):(b)

double max(double a, double b)
{
    return a < b?b:a;
}

double min(double a, double b)
{
    return a < b?a:b;
}

Point::Point(int _x, int _y):x(_x),y(_y){}

Scalar::Scalar(double b, double g, double r)
{
    val[0] = b;
    val[1] = g;
    val[2] = r;
}

Scalar::Scalar(uint8_t *bgr)
{
    val[0] = bgr[0];
    val[1] = bgr[1];
    val[2] = bgr[2];
}

void Scalar::operator += (const Scalar &a)
{
    val[0] += a.val[0];
    val[1] += a.val[1];
    val[2] += a.val[2];
}

void Scalar::operator /= (double a)
{
    val[0] /= a;
    val[1] /= a;
    val[2] /= a;
}

SharedMatting::SharedMatting()
{
    kI = 10;
    kC = 5.0;
    kG = 4;  //each unknown p gathers at most kG forground and background samples
}

SharedMatting::~SharedMatting()
{
    delete[] unknownIndex;
}

void SharedMatting::loadImage(unsigned char *im, unsigned char *trimap, int64_t w, int64_t h)
{
    height     = h;
    width      = w;
    channels   = 3;
    step       = w * channels;
    data       = im;
    m_trimap = trimap;
    unknownIndex  = new int[height * width];
}

void SharedMatting::expandKnown(unsigned char *alpha)
{
    int kc2 = kC * kC;

    memcpy(alpha, m_trimap, width * height);

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            if (m_trimap[i * width + j] != 0 && m_trimap[i * width + j] != 255)
            {
                int label = -1;
                double dmin = 10000.0;
                Scalar p = Scalar(&data[i * step + j * channels]);

                int i1 = max(0, i - kI);
                int i2 = min(i + kI, height - 1);
                int j1 = max(0, j - kI);
                int j2 = min(j + kI, width - 1);

                for (int k = i1; k <= i2; ++k)
                {
                 for (int l = j1; l <= j2; ++l)
                 {
                     int temp = m_trimap[k * width + l];
                     if (temp != 0 && temp != 255)
                     {
                         continue;
                     }
                     double dis = dP(Point(i, j), Point(k, l));
                     if (dis > dmin)
                     {
                         continue;
                     }

                     Scalar q = Scalar(&data[k * step + l * channels]);
                     double distanceColor = distanceColor2(p, q);

                     if (distanceColor <= kc2)
                     {
                         dmin = dis;
                         label = temp;
                     }
                 }
                }

                if (label != -1)
                {
                    alpha[i * width + j] = label;
                }
            }
        }
    }

    // m_trimap = alpha;
    memcpy(m_trimap, alpha, width * height);
}

double SharedMatting::comalpha(Scalar c, Scalar f, Scalar b)
{
    // 已知：c, f, b
    // 求解：argmin |c - (a * f + (1 - a) * b)|
    // 约束：0 <= a <= 1
    double alpha = ((c.val[0] - b.val[0]) * (f.val[0] - b.val[0]) +
                    (c.val[1] - b.val[1]) * (f.val[1] - b.val[1]) +
                    (c.val[2] - b.val[2]) * (f.val[2] - b.val[2]))
                 / ((f.val[0] - b.val[0]) * (f.val[0] - b.val[0]) +
                    (f.val[1] - b.val[1]) * (f.val[1] - b.val[1]) +
                    (f.val[2] - b.val[2]) * (f.val[2] - b.val[2]) + 0.0000001);
    return min(1.0, max(0.0, alpha));
}

double SharedMatting::mP(int i, int j, Scalar f, Scalar b)
{
    Scalar c = Scalar(&data[i * step + j * channels]);

    double alpha = comalpha(c, f, b);

    double result = sqrt((c.val[0] - alpha * f.val[0] - (1 - alpha) * b.val[0]) * (c.val[0] - alpha * f.val[0] - (1 - alpha) * b.val[0]) +
                         (c.val[1] - alpha * f.val[1] - (1 - alpha) * b.val[1]) * (c.val[1] - alpha * f.val[1] - (1 - alpha) * b.val[1]) +
                         (c.val[2] - alpha * f.val[2] - (1 - alpha) * b.val[2]) * (c.val[2] - alpha * f.val[2] - (1 - alpha) * b.val[2]));
    return result / 255.0;
}

double SharedMatting::nP(int i, int j, Scalar f, Scalar b)
{
    int i1 = max(0, i - 1);
    int i2 = min(i + 1, height - 1);
    int j1 = max(0, j - 1);
    int j2 = min(j + 1, width - 1);

    double  result = 0;

    for (int k = i1; k <= i2; ++k)
    {
        for (int l = j1; l <= j2; ++l)
        {
            double m = mP(k, l, f, b);
            result += m * m;
        }
    }

    return result;
}

double SharedMatting::eP(int i1, int j1, int i2, int j2)
{
    double ci = i2 - i1;
    double cj = j2 - j1;
    double z  = sqrt(ci * ci + cj * cj);

    double ei = ci / (z + 0.0000001);
    double ej = cj / (z + 0.0000001);

    double stepinc = min(1 / (abs(ei) + 1e-10), 1 / (abs(ej) + 1e-10));

    double result = 0;

    Scalar pre = Scalar(&data[i1 * step + j1 * channels]);

    int ti = i1;
    int tj = j1;

    for (double t = 1; ;t += stepinc)
    {
        double inci = ei * t;
        double incj = ej * t;
        int i = int(i1 + inci + 0.5);
        int j = int(j1 + incj + 0.5);

        double z = 1;

        Scalar cur = Scalar(&data[i * step + j * channels]);

        if (ti - i > 0 && tj - j == 0)
        {
            z = ej;
        }
        else if(ti - i == 0 && tj - j > 0)
        {
            z = ei;
        }

        result += ((cur.val[0] - pre.val[0]) * (cur.val[0] - pre.val[0]) +
                   (cur.val[1] - pre.val[1]) * (cur.val[1] - pre.val[1]) +
                   (cur.val[2] - pre.val[2]) * (cur.val[2] - pre.val[2])) * z;
        pre = cur;

        ti = i;
        tj = j;

        if(abs(ci) >= abs(inci) || abs(cj) >= abs(incj))
            break;
    }

    return result;
}

double SharedMatting::pfP(Point p, vector<Point>& f, vector<Point>& b)
{
    double fmin = 1e10;
    vector<Point>::iterator it;
    for (it = f.begin(); it != f.end(); ++it)
    {
        double fp = eP(p.x, p.y, it->x, it->y);
        if (fp < fmin)
        {
            fmin = fp;
        }
    }

    double bmin = 1e10;
    for (it = b.begin(); it != b.end(); ++it)
    {
        double bp = eP(p.x, p.y, it->x, it->y);
        if (bp < bmin)
        {
            bmin = bp;
        }
    }
    return bmin / (fmin + bmin + 1e-10);
}

double SharedMatting::aP(int i, int j, double pf, Scalar f, Scalar b)
{
    Scalar c = Scalar(&data[i * step + j * channels]);

    double alpha = comalpha(c, f, b);

    return pf + (1 - 2 * pf) * alpha;
}

double SharedMatting::dP(Point s, Point d)
{
    return sqrt(double((s.x - d.x) * (s.x - d.x) + (s.y - d.y) * (s.y - d.y)));
}

double SharedMatting::gP(Point p, Point fp, Point bp, double pf)
{
    Scalar f = Scalar(&data[fp.x * step + fp.y * channels]);
    Scalar b = Scalar(&data[bp.x * step + bp.y * channels]);

    double tn = pow(nP(p.x, p.y, f, b), 3);
    double ta = pow(aP(p.x, p.y, pf, f, b), 2);
    double tf = dP(p, fp);
    double tb = pow(dP(p, bp), 4);

    return tn * ta * tf * tb;
}

double SharedMatting::sigma2(Point p)
{
    int xi = p.x;
    int yj = p.y;
    Scalar pc = Scalar(&data[xi * step + yj * channels]);

    int i1 = max(0, xi - 2);
    int i2 = min(xi + 2, height - 1);
    int j1 = max(0, yj - 2);
    int j2 = min(yj + 2, width - 1);

    double result = 0;
    int    num    = 0;

    for (int i = i1; i <= i2; ++i)
    {
        for (int j = j1; j <= j2; ++j)
        {
            Scalar temp = Scalar(&data[i * step + j * channels]);
            result += distanceColor2(pc, temp);
            ++num;
        }
    }

    return result / num;
}

double SharedMatting::distanceColor2(Scalar cs1, Scalar cs2)
{
    return (cs1.val[0] - cs2.val[0]) * (cs1.val[0] - cs2.val[0]) +
           (cs1.val[1] - cs2.val[1]) * (cs1.val[1] - cs2.val[1]) +
           (cs1.val[2] - cs2.val[2]) * (cs1.val[2] - cs2.val[2]);
}

void SharedMatting::sample(Point p, std::vector<Point> &f, std::vector<Point> &b)
{
    int i = p.x;
    int j = p.y;

    double inc   = 360.0 / kG;
    double ca    = inc / 9;
    double angle = (i % 3 * 3 + j % 9) * ca;
    for (int k = 0; k  < kG; ++k)
    {
        bool flagf = false;
        bool flagb = false;

        double z  = (angle + k * inc) / 180 * 3.1415926;
        double ei = sin(z);
        double ej = cos(z);

        double step = min(1.0 / (abs(ei) + 1e-10), 1.0 / (abs(ej) + 1e-10));

        for (double t = 1; ;t += step)
        {
            int ti = int(i + ei * t + 0.5);
            int tj = int(j + ej * t + 0.5);

            if(ti >= height || ti < 0 || tj >= width || tj < 0)
            {
                break;
            }
            int gray = m_trimap[ti * width + tj];

            if (!flagf && gray == 255)
            {
                Point tp = Point(ti, tj);
                f.push_back(tp);
                flagf = true;
            }
            else if (!flagb && gray == 0)
            {
                Point tp = Point(ti, tj);
                b.push_back(tp);
                flagb = true;
            }
            if (flagf && flagb)
            {
                break;
            }
        }
    }
}

void SharedMatting::gathering()
{
    vector<Point> f;
    vector<Point> b;
    vector<Point>::iterator it1;
    vector<Point>::iterator it2;

    int index = 0;

    for(int i = 0; i < height; i++)
    for(int j = 0; j < width; j++)
    {
        if(m_trimap[i * width + j] == 0 || m_trimap[i * width + j] == 255)
            continue;

        f.clear();
        b.clear();
        sample(Point(i, j), f, b);

        double pfp = pfP(Point(i, j), f, b);

        double gmin = 1.0e10;

        Point tf;
        Point tb;

        bool flag = false;
        bool first = true;

        for (it1 = f.begin(); it1 != f.end(); ++it1)
        {
            for (it2 = b.begin(); it2 < b.end(); ++it2)
            {
                double gp = gP(Point(i, j), *(it1), *(it2), pfp);
                if (gp < gmin)
                {
                    gmin = gp;
                    tf   = *it1;
                    tb   = *it2;
                    flag = true;
                }
            }
        }

        struct Tuple st;
        st.flag = -1;
        if (flag)
        {
            st.flag   = 1;
            st.f      = Scalar(&data[tf.x * step +  tf.y * channels]);
            st.b      = Scalar(&data[tb.x * step +  tb.y * channels]);
            st.sigmaf = sigma2(tf);
            st.sigmab = sigma2(tb);
        }

        tuples.push_back(st);
        unknownIndex[i * width + j] = index;
        ++index;
    }
}

void SharedMatting::refineSample(unsigned char *alpha)
{
    ftuples.resize(width * height);
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            Scalar c = Scalar(&data[i * step +  j * channels]);
            int indexf = i * width + j;
            int gray = m_trimap[i * width + j];
            if (gray == 0)
            {
                ftuples[indexf].f = c;
                ftuples[indexf].b = c;
                ftuples[indexf].alphar = 0;
                ftuples[indexf].confidence = 1;
            }
            else if (gray == 255)
            {
                ftuples[indexf].f = c;
                ftuples[indexf].b = c;
                ftuples[indexf].alphar = 1;
                ftuples[indexf].confidence = 1;
            }
        }
    }
    for(int xi = 0; xi < height; xi++)
    for(int yj = 0; yj < width; yj++)
    {
        if(m_trimap[xi * width + yj] == 0 || m_trimap[xi * width + yj] == 255)
            continue;
        int i1 = max(0, xi - 5);
        int i2 = min(xi + 5, height - 1);
        int j1 = max(0, yj - 5);
        int j2 = min(yj + 5, width - 1);

        double minvalue[3] = {1e10, 1e10, 1e10};
        Point p[3];
        int num = 0;
        for (int k = i1; k <= i2; ++k)
        {
            for (int l = j1; l <= j2; ++l)
            {
                int temp = m_trimap[k * width + l];

                if (temp == 0 || temp == 255)
                {
                    continue;
                }

                int index = unknownIndex[k * width + l];
                Tuple t   = tuples[index];
                if (t.flag == -1)
                {
                    continue;
                }

                double m  = mP(xi, yj, t.f, t.b);

                if (m >= minvalue[2])
                {
                    continue;
                }

                if (m < minvalue[0])
                {
                    minvalue[2] = minvalue[1];
                    p[2]   = p[1];

                    minvalue[1] = minvalue[0];
                    p[1]   = p[0];

                    minvalue[0] = m;
                    p[0].x = k;
                    p[0].y = l;
                }
                else if (m < minvalue[1])
                {
                    minvalue[2] = minvalue[1];
                    p[2]   = p[1];

                    minvalue[1] = m;
                    p[1].x = k;
                    p[1].y = l;
                }
                else if (m < minvalue[2])
                {
                    minvalue[2] = m;
                    p[2].x = k;
                    p[2].y = l;
                }
                ++num;
            }
        }

        num = min(num, 3);

        Scalar fc;
        Scalar bc;
        double sf = 0;
        double sb = 0;

        for (int k = 0; k < num; ++k)
        {
            int i  = unknownIndex[p[k].x * width + p[k].y];
            fc += tuples[i].f;
            bc += tuples[i].b;
            sf += tuples[i].sigmaf;
            sb += tuples[i].sigmab;
        }

        fc /= (num + 1e-10);
        bc /= (num + 1e-10);
        sf /= (num + 1e-10);
        sb /= (num + 1e-10);

        Scalar pc = Scalar(&data[xi * step +  yj * channels]);
        double df = distanceColor2(pc, fc);
        double db = distanceColor2(pc, bc);
        Scalar tf = fc;
        Scalar tb = bc;

        int index = xi * width + yj;
        if (df < sf)
        {
            fc = pc;
        }
        if (db < sb)
        {
            bc = pc;
        }
        if (fc.val[0] == bc.val[0] && fc.val[1] == bc.val[1] && fc.val[2] == bc.val[2])
        {
            ftuples[index].confidence = 0.00000001;
        }
        else
        {
            ftuples[index].confidence = exp(-10 * mP(xi, yj, tf, tb));
        }

        ftuples[index].f = fc;
        ftuples[index].b = bc;

        ftuples[index].alphar = max(0.0, min(1.0,comalpha(pc, fc, bc)));
    }
    tuples.clear();
}

void SharedMatting::localSmooth(unsigned char *alpha)
{
    // http://www.ruanyifeng.com/blog/2012/11/gaussian_blur.html
    double sig2 = 100.0 / (9 * 3.1415926);
    double r = 3 * sqrt(sig2);
    for(int xi = 0; xi < height; xi++)
    for(int yj = 0; yj < width; yj++)
    {
        if(m_trimap[xi * width + yj] == 0 || m_trimap[xi * width + yj] == 255)
            continue;

        int i1 = max(0, int(xi - r));
        int i2 = min(int(xi + r), height - 1);
        int j1 = max(0, int(yj - r));
        int j2 = min(int(yj + r), width - 1);

        Ftuple ptuple = ftuples[xi * width + yj];

        Scalar wcfsumup = Scalar();
        Scalar wcbsumup = Scalar();
        double wcfsumdown = 0;
        double wcbsumdown = 0;
        double wfbsumup   = 0;
        double wfbsundown = 0;
        double wasumup    = 0;
        double wasumdown  = 0;

        for (int k = i1; k <= i2; ++k)
        {
            for (int l = j1; l <= j2; ++l)
            {
                Ftuple qtuple = ftuples[k * width + l];

                double d = dP(Point(xi, yj), Point(k, l));

                if (d > r)
                {
                    continue;
                }

                double wc;
                if (d == 0)
                {
                    wc = qtuple.confidence;
                }
                else
                {
                    wc = exp(-(d * d) / sig2) * qtuple.confidence * abs(qtuple.alphar - ptuple.alphar);
                }
                wcfsumdown += wc * qtuple.alphar;
                wcbsumdown += wc * (1 - qtuple.alphar);

                wcfsumup.val[0] += wc * qtuple.alphar * qtuple.f.val[0];
                wcfsumup.val[1] += wc * qtuple.alphar * qtuple.f.val[1];
                wcfsumup.val[2] += wc * qtuple.alphar * qtuple.f.val[2];

                wcbsumup.val[0] += wc * (1 - qtuple.alphar) * qtuple.b.val[0];
                wcbsumup.val[1] += wc * (1 - qtuple.alphar) * qtuple.b.val[1];
                wcbsumup.val[2] += wc * (1 - qtuple.alphar) * qtuple.b.val[2];

                double wfb = qtuple.confidence * qtuple.alphar * (1 - qtuple.alphar);
                wfbsundown += wfb;
                wfbsumup   += wfb * sqrt(distanceColor2(qtuple.f, qtuple.b));

                double delta = 0;
                double wa;
                if (m_trimap[k * width + l] == 0 || m_trimap[k * width + l] == 255)
                {
                    delta = 1;
                }
                wa = qtuple.confidence * exp(-(d * d) / sig2) + delta;
                wasumdown += wa;
                wasumup   += wa * qtuple.alphar;
            }
        }

        Scalar cp = Scalar(&data[xi * step +  yj * channels]);
        Scalar fp;
        Scalar bp;

        double dfb;
        double conp;
        double alp;

        bp.val[0] = min(255.0, max(0.0,wcbsumup.val[0] / (wcbsumdown + 1e-200)));
        bp.val[1] = min(255.0, max(0.0,wcbsumup.val[1] / (wcbsumdown + 1e-200)));
        bp.val[2] = min(255.0, max(0.0,wcbsumup.val[2] / (wcbsumdown + 1e-200)));

        fp.val[0] = min(255.0, max(0.0,wcfsumup.val[0] / (wcfsumdown + 1e-200)));
        fp.val[1] = min(255.0, max(0.0,wcfsumup.val[1] / (wcfsumdown + 1e-200)));
        fp.val[2] = min(255.0, max(0.0,wcfsumup.val[2] / (wcfsumdown + 1e-200)));

        //double tempalpha = comalpha(cp, fp, bp);
        dfb  = wfbsumup / (wfbsundown + 1e-200);

        conp = min(1.0, sqrt(distanceColor2(fp, bp)) / dfb) * exp(-10 * mP(xi, yj, fp, bp));
        alp  = wasumup / (wasumdown + 1e-200);

        double alpha_t = conp * comalpha(cp, fp, bp) + (1 - conp) * max(0.0, min(alp, 1.0));

        alpha[xi * width + yj] = alpha_t * 255;
    }
    ftuples.clear();
}

//主干方法
void SharedMatting::solveAlpha(unsigned char *alpha)
{
    clock_t start, finish;
    //expandKnown()
    start = clock();
    cout << "Expanding...";
    expandKnown(alpha);
    cout << "    over!!!" << endl;
    finish = clock();
    cout <<  double(finish - start) / (CLOCKS_PER_SEC * 2.5) << endl;

    //gathering()
    start = clock();
    cout << "Gathering...";
    gathering();
    cout << "    over!!!" << endl;
    finish = clock();
    cout <<  double(finish - start) / (CLOCKS_PER_SEC * 2.5) << endl;

    //refineSample()
    start = clock();
    cout << "Refining...";
    refineSample(alpha);
    cout << "    over!!!" << endl;
    finish = clock();
    cout <<  double(finish - start) / (CLOCKS_PER_SEC * 2.5) << endl;

    //localSmooth()
    start = clock();
    cout << "LocalSmoothing...";
    localSmooth(alpha);
    cout << "    over!!!" << endl;
    finish = clock();
    cout <<  double(finish - start) / (CLOCKS_PER_SEC * 2.5) << endl;
}
