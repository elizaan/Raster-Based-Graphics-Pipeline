#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <initializer_list>
#include "bitmap_image.hpp"
using namespace std;
#define pi (2 * acos(0.0))

struct Point
{
    double x, y, z;
};

struct Triangle
{
    struct Point points[3];
    double color[3];
};

class Matrix
{

public:
    double **arr;
    int row, column;
    Matrix();
    Matrix(double *, int r, int c);
    Matrix(int r, int c);
    Matrix(const Matrix &p1)
    {
        row = p1.row;
        column = p1.column;
        arr = new double *[row];
        for (int i = 0; i < row; ++i)
        {
            arr[i] = new double[column];
        }
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < column; j++)
            {
                arr[i][j] = p1.arr[i][j];
            }
        }
    }

    ~Matrix();

    void setMatrix(double *, int newrow, int newcolumn);
    Matrix getMatrix();
    Matrix &operator=(const Matrix &m)
    {
        if (this == &m)
        {
            return *this;
        }

        if (row != m.row || column != m.column)
        {
            for (int i = 0; i < row; ++i)
            {
                delete[] arr[i];
            }
            delete[] arr;

            row = m.row;
            column = m.column;
            arr = new double *[row];
            for (int i = 0; i < row; ++i)
            {
                arr[i] = new double[column];
            }
        }

        for (int i = 0; i < row; ++i)
        {
            for (int j = 0; j < column; ++j)
            {
                arr[i][j] = m.arr[i][j];
            }
        }
        return *this;
    }

    void setArray(double **, int r, int c);
    void setRow(int row);
    void setColumn(int column);

    double **getArray();
    int getRow();
    int getColumn();

    Matrix add(const Matrix otherMatrix);

    Matrix subtract(const Matrix otherMatrix);

    Matrix *multiply(const Matrix otherMatrix);

    void printMatrix();
    double element(int row, int column);

    // private:
    //     double **arr;
    //     int row, column;
};

Matrix::Matrix()
{
}

Matrix::Matrix(double *p, int row, int column)

{
    this->row = row;
    this->column = column;
    this->arr = new double *[row];
    for (int i = 0; i < row; ++i)
    {
        this->arr[i] = new double[column];
    }

    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < column; ++j)
        {
            int index = i * column + j;
            this->arr[i][j] = p[index];
        }
    }
}

Matrix::Matrix(int r, int c)
{
    this->row = r;
    this->column = c;
}
void Matrix::setMatrix(double *newarr, int newR, int newC)

{
    this->row = newR;
    this->column = newC;

    this->arr = new double *[newR];
    for (int i = 0; i < newR; ++i)
    {
        this->arr[i] = new double[newC];
    }

    for (int i = 0; i < newR; ++i)
    {
        for (int j = 0; j < newC; ++j)
        {
            int index = i * newC + j;
            this->arr[i][j] = newarr[index];
        }
    }
}

Matrix Matrix::getMatrix()
{
    Matrix gm({}, 0, 0);
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
        {
            gm.arr[i][j] = this->arr[i][j];
        }
    }
    gm.row = this->row;
    gm.column = this->column;
    return gm;
}
void Matrix::setRow(int r)
{
    this->row = r;
}

void Matrix::setColumn(int c)
{
    this->column = c;
}

void Matrix::setArray(double **ax, int r, int c)
{

    this->setRow(r);
    this->setColumn(c);

    arr = new double *[r];
    for (int i = 0; i < r; ++i)
    {
        arr[i] = new double[c];
    }
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {

            this->arr[i][j] = ax[i][j];
            // cout << this->arr[i][j] << endl;
        }
    }
}

int Matrix::getRow()
{
    return this->row;
}

int Matrix::getColumn()
{
    return this->column;
}

double **Matrix::getArray()
{

    double **retArr;
    retArr = new double *[this->getRow()];

    for (int h = 0; h < this->getRow(); h++)
    {
        retArr[h] = new double[this->getColumn()];
    }

    for (int h = 0; h < this->getRow(); h++)
    {

        for (int w = 0; w < this->getColumn(); w++)
        {
            // fill in some initial values
            // (filling in zeros would be more logic, but this is just for the example)
            retArr[h][w] = this->arr[h][w];
        }
    }

    return retArr;
}
Matrix Matrix::add(const Matrix p)
{

    Matrix r;

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
        {
            r.arr[i][j] = arr[i][j] + p.arr[i][j];
        }
    }

    return r;
}

Matrix Matrix::subtract(const Matrix p)
{
    Matrix r;

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
        {
            r.arr[i][j] = arr[i][j] - p.arr[i][j];
        }
    }

    return r;
}

Matrix *Matrix::multiply(const Matrix p)
{

    double **marr;
    marr = new double *[this->row];
    for (int i = 0; i < this->row; ++i)
    {
        marr[i] = new double[p.column];
    }

    for (int i = 0; i < this->row; i++)
    {
        for (int j = 0; j < p.column; j++)
        {
            marr[i][j] = 0.0;
            // cout << marr[i][j] << endl;
        }
    }

    Matrix *m = new Matrix(this->row, p.column);
    // Matrix m();
    m->setArray(marr, this->row, p.column);
    // cout << "initialized inside multiply function" << endl;
    // m.printMatrix();

    for (int i = 0; i < this->row; i++)
    {
        for (int j = 0; j < p.column; j++)
        {
            for (int k = 0; k < this->column; k++)
            {
                m->arr[i][j] = m->arr[i][j] + (this->arr[i][k] * p.arr[k][j]);
                // cout << m.arr[i][j] << endl;
            }
        }
    }

    return m;
}

void Matrix::printMatrix()
{
    cout << "The Matrix is: " << endl;
    for (int i = 0; i < this->row; i++)
    {
        for (int j = 0; j < this->column; j++)
        {
            cout << this->arr[i][j] << " ";
            // cout << "ki shmossa" << endl;
        }
        cout << endl;
    }
}

double Matrix::element(int r, int c)
{
    return this->arr[r][c];
}
Matrix::~Matrix()
{
    for (int i = 0; i < row; ++i)
    {
        delete[] arr[i];
    }
    delete[] arr;
    arr = 0;
}

int main(int argc, char const *argv[])
{

    fstream newfile;
    int stackSize;
    double eyeX, eyeY, eyeZ;
    double lookX, lookY, lookZ;
    double upX, upY, upZ;
    double fovY, aspectratio, near, far;
    stack<Matrix> S;
    double I[] = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    Matrix identity(I, 4, 4);
    cout << identity.arr[1][0] << endl;
    // identity.setMatrix(I, 4, 4);
    // identity.setMatrix(I, 4, 4);
    // identity.printMatrix();
    S.push(identity);
    S.top().printMatrix();
    Matrix *Projection = new Matrix();

    vector<Triangle> triangles;

    newfile.open("../4/scene.txt", ios::in);
    ofstream stage1("stage1.txt");
    ofstream stage2("stage2.txt");
    fstream stage3;
    stage3.open("stage3.txt", ios::out);
    // cout << "hi" << endl;

    if (newfile.is_open())
    { //checking whether the file is open
        string lines;

        for (int lineno = 0; getline(newfile, lines); lineno++)
        {
            if (lineno == 0)
            {
                sscanf(lines.c_str(), "%lf%lf%lf", &eyeX, &eyeY, &eyeZ);
            }
            else if (lineno == 1)
            {
                sscanf(lines.c_str(), "%lf%lf%lf", &lookX, &lookY, &lookZ);
            }
            else if (lineno == 2)
            {
                sscanf(lines.c_str(), "%lf%lf%lf", &upX, &upY, &upZ);
            }
            else if (lineno == 3)
            {
                sscanf(lines.c_str(), "%lf%lf%lf%lf", &fovY, &aspectratio, &near, &far);

                double fovX, t1, r1;

                fovX = fovY * aspectratio;
                fovX = fovX * (pi / 180.0);
                fovY = fovY * (pi / 180.0);
                t1 = near * tan(fovY / 2);
                r1 = near * tan(fovX / 2);
                double near1 = near / r1;
                double near2 = near / t1;
                double near3 = -(far + near) / (far - near);
                double near4 = -(2 * far * near) / (far - near);
                double pr[] = {near1, 0.0, 0.0, 0.0, 0.0, near2, 0.0, 0.0, 0.0, 0.0, near3, near4, 0.0, 0.0, -1.0, 0.0};
                Projection->setMatrix(pr, 4, 4);
            }
            else
            {

                if (lines == string("triangle"))
                {
                    double a, b, c;
                    // Matrix point({}, 0, 0);
                    Matrix *point = new Matrix();
                    Matrix *newPoint = new Matrix();
                    Matrix *Translation = new Matrix();
                    Matrix *Rotation = new Matrix();
                    Matrix *View = new Matrix();
                    Matrix *newPoint2 = new Matrix();

                    Matrix *newPoint3 = new Matrix();

                    Triangle triangle = {};

                    double lx, ly, lz;
                    double rx, ry, rz;
                    double ux, uy, uz;
                    double modulusl, modulusr;

                    cout << "triangle got" << endl;

                    for (int n = 0; n < 3; n++)
                    {
                        getline(newfile, lines);
                        lineno++;
                        sscanf(lines.c_str(), "%lf%lf%lf", &a, &b, &c);
                        double p[] = {a, b, c, 1.0};
                        point->setMatrix(p, 4, 1);
                        cout << " The " << n << "th point is: " << endl;
                        point->printMatrix();
                        cout << "After the multiplication: " << endl;
                        newPoint = S.top().multiply(*point);
                        newPoint->printMatrix();

                        stage1 << setprecision(7) << fixed << newPoint->element(0, 0) / newPoint->element(3, 0) << " " << newPoint->element(1, 0) / newPoint->element(3, 0) << " " << newPoint->element(2, 0) / newPoint->element(3, 0) << endl;

                        lx = lookX - eyeX;
                        ly = lookY - eyeY;
                        lz = lookZ - eyeZ;

                        modulusl = sqrt(lx * lx + ly * ly + lz * lz);

                        lx = lx / modulusl;
                        ly = ly / modulusl;
                        lz = lz / modulusl;
                        cout << lx << " " << ly << " " << lz << endl;

                        rx = ly * upZ - lz * upY;
                        ry = lz * upX - lx * upZ;
                        rz = lx * upY - ly * upX;

                        modulusr = sqrt(rx * rx + ry * ry + rz * rz);

                        rx = rx / modulusr;
                        ry = ry / modulusr;
                        rz = rz / modulusr;

                        ux = ry * lz - rz * ly;
                        uy = rz * lx - rx * lz;
                        uz = rx * ly - ry * lx;

                        double r[] = {rx, ry, rz, 0.0, ux, uy, uz, 0.0, -lx, -ly, -lz, 0.0, 0.0, 0.0, 0.0, 1.0};
                        Rotation->setMatrix(r, 4, 4);

                        double t[] = {1.0, 0.0, 0.0, -eyeX, 0.0, 1.0, 0.0, -eyeY, 0.0, 0.0, 1.0, -eyeZ, 0.0, 0.0, 0.0, 1.0};
                        Translation->setMatrix(t, 4, 4);

                        View = Rotation->multiply(*Translation);
                        cout << "The view matrix is: " << endl;
                        View->printMatrix();
                        newPoint2 = View->multiply(*newPoint);
                        cout << "The newpoint after view transformation is: " << endl;
                        newPoint2->printMatrix();

                        stage2 << setprecision(7) << fixed << newPoint2->element(0, 0) / newPoint2->element(3, 0) << " " << newPoint2->element(1, 0) / newPoint2->element(3, 0) << " " << newPoint2->element(2, 0) / newPoint2->element(3, 0) << endl;

                        // fovY = fovY * (pi / 180.0);

                        cout << "The projection matrix is: " << endl;
                        Projection->printMatrix();
                        newPoint3 = Projection->multiply(*newPoint2);
                        cout << "The newpoint after projection transformation is: " << endl;
                        newPoint3->printMatrix();

                        stage3 << setprecision(7) << fixed << newPoint3->element(0, 0) / newPoint3->element(3, 0) << " " << newPoint3->element(1, 0) / newPoint3->element(3, 0) << " " << newPoint3->element(2, 0) / newPoint3->element(3, 0) << endl;

                        triangle.points[n].x = newPoint3->element(0, 0) / newPoint3->element(3, 0);
                        triangle.points[n].y = newPoint3->element(1, 0) / newPoint3->element(3, 0);
                        triangle.points[n].z = newPoint3->element(2, 0) / newPoint3->element(3, 0);
                        triangle.color[n] = rand() % 255;
                    }
                    triangles.push_back(triangle);

                    stage1 << endl;
                    stage2 << endl;
                    stage3 << endl;
                }
                else if (lines == string("push"))
                {
                    // Matrix *pushMatrix = new Matrix();
                    // pushMatrix = (Matrix *)&S.top();

                    stackSize = S.size();
                    cout << "The stack size is: " << endl;
                    cout << stackSize << endl;
                    Matrix pushMatrix = S.top();
                    // pushMatrix.printMatrix();
                    S.pop();
                    S.push(pushMatrix);
                    S.top().printMatrix();

                    cout << S.top().element(0, 0) << endl;

                    cout << "push got" << endl;
                }
                else if (lines == string("pop"))
                {
                    cout << stackSize << endl;
                    cout << S.size() - stackSize << endl;
                    while (S.size() != stackSize)
                    {
                        S.pop();
                    }

                    cout << S.size() << endl;

                    cout << "pop got" << endl;
                }
                else if (lines == string("scale"))
                {

                    double sx, sy, sz;
                    Matrix *sm = new Matrix();
                    Matrix *afterS = new Matrix();
                    cout << "scale got" << endl;

                    getline(newfile, lines);
                    lineno++;
                    sscanf(lines.c_str(), "%lf%lf%lf", &sx, &sy, &sz);
                    double s[] = {sx, 0.0, 0.0, 0.0, 0.0, sy, 0.0, 0.0, 0.0, 0.0, sz, 0.0, 0.0, 0.0, 0.0, 1.0};
                    sm->setMatrix(s, 4, 4);
                    cout << " The scale matrix is: " << endl;
                    sm->printMatrix();
                    cout << "After the multiplication with top: " << endl;
                    afterS = S.top().multiply(*sm);
                    afterS->printMatrix();
                    S.push(*afterS);
                    S.top().printMatrix();
                }
                else if (lines == string("translate"))
                {
                    double tx, ty, tz;
                    Matrix *tm = new Matrix();
                    Matrix *afterT = new Matrix();
                    cout << "translate got" << endl;

                    getline(newfile, lines);
                    lineno++;
                    sscanf(lines.c_str(), "%lf%lf%lf", &tx, &ty, &tz);
                    double t[] = {1.0, 0.0, 0.0, tx, 0.0, 1.0, 0.0, ty, 0.0, 0.0, 1.0, tz, 0.0, 0.0, 0.0, 1.0};
                    tm->setMatrix(t, 4, 4);
                    cout << " The translation matrix is: " << endl;
                    tm->printMatrix();
                    cout << "After the multiplication with top: " << endl;
                    afterT = S.top().multiply(*tm);
                    afterT->printMatrix();
                    S.push(*afterT);
                    S.top().printMatrix();
                }
                else if (lines == string("rotate"))
                {
                    double angle, ax, ay, az, radian, modulus;
                    double c1x, c1y, c1z;
                    double c2x, c2y, c2z;
                    double c3x, c3y, c3z;

                    Matrix *rm = new Matrix();
                    Matrix *afterR = new Matrix();
                    cout << "rotate got" << endl;

                    getline(newfile, lines);
                    lineno++;
                    sscanf(lines.c_str(), "%lf%lf%lf%lf", &angle, &ax, &ay, &az);

                    modulus = sqrt(ax * ax + ay * ay + az * az);
                    ax = ax / modulus;
                    ay = ay / modulus;
                    az = az / modulus;

                    radian = angle * (pi / 180.0);
                    cout << radian << ax << ay << az << endl;

                    c1x = cos(radian) + (ax * ax) - (ax * ax * cos(radian));
                    c1y = (ax * ay * (1.0 - cos(radian))) + (az * sin(radian));
                    c1z = (ax * az * (1.0 - cos(radian))) - (ay * sin(radian));

                    c2x = ((1.0 - cos(radian)) * ay * ax) - (az * sin(radian));
                    c2y = cos(radian) + ((1.0 - cos(radian)) * ay * ay);
                    c2z = ((1.0 - cos(radian)) * ay * az) + (ax * sin(radian));

                    c3x = ((1.0 - cos(radian)) * ax * az) + (ay * sin(radian));
                    c3y = ((1.0 - cos(radian)) * ay * az) - (ax * sin(radian));
                    c3z = cos(radian) + ((1.0 - cos(radian)) * az * az);

                    double a[] = {c1x, c2x, c3x, 0.0, c1y, c2y, c3y, 0.0, c1z, c2z, c3z, 0.0, 0.0, 0.0, 0.0, 1.0};
                    rm->setMatrix(a, 4, 4);
                    cout << " The rotation matrix is: " << endl;
                    rm->printMatrix();
                    cout << "After the multiplication with top: " << endl;
                    afterR = S.top().multiply(*rm);
                    afterR->printMatrix();
                    S.push(*afterR);
                    S.top().printMatrix();
                }
                else if (lines == string("end"))
                {
                    cout << "end got" << endl;
                    break;
                }
                else
                {
                    cout << "lines" << endl;
                }
            }
        }
        // while (getline(newfile, tp))
        // {                       //read data from file object and put it into string.
        //     cout << tp << "\n"; //print the data of the string
        // }
        newfile.close(); //close the file object.
    }
    stage1.close();
    stage2.close();
    stage3.close();
    cout << " hi1" << endl;

    ifstream configure;

    int screen_width, screen_height;
    double Left, Bottom;
    double front, rear;
    // cout << " hi2" << endl;

    double dx, dy, Left_X, Top_Y, Bottom_Y, Right_X;

    // stage3.open("stage3.txt", ios::in);
    // cout << " hi3" << endl;
    configure.open("../4/config.txt");
    ofstream stage4("z_buffer.txt");
    // cout << " hi4" << endl;

    if (configure.is_open())
    {
        string l;
        for (int lineno = 0; getline(configure, l); lineno++)
        {
            if (lineno == 0)
            {
                sscanf(l.c_str(), "%d%d", &screen_height, &screen_width);
                cout << screen_height << " " << screen_width << endl;
            }
            else if (lineno == 1)
            {
                sscanf(l.c_str(), "%lf", &Left);
                cout << Left << endl;
            }
            else if (lineno == 2)
            {
                sscanf(l.c_str(), "%lf", &Bottom);
                cout << Bottom << endl;
            }
            else if (lineno == 3)
            {
                sscanf(l.c_str(), "%lf%lf", &front, &rear);
                cout << front << " " << rear << endl;
            }
        }

        configure.close();
    }

    double z[screen_height][screen_width];
    bitmap_image image(screen_width, screen_height);

    cout << triangles.size() << endl;
    for (int i = 0; i < triangles.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << triangles[i].points[j].x << " " << triangles[i].points[j].y << " " << triangles[i].points[j].z << endl;
            cout << triangles[i].color[j] << endl;
        }
        cout << endl;
    }

    for (int i = 0; i < screen_height; i++)
    {
        for (int j = 0; j < screen_width; j++)
        {
            z[i][j] = rear;
        }
    }

    for (int i = 0; i < screen_height; i++)
    {
        for (int j = 0; j < screen_width; j++)
        {
            image.set_pixel(j, i, 0, 0, 0);
        }
    }

    dx = (0.0 - Left - Left) / screen_width;
    dy = (0.0 - Bottom - Bottom) / screen_height;

    Top_Y = 0.0 - Bottom - (dy / 2);
    Left_X = Left + (dx / 2);
    Bottom_Y = Bottom + dy / 2;
    Right_X = 0.0 - Left - (dx / 2);

    double ymax, ymin, ys, xa, xb, za, zb, xp, zp;
    int rownomax, rownomin;

    for (int i = 0; i < triangles.size(); i++)
    {
        //clipping to checking
        ymax = max({triangles[i].points[0].y, triangles[i].points[1].y, triangles[i].points[2].y});
        // cout << ymax << endl;

        if (ymax > Top_Y)
        {
            rownomin = 0;
        }
        else
        {
            rownomin = round((Top_Y - ymax) / dy);
        }
        ymin = min({triangles[i].points[0].y, triangles[i].points[1].y, triangles[i].points[2].y});
        // cout << ymin << endl;

        if (ymin < Bottom_Y)
        {
            rownomax = screen_height;
        }
        else
        {
            rownomax = round((Top_Y - ymin) / dy);
        }

        // cout << rownomax << " " << rownomin << endl;

        //clipping to checking
        for (int j = rownomin; j <= rownomax; j++)
        {
            ys = Top_Y - (j * dy);
            double x21 = max({triangles[i].points[1].x, triangles[i].points[2].x});
            double x31 = max({triangles[i].points[1].x, triangles[i].points[3].x});
            double x32 = max({triangles[i].points[3].x, triangles[i].points[2].x});

            double xa21 = triangles[i].points[1].x + (((ys - triangles[i].points[1].y) * (triangles[i].points[2].x - triangles[i].points[1].x)) / (triangles[i].points[2].y - triangles[i].points[1].y));
            double xa31 = triangles[i].points[1].x + (((ys - triangles[i].points[1].y) * (triangles[i].points[3].x - triangles[i].points[1].x)) / (triangles[i].points[3].y - triangles[i].points[1].y));
            double xa32 = triangles[i].points[2].x + (((ys - triangles[i].points[2].y) * (triangles[i].points[3].x - triangles[i].points[2].x)) / (triangles[i].points[3].y - triangles[i].points[2].y));

            double za21 = triangles[i].points[1].z + (((ys - triangles[i].points[1].y) * (triangles[i].points[2].z - triangles[i].points[1].z)) / (triangles[i].points[2].y - triangles[i].points[1].y));
            double za31 = triangles[i].points[1].z + (((ys - triangles[i].points[1].y) * (triangles[i].points[3].z - triangles[i].points[1].z)) / (triangles[i].points[3].y - triangles[i].points[1].y));
            double za32 = triangles[i].points[2].z + (((ys - triangles[i].points[2].y) * (triangles[i].points[3].z - triangles[i].points[2].z)) / (triangles[i].points[3].y - triangles[i].points[2].y));

            double l21x = triangles[i].points[2].x - triangles[i].points[1].x;
            double l31x = triangles[i].points[3].x - triangles[i].points[1].x;
            double l32x = triangles[i].points[3].x - triangles[i].points[2].x;

            double l21y = triangles[i].points[2].y - triangles[i].points[1].y;
            double l31y = triangles[i].points[3].y - triangles[i].points[1].y;
            double l32y = triangles[i].points[3].y - triangles[i].points[2].y;

            double l21 = sqrt(l21x * l21x + l21y * l21y);
            double l31 = sqrt(l31x * l31x + l31y * l31y);
            double l32 = sqrt(l32x * l32x + l32y * l32y);

            double lxa212x = triangles[i].points[2].x - xa21;
            double lxa211x = triangles[i].points[1].x - xa21;
            double lxa212y = triangles[i].points[2].y - ys;
            double lxa211y = triangles[i].points[1].y - ys;
            double lxa212 = sqrt(lxa212x * lxa212x + lxa212y * lxa212y);
            double lxa211 = sqrt(lxa211x * lxa211x + lxa211y * lxa211y);

            double lxa313x = triangles[i].points[3].x - xa31;
            double lxa311x = triangles[i].points[1].x - xa31;
            double lxa313y = triangles[i].points[3].y - ys;
            double lxa311y = triangles[i].points[1].y - ys;
            double lxa313 = sqrt(lxa313x * lxa313x + lxa313y * lxa313y);
            double lxa311 = sqrt(lxa311x * lxa311x + lxa311y * lxa311y);

            double lxa323x = triangles[i].points[3].x - xa32;
            double lxa322x = triangles[i].points[2].x - xa32;
            double lxa323y = triangles[i].points[3].y - ys;
            double lxa322y = triangles[i].points[2].y - ys;
            double lxa323 = sqrt(lxa323x * lxa323x + lxa323y * lxa323y);
            double lxa322 = sqrt(lxa322x * lxa322x + lxa322y * lxa322y);

            vector<pair<double, double>> x;
            cout << xa21 << " " << xa31 << " " << xa32 << " " << endl;
            if (lxa323 + lxa322 - l32 < 0.1)
            {
                cout << xa32 << endl;

                x.push_back(pair<double, double>(xa32, za32));
            }
            if (lxa313 + lxa311 - l31 < 0.1)
            {

                cout << xa31 << endl;
                x.push_back(pair<double, double>(xa31, za31));
            }

            if (lxa212 + lxa211 - l21 < 0.1)
            {
                cout << xa21 << endl;

                x.push_back(pair<double, double>(xa21, za21));
            }

            // another try
            // if (xa21 >= Left_X and xa21 <= Right_X)
            // {
            //     x.push_back(pair<double, double>(xa21, za21));
            //     // x.push_back(za21);
            // }
            // if (xa31 >= Left_X and xa31 <= Right_X)
            // {
            //     x.push_back(pair<double, double>(xa31, za31));
            //     // x.push_back(xa31);
            //     // z.push_back(za31);
            // }
            // if (xa32 >= Left_X and xa32 <= Right_X)
            // {
            //     x.push_back(pair<double, double>(xa32, za32));
            //     // x.push_back(xa32);
            //     // z.push_back(za32);
            // }

            sort(x.begin(), x.end());
            cout << x.size() << endl;
            int colnomin, colnomax;
            //clipping checking
            xa = x[0].first;
            if (xa < Left_X)
            {
                colnomin = 0;
            }
            else
            {
                colnomin = round((xa - Left_X) / dx);
            }
            xb = x[1].first;

            if (xb > Right_X)
            {
                colnomax = screen_width;
            }
            else
            {

                colnomax = round((xb - Left_X) / dx);
            }

            // cout << colnomax << " " << colnomin << endl;

            //clipping checking

            for (int k = colnomin; k <= colnomax; k++)
            {
                xp = Left_X + k * dx;
                zp = x[0].second + (((xp - xa) * (x[1].second - x[0].second)) / (xb - xa));

                if ((zp >= front) && (zp < z[j][k]))
                {
                    z[j][k] = zp;
                    stage4 << setprecision(6) << fixed << z[j][k] << " ";
                    image.set_pixel(k, j, triangles[i].color[0], triangles[i].color[1], triangles[i].color[2]);
                }
                else
                {
                    stage4 << " ";
                    // z[j][k] = " ";
                }
            }

            stage4 << endl;
        }
    }

    stage4.close();
    cout << "hi" << endl;
    image.save_image("output.bmp");

    // image.save_image("test.bmp");

    return 0;
}