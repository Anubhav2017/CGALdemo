#include <QCoreApplication>
#include<iostream>
#include <CGAL/Simple_cartesian.h>
#include<CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K:: Vector_3 Vector;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef Mesh::vertex_index vi;
typedef Mesh::face_index fi;
using namespace std;
using namespace CGAL;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    Point p(1,3,2), q(10,5,1), r(4,-2,8);
    Vector vector1(p-q);
    Mesh m;
    vi v1=m.add_vertex(p);
    vi v2=m.add_vertex(q);
    vi v3=m.add_vertex(r);
    cout<<vector1.x()<<'\n';
    cout<<p.x()<<'\n';
    fi f1=m.add_face(v1,v2,v3);
    CGAL::draw(m);


    return a.exec();
}
