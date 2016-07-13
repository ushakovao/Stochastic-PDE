#include <feel/feel.hpp>

using namespace Feel;
inline
po::options_description
makeOptions()
{
po::options_description blackscholesoptions ("BS options");
blackscholesoptions.add_options()
        ("T", po::value<double>()->default_value(1.), "Time to maturity")
        ("K", po::value<double>()->default_value(1.), "Strike price")
        ("r", po::value<double>()->default_value(1.),"Interest rate")
        ("sigmax", po::value<double>()->default_value(1.),"sigma x")
        ("sigmay", po::value<double>()->default_value(1.),"sigma y")
        ("rho", po::value<double>()->default_value(1.),"rho")
        ("mu", po::value<double>()->default_value(1.),"mu")
        ("dt", po::value<double>()->default_value(1.),"dt")
        ;
        return  blackscholesoptions.add( backend_options("bs"));
}

int main (int argc, char* argv[])
{

        Environment env(_argc = argc, _argv=argv,
        _desc=makeOptions(),
        _about=about(_name="mymesh",
                     _author="Feel++ Consortium",
                     _email="feelpp-devel@feelpp.org"));

//loading a created mesh

        auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);

//mesh adaptation - TODO

//determine space

        auto Xh = Pch<1>(mesh);

//initialisation of elements

        auto u=Xh->element("u");
        auto uold=Xh->element("uold");
        auto v = Xh->element("v");

//initialisation of parameters

        double T=doption(_name="T"); //Time to maturiry
        double K=doption(_name="K"); //Strike price
        double r=doption(_name="r"); //Interest rate
        double sigmax=doption(_name="sigmax"); //Volatility for option A
        double sigmay=doption(_name="sigmay"); //Volatility for option B
        double rho = doption(_name="rho"); //Correlation between A and B
        double mu = doption(_name="mu"); //Drift rate
        double dt = doption(_name="dt"); //Time step

//additional functions

        auto u0=max(K-max( Px(), Py() ),0);
        auto xvel = -Px()*r + Px()*sigmax*sigmax+Px()*rho*sigmax*sigmay/2;
        auto yvel = -Py()*r + Py()*sigmay*sigmay+Py()*rho*sigmax*sigmay/2;

//initialisation of forms
        auto l = form1(_test=Xh);
        auto a = form2(_trial=Xh, _test=Xh);

//export
        auto e = exporter(_mesh = mesh);



//iteration

        for (double t=dt; t<T; t+=dt){
        l.zero();
        a=integrate(_range=elements(mesh),_expr =( (gradt(u)*vec(xvel,yvel))*id(v)));
        a+=integrate(_range=elements(mesh), _expr= ((idt(u)*id(v)*(r+(1/dt)))+dxt(u)*dx(v)*(sigmax*Px())*(sigmax*Px())/2+dyt(u)*dy(v)*(sigmay*Py())*(sigmay*Py())/2+(dyt(u)*dx(v) + dxt(u)*dy(v))*rho*sigmax*sigmay*Px()*Py()/2) );

        a+=on(_range=markedfaces(mesh, "out"), _rhs=l, _element=u,_expr=cst(0) );
        a+=on(_range=markedfaces(mesh,"in"), _rhs=l, _element=u,_expr=u0);
        a.solve(_solution=u, _rhs=l, _name="bs");
        uold = u;
 	e->step(t)->add("u",u);
        e->save();



};
}



