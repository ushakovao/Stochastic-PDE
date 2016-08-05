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
        ("dt", po::value<double>()->default_value(1.),"dt")
        ;
        return  blackscholesoptions.add( backend_options("bs1d"));
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

    //    auto u=Xh->element("u");
        auto uold=Xh->element("uold");
        auto v = Xh->element("v");
	auto u = Xh->element("u");
	

//initialisation of parameters

        double T=doption(_name="T"); //Time to maturiry
        double K=doption(_name="K"); //Strike price
        double r=doption(_name="r"); //Interest rate
        double sigmax=doption(_name="sigmax"); //Volatility for option A
        double dt = doption(_name="dt"); //Time step

//additional functions uold!!!!

        u.on(_range=elements(mesh), _expr=max(K- Px(),0));
	

//initialisation of forms
        auto l = form1(_test=Xh);
        auto a = form2(_trial=Xh, _test=Xh);

//export
        auto e = exporter(_mesh = mesh);

//iteration


	for (double t=dt; t<T; t+=dt){
	uold = u;
	l.zero();
	a=integrate(_range=elements(mesh), _expr= ((idt(u)*id(v)*(r+(1/dt)))+dxt(u)*dx(v)*(sigmax*Px())*(sigmax*Px())/2 -Px()*(r-(sigmax*sigmax)/2)dxt(u)id(v);
	a.solve(_solution=u, _rhs=l, _name="bs1d");
	e->step(t)->add("u",u);
  	e->save();
};

}
