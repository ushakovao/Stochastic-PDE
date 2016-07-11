#include <feel/feel.hpp>

int main(int argc, char* argv[])
{
  using namespace Feel;
  
  Environment env( _argc=argc, _argv=argv,
    _der(_meshsc=feel_options(),
    _about=about( _name="mymesh",
                  _author="Feel++ Consortium",
                  _email="feelpp-devel@feelpp.org") );


  // create mesh and exporting it (ADD GEO IN PROJECT FOLDER!!! )
	
  auto mesh = loadMesh( _mesh=new Mesh<Simplex<DIM>>,
                          _filename="square.geo" );


	//Mesh Adaptation TODO
	int adapt_iter = 0;
	int max_iter = 10;


	//determine space
  auto Xh = Pch<1>(mesh);
	
	//initialisation of elements
  auto u = Xh->element("u");
	auto uold=Xh->element("u");
  auto v = Xh->element("v");
	auto vold = Xh ->element("v");

  // parameters
  double T=doption(_name="T");//Time to maturity
  int K= doption(_name="K");// Strike price
  double r=doption(_name="r"); //Risk-free rate
  double sigmax = doption(_name="sigmax");//Volatility for optionA
	double sigmay = doption(_name="sigmay");//Volaility for optionB
	double rho = doption(_name="rho");//Correlation 
  double mu = doption(_name="mu"); //Drift rate 
  double dt =doption(_name="dt");

  double alpha = doption(_name="alpha");
  double beta = doption(_name="beta");
  double gamma = doption(_name="gamma");
  double eta = doption(_name="eta");

  int Nmax = doption(_name="Nmax");
  double scale=doption(_name="scale");
  
  auto u0 = max(K-max(Px(),Py()), 0); 
	
	auto xvel = -Px()*r + Px()*sigmax^2+Px()*rho*sigmax*sigmay/2;
	auto yvel = -Py()*r + Py()*sigmay^2+Py()*rho*sigmax*sigmay/2;


//Introduction of forms
	auto l = form1(_test=Xh);
	auto a = form2(_trial=Xh, _test=Xh);

//iterating
	for (int i=0, i<0.5/dt; i++){
		l = u0;
	
		a =  integrate (_range=elements(mesh),_expr = ((idt(u)*id(v)*(r+(1/dt))) + 
																									dx(u)*dx(v) *(sigmax*Px())^2)/2+ 
																									dy(u)*dy(v)*(sigmay*Py())^2)/2 + 
																									(dy(u)*dx(v) + dx(u)*dy(v))*rho*sigmax*sigmay*Px()*Py()/2) - 
																									((gradt(u)*vec(xvel, yvel))*id(v)/dt;
		a+=on(_range=markedfaces(mesh,"out"), _rhs=l1, _element=u, _expr=cst(0.) );

		a.solve( _solution=u, _rhs=l, _name="BS" );

		uold = u;
	}
}
