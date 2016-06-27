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
	
  int Lx = doption(_name="Lx");
  int Ly = doption(_name="Ly");

	mesh_ptrtype mesh = createGMSHMesh(_mesh = new mesh_type,
																		_desc = createBSGeo( meshsize, Lx, Ly),
																		_update=MESH_UPDATE_FACES | MESH_UPDATE_EDGES);

	int adapt_iter = 0;
	int max_iter = 10;

	//Mesh Adaptation
/*
	MeshAdapt mash_adaption;
	do{
		 space_ptrtype Xh = space_type::New( mesh );
        u = Xh->element();
        v = Xh->element();


*/


/*
  auto mesh = unitSquare();
  auto exp=exporter(_mesh=mesh);
  e->addRegions();
  e->save();
*/


	



	//determine space
  auto Xh = Pch<1>(mesh);
	
	//initialisation of forms
	auto a  = form2( _trial=Xh,_test=Xh);
	auto l = form1(_test = Xh);

  auto u = Xh->element("u");
  auto v = Xh->element("v");
	

  // parameters
  double T=doption(_name="T");
  int K= doption(_name="K");// Strike price

  double r=doption(_name="r") //risk-free rate
  double sigmax = doption(_name="sigmax");
	double sigmay = doption(_name="sigmay");
	double rho = doption(_name="rho");
	double dt =T/Nmax;
  double mu = doption(_name="mu"); //drift rate 
  

  double alpha = doption(_name="alpha");
  double beta = doption(_name="beta");
  double gamma = doption(_name="gamma");
  double eta = doption(_name="eta");

  int Nmax = doption(_name="Nmax");
  double scale=doption(_name="scale");
  
  auto u0 = max(K-x, 0); 



	auto l = form1(_test=Xh);
	l = integrate (_range = elements(mesh), _expr= (idt(x)*id(x)*(r+1/dt));

	auto a  =form2(_trial=Xh, _test=Xh);
	a = integrate (_range=elements(mesh),_expr = ((sigma*Px()/2)^2)*dx(u)*dx(v) + Px()*dx(u)*id(v)*((sigma^2)-r)+r*idt(u)*id(v));

	a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=cst(0.));


	auto xvel = -Px()*r + Px()*sigmax^2+Px()*rho*sigmax*sigmay/2;
	auto yvel = -Py()*r + Py()*sigmay^2+Py()*rho*sigmax*sigmay/2;



