#include "NcReader.h"

NcReader::NcReader( const std::string& filePath, const std::string& cfgPath ):
					_filePath(filePath) 
{

	cfg = SimConfig(cfgPath);
	cfg.getGridResolution(_n, _m, _l);
	cfg.getSystemSize(_dx, _dy, _dz);

	num_particles = getDim("num_particles");
	timestep      = readScalar<unsigned>("timestep");

	particles = new Particle[num_particles];
	MACGrid   = new Mac3d(_n, _m, _l, _dx, _dy, _dz);

	unsigned cacheBlockSize = 64;
	x     = (double*) aligned_alloc(cacheBlockSize, num_particles    * sizeof(double));
	y     = (double*) aligned_alloc(cacheBlockSize, num_particles    * sizeof(double));
	z     = (double*) aligned_alloc(cacheBlockSize, num_particles    * sizeof(double));
	u     = (double*) aligned_alloc(cacheBlockSize, num_particles    * sizeof(double));
	v     = (double*) aligned_alloc(cacheBlockSize, num_particles    * sizeof(double));
	w     = (double*) aligned_alloc(cacheBlockSize, num_particles    * sizeof(double));
	uMAC  = (double*) aligned_alloc(cacheBlockSize, (_n+1) * _m * _l * sizeof(double));
	vMAC  = (double*) aligned_alloc(cacheBlockSize, _n * (_m+1) * _l * sizeof(double));
	wMAC  = (double*) aligned_alloc(cacheBlockSize, _n * _m * (_l+1) * sizeof(double));
	uStar = (double*) aligned_alloc(cacheBlockSize, (_n+1) * _m * _l * sizeof(double));
	vStar = (double*) aligned_alloc(cacheBlockSize, _n * (_m+1) * _l * sizeof(double));
	wStar = (double*) aligned_alloc(cacheBlockSize, _n * _m * (_l+1) * sizeof(double));
	pMAC  = (double*) aligned_alloc(cacheBlockSize, _n * _m * _l     * sizeof(double));

	fluid_cells = (bool*) aligned_alloc(cacheBlockSize, _n * _m * _l * sizeof(bool));
	solid_cells = (bool*) aligned_alloc(cacheBlockSize, _n * _m * _l * sizeof(bool));
}


NcReader::~NcReader(){
	
	delete[] particles;
	delete MACGrid;
	free(x);
	free(y);
	free(z);
	free(u);
	free(v);
	free(w);
	free(uMAC);
	free(vMAC);
	free(wMAC);
	free(uStar);
	free(vStar);
	free(wStar);
	free(pMAC);
	free(fluid_cells);
	free(solid_cells);
}


void NcReader::toOldStruct(){

	for(unsigned i = 0; i < num_particles; ++i){
		
		particles[i].set_position(x[i], y[i], z[i]);
		particles[i].set_velocity(u[i], v[i], w[i]);
	}
	
	for(unsigned i = 0; i < _n+1; ++i){
		for(unsigned j = 0; j < _m+1; ++j){
			for(unsigned k = 0; k < _l+1; ++k){

				if( j < _m && k < _l ){
					MACGrid->set_u(i, j, k, uMAC[i + (_n+1) * j + (_n+1) *  _m    * k]);
					MACGrid->set_u_star(i, j, k, uStar[i + (_n+1) * j + (_n+1) *  _m    * k]);
				}
				if( i < _n && k < _l ){
					MACGrid->set_v(i, j, k, vMAC[i +  _n    * j +  _n    * (_m+1) * k]);
					MACGrid->set_v_star(i, j, k, vStar[i +  _n    * j +  _n    * (_m+1) * k]);
				}
				if( i < _n && j < _m ){
					MACGrid->set_w(i, j, k, wMAC[i +  _n    * j +  _n    *  _m    * k]);
					MACGrid->set_w_star(i, j, k, wStar[i +  _n    * j +  _n    *  _m    * k]);
				}
				
				if( i < _n && j < _m && k < _l ){
					MACGrid->set_pressure(i, j, k, pMAC[i + _n*j + _n*_m*k]);
					if(fluid_cells[i + _n*j + _n*_m*k]) MACGrid->set_fluid(i, j, k);
					if(solid_cells[i + _n*j + _n*_m*k]) MACGrid->set_solid(i, j, k);
				}
			}
		}
	}
}


void NcReader::readAll( unsigned breakPt ){
	
	read(breakPt, "x", x);
	read(breakPt, "y", y);
	read(breakPt, "z", z);
	read(breakPt, "u", u);
	read(breakPt, "v", v);
	read(breakPt, "w", w);
	read(breakPt, "uMAC", uMAC);
	read(breakPt, "vMAC", vMAC);
	read(breakPt, "wMAC", wMAC);
	read(breakPt, "uStar", uStar);
	read(breakPt, "vStar", vStar);
	read(breakPt, "wStar", wStar);
	read(breakPt, "pMAC", pMAC);
	read(breakPt, "fluid_cells", fluid_cells);
	read(breakPt, "solid_cells", solid_cells);
}


double NcReader::rErr( const double ref, const double value ){

	if( ref == value ){
		return 0.;
	}
	else if( ref != 0. ){
		return std::abs(value - ref) / ref;
	}
	else{
		return 1.;
	}
}


void NcReader::outputMessage( const bool flag, const std::string& varName, const double maxErr, const double tol ){

	if(!flag) std::cout << ">> " << varName << " does not validate! The maximum relative error is " << maxErr << " >= " << tol << std::endl;
}


void NcReader::validate(){

	double tol = 1e-12;
	double maxErrX     = 0.;
	double maxErrY     = 0.;
	double maxErrZ     = 0.;
	double maxErrU     = 0.;
	double maxErrV     = 0.;
	double maxErrW     = 0.;
	double maxErrUMAC  = 0.;
	double maxErrVMAC  = 0.;
	double maxErrWMAC  = 0.;
	double maxErrUStar = 0.;
	double maxErrVStar = 0.;
	double maxErrWStar = 0.;
	double maxErrPMAC  = 0.;

	bool fFlag = true;
	bool sFlag = true;

	Eigen::Vector3d pos, vel;
	for(unsigned i = 0; i < num_particles; ++i){

		pos = particles[i].get_position();
		vel = particles[i].get_velocity();

		maxErrX = std::max(maxErrX, rErr(x[i], pos(0)));
		maxErrY = std::max(maxErrY, rErr(y[i], pos(1)));
		maxErrZ = std::max(maxErrZ, rErr(z[i], pos(2)));
		maxErrU = std::max(maxErrU, rErr(u[i], vel(0)));
		maxErrV = std::max(maxErrV, rErr(v[i], vel(1)));
		maxErrW = std::max(maxErrW, rErr(w[i], vel(2)));
	}

	for(unsigned i = 0; i < _n+1; ++i){
		for(unsigned j = 0; j < _m+1; ++j){
			for(unsigned k = 0; k < _l+1; ++k){

				if( j < _m && k < _l ){
					maxErrUMAC  = std::max(maxErrUMAC,  rErr(uMAC [i + (_n+1) * j + (_n+1) *  _m    * k], MACGrid->get_u(i, j, k)));
					maxErrUStar = std::max(maxErrUStar, rErr(uStar[i + (_n+1) * j + (_n+1) *  _m    * k], MACGrid->get_u_star(i, j, k)));
				}
				if( i < _n && k < _l ){
					maxErrVMAC  = std::max(maxErrVMAC,  rErr(vMAC [i +  _n    * j +  _n    * (_m+1) * k], MACGrid->get_v(i, j, k)));
					maxErrVStar = std::max(maxErrVStar, rErr(vStar[i +  _n    * j +  _n    * (_m+1) * k], MACGrid->get_v_star(i, j, k)));
				}
				if( i < _n && j < _m ){
					maxErrWMAC  = std::max(maxErrWMAC,  rErr(wMAC [i +  _n    * j +  _n    *  _m    * k], MACGrid->get_w(i, j, k)));
					maxErrWStar = std::max(maxErrWStar, rErr(wStar[i +  _n    * j +  _n    *  _m    * k], MACGrid->get_w_star(i, j, k)));
				}
				
				if( i < _n && j < _m && k < _l ) {
					maxErrPMAC  = std::max(maxErrPMAC,  rErr(pMAC [i +  _n    * j +  _n    *  _m    * k], MACGrid->get_pressure(i, j, k)));
					fFlag = fFlag && (fluid_cells[i + _n*j + _n*_m*k] == MACGrid->is_fluid(i, j, k));
					sFlag = sFlag && (solid_cells[i + _n*j + _n*_m*k] == MACGrid->is_solid(i, j, k));
				}
			}
		}
	}

	bool xFlag     = maxErrX < tol;
	bool yFlag     = maxErrY < tol;
	bool zFlag     = maxErrZ < tol;
	bool uFlag     = maxErrU < tol;
	bool vFlag     = maxErrV < tol;
	bool wFlag     = maxErrW < tol;
	bool uMACFlag  = maxErrUMAC < tol;
	bool vMACFlag  = maxErrVMAC < tol;
	bool wMACFlag  = maxErrWMAC < tol;
	bool uStarFlag = maxErrUStar < tol;
	bool vStarFlag = maxErrVStar < tol;
	bool wStarFlag = maxErrWStar < tol;
	bool pMACFlag  = maxErrPMAC < tol;

	outputMessage(xFlag,     "x",     maxErrX,     tol);
	outputMessage(yFlag,     "y",     maxErrY,     tol);
	outputMessage(zFlag,     "z",     maxErrZ,     tol);
	outputMessage(uFlag,     "u",     maxErrU,     tol);
	outputMessage(vFlag,     "v",     maxErrV,     tol);
	outputMessage(wFlag,     "w",     maxErrW,     tol);
	outputMessage(uMACFlag,  "uMAC",  maxErrUMAC,  tol);
	outputMessage(vMACFlag,  "vMAC",  maxErrVMAC,  tol);
	outputMessage(wMACFlag,  "wMAC",  maxErrWMAC,  tol);
	outputMessage(uStarFlag, "uStar", maxErrUStar, tol);
	outputMessage(vStarFlag, "vStar", maxErrVStar, tol);
	outputMessage(wStarFlag, "wStar", maxErrWStar, tol);
	outputMessage(pMACFlag,  "pMAC",  maxErrPMAC,  tol);

	assert(xFlag);
	assert(yFlag);
	assert(zFlag);
	assert(uFlag);
	assert(vFlag);
	assert(wFlag);
	assert(uMACFlag);
	assert(vMACFlag);
	assert(wMACFlag);
	assert(uStarFlag);
	assert(vStarFlag);
	assert(wStarFlag);
	assert(pMACFlag);
	assert(fFlag);
	assert(sFlag);
}


unsigned NcReader::getDim( const std::string& dimName ){
	
	try{
		
		// Open the file to read
		NcFile dataFile(_filePath, NcFile::read);
		
		return dataFile.getDim(dimName).getSize();

	} catch(NcException& e){
		std::cout << "--- GET DIM FAILURE: " + dimName + " ---" << std::endl;
		std::cout << e.what() << std::endl;

		return 0;
	}
}