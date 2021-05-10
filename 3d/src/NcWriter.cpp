#include "NcWriter.h"

NcWriter::NcWriter( const std::string& filePath, 
					const unsigned break_points, 
					const unsigned num_particles, 
					const unsigned n, 
					const unsigned m, 
					const unsigned l ):
					_filePath(filePath)
{

	try{
		
		NcFile dataFile(_filePath, NcFile::replace);

		dataFile.addDim("break_points", break_points);
		dataFile.addDim("num_particles", num_particles);
		dataFile.addDim("x_mac_faces", (n+1) *  m    *  l   );
		dataFile.addDim("y_mac_faces",  n    * (m+1) *  l   );
		dataFile.addDim("z_mac_faces",  n    *  m    * (l+1));
		dataFile.addDim("mac_centers",  n * m * l);

		dataFile.close();

	} catch(NcException& e){
		std::cout << "--- NCWRITER CONSTRUCTOR FAILURE ---" << std::endl;
		std::cout << e.what() << std::endl;
	}

}


NcWriter::~NcWriter(){
	
}


void NcWriter::addVar( const std::string& varName, const std::string& dimName, NcType ncType ){
	
	try{
		
		// Open the file to write
		NcFile dataFile(_filePath, NcFile::write);
		
		// Get the dimensions
		NcDim bpDim = dataFile.getDim("break_points");
		NcDim nDim  = dataFile.getDim(dimName);
		
		// Store the dimensions into a vector in the correct order
		std::vector<NcDim> dims;
		dims.push_back(bpDim);
		dims.push_back(nDim);
		
		// Add the new variable to the file
		NcVar var = dataFile.addVar(varName, ncType, dims);
		
		// Define the attributes of the new variable
		// var.putAtt("units", "m");
		
	} catch(NcException& e){
		std::cout << "--- ADD VAR FAILURE: " + varName + " ---" << std::endl;
		std::cout << e.what() << std::endl;
	}
}


void NcWriter::toLinArrays( Particle* particles, 
							unsigned num_particles, 
							Mac3d* MACGrid, 
							unsigned n, 
							unsigned m, 
							unsigned l, 
							double* x, 
							double* y, 
							double* z, 
							double* u, 
							double* v, 
							double* w, 
							double* uMAC, 
							double* vMAC, 
							double* wMAC, 
							double* uStar, 
						 	double* vStar, 
						 	double* wStar, 
							double* pMAC, 
							bool* fluid_cells, 
							bool* solid_cells )
{
	
	Eigen::Vector3d pos;
	Eigen::Vector3d vel;

	for(unsigned i = 0; i < num_particles; ++i){

		pos = (particles + i)->get_position();
		vel = (particles + i)->get_velocity();

		x[i] = pos(0);
		y[i] = pos(1);
		z[i] = pos(2);
		u[i] = vel(0);
		v[i] = vel(1);
		w[i] = vel(2);
	}

	for(unsigned i = 0; i < n+1; ++i){
		for(unsigned j = 0; j < m+1; ++j){
			for(unsigned k = 0; k < l+1; ++k){

				if( j < m && k < l ){
					uMAC [i + (n+1) * j + (n+1) *  m    * k] = MACGrid->get_u(i, j, k);
					uStar[i + (n+1) * j + (n+1) *  m    * k] = MACGrid->get_u_star(i, j, k);
				}
				if( i < n && k < l ){
					vMAC [i +  n    * j +  n    * (m+1) * k] = MACGrid->get_v(i, j, k);
					vStar[i +  n    * j +  n    * (m+1) * k] = MACGrid->get_v_star(i, j, k);
				}
				if( i < n && j < m ){
					wMAC [i +  n    * j +  n    *  m    * k] = MACGrid->get_w(i, j, k);
					wStar[i +  n    * j +  n    *  m    * k] = MACGrid->get_w_star(i, j, k);
				}
				
				if( i < n && j < m && k < l ){
					pMAC[i + n*j + n*m*k] = MACGrid->get_pressure(i, j, k);
					fluid_cells[i + n*j + n*m*k] = MACGrid->is_fluid(i, j, k);
					solid_cells[i + n*j + n*m*k] = MACGrid->is_solid(i, j, k);
				}
			}
		}
	}

}


void NcWriter::writeAll( unsigned breakPt, 
						 double* x, 
						 double* y, 
						 double* z, 
						 double* u, 
						 double* v, 
					 	 double* w, 
				 		 double* uMAC, 
						 double* vMAC, 
						 double* wMAC, 
						 double* uStar, 
						 double* vStar, 
						 double* wStar, 
						 double* pMAC, 
						 bool* fluid_cells, 
						 bool* solid_cells )
{

	write(breakPt, "x", x);
	write(breakPt, "y", y);
	write(breakPt, "z", z);
	write(breakPt, "u", u);
	write(breakPt, "v", v);
	write(breakPt, "w", w);
	write(breakPt, "uMAC", uMAC);
	write(breakPt, "vMAC", vMAC);
	write(breakPt, "wMAC", wMAC);
	write(breakPt, "uStar", uStar);
	write(breakPt, "vStar", vStar);
	write(breakPt, "wStar", wStar);
	write(breakPt, "pMAC", pMAC);
	write(breakPt, "fluid_cells", fluid_cells);
	write(breakPt, "solid_cells", solid_cells);

}
