#include "include/Mac2d.cpp"

int main(){
	Mac2d *grid = new Mac2d(4,3,0.1,0.1);	
	std::cout << (grid->index_from_coord(5.6, 7.8)).first << "     " << (grid->index_from_coord(5.6, 7.8)).second << std::endl;
	
	return 0;
}
