#include "SimConfig.h"

#include <iomanip>

bool SimConfig::writeFile(const std::string& json_file) {
	std::ofstream o(json_file);
	if (o.fail()) {
		std::cout << "*** Warning: could not write configuration file '"
		          << json_file << "'.\n" << std::endl;
		return false;
	} else {
		std::cout << "Writing configuration file '"
		          << json_file << "'.\n" << std::endl;
		o << std::setw(4) << m_config << std::endl;
		return true;
	}
}

bool SimConfig::readFile(const std::string &json_file) {
	std::ifstream infile(json_file);
	if (infile.fail()) {
		std::cout << "*** Warning: could not read configuration file '"
		          << json_file << "'. Reverting to defaults.\n" << std::endl;
		setDefaults(true);
		return false;
	} else {
		std::cout << "Reading configuration file '"
		          << json_file << "'.\n" << std::endl;
		m_config = json();
		infile >> m_config;
		setDefaults();

		return true;
	}
}

void SimConfig::setDefaults(bool hard) {
	if (hard)
		m_config = json();

	if (!m_config.contains("exportMeshes"))
		setExportMeshes(false);
	if (!m_config.contains("systemSize"))
		setSystemSize(120, 120, 120);
	if (!m_config.contains("gridResolution"))
		setGridResolution(15, 15, 15);
	if (!m_config.contains("jitterParticles"))
		setJitterParticles(true);
	if (!m_config.contains("timeStep"))
		setTimeStep(0.025);
	if (!m_config.contains("alpha"))
		setAlpha(0.01);
	if (!m_config.contains("density"))
		setDensity(1000.0);
	if (!m_config.contains("gravity"))
		setGravity(9.81);
	if (!m_config.contains("displayGrid"))
		setDisplayGrid(false);
	if (!m_config.contains("displayMeshes"))
		setDisplayMeshes(true, false);
	if (!m_config.contains("maxParticlesDisplay"))
		setMaxParticlesDisplay(424242);
	if (!m_config.contains("maxSteps"))
		setMaxSteps(-1);
	if (!m_config.contains("applyMeteorForce"))
		setApplyMeteorForce(false);
	if (!m_config.contains("fluidRegion"))
		setFluidRegion(22, 22, 22, 90, 110, 90);
	if (!m_config.contains("randomSeed"))
		setRandomSeed(-1);
}

void SimConfig::setExportMeshes(bool v) {
	m_config["exportMeshes"] = v;
}

bool SimConfig::getExportMeshes() const {
	return m_config["exportMeshes"];
}

void SimConfig::setSystemSize(double x, double y, double z) {
	m_config["systemSize"] = {x, y, z};
}

void SimConfig::getSystemSize(double &x, double &y, double &z) const {
	x = m_config["systemSize"][0];
	y = m_config["systemSize"][1];
	z = m_config["systemSize"][2];
}

void SimConfig::setGridResolution(int x, int y, int z) {
	m_config["gridResolution"] = {x, y, z};
}

void SimConfig::getGridResolution(int &x, int &y, int &z) const {
	x = m_config["gridResolution"][0];
	y = m_config["gridResolution"][1];
	z = m_config["gridResolution"][2];
}

void SimConfig::setJitterParticles(bool jitter) {
	m_config["jitterParticles"] = jitter;
}

bool SimConfig::getJitterParticles() const {
	return m_config["jitterParticles"];
}

void SimConfig::setTimeStep(double dt) {
	m_config["timeStep"] = dt;
}

double SimConfig::getTimeStep() const {
	return m_config["timeStep"];
}

void SimConfig::setAlpha(double alpha) {
	m_config["alpha"] = alpha;
}

double SimConfig::getAlpha() const {
	return m_config["alpha"];
}

void SimConfig::setDensity(double density) {
	m_config["density"] = density;
}

double SimConfig::getDensity() const {
	return m_config["density"];
}

void SimConfig::setGravity(double gravity) {
	m_config["gravity"] = gravity;
}

double SimConfig::getGravity() const {
	return m_config["gravity"];
}

void SimConfig::setDisplayGrid(bool display) {
	m_config["displayGrid"] = display;
}

bool SimConfig::getDisplayGrid() const {
	return m_config["displayGrid"];
}

void SimConfig::setDisplayMeshes(bool displayEdges, bool displayFaces) {
	m_config["displayMeshes"] = {displayEdges, displayFaces};
}

bool SimConfig::getDisplayMeshes() const {
	return m_config["displayMeshes"][0] || m_config["displayMeshes"][1];
}

bool SimConfig::getDisplayMeshEdges() const {
	return m_config["displayMeshes"][0];
}

bool SimConfig::getDisplayMeshFaces() const {
	return m_config["displayMeshes"][1];
}

void SimConfig::setMaxParticlesDisplay(int maxParticles) {
	m_config["maxParticlesDisplay"] = maxParticles;
}

int SimConfig::getMaxParticlesDisplay() const {
	return m_config["maxParticlesDisplay"];
}

void SimConfig::setMaxSteps(int maxSteps) {
	m_config["maxSteps"] = maxSteps;
}

int SimConfig::getMaxSteps() const {
	return m_config["maxSteps"];
}

void SimConfig::setApplyMeteorForce(bool explode) {
	m_config["applyMeteorForce"] = explode;
}

bool SimConfig::getApplyMeteorForce() const {
	return m_config["applyMeteorForce"];
}

void SimConfig::setFluidRegion(double from_x, double from_y, double from_z, double to_x,
							   double to_y, double to_z) {
	simpleSort(from_x, to_x);
	simpleSort(from_y, to_y);
	simpleSort(from_z, to_z);
	m_config["fluidRegion"] = {{from_x, from_y, from_z}, {to_x, to_y, to_z}};
}

void SimConfig::getFluidRegion(double &from_x, double &from_y, double &from_z,
							   double &to_x, double &to_y, double &to_z) const {
	from_x = m_config["fluidRegion"][0][0];
	from_y = m_config["fluidRegion"][0][1];
	from_z = m_config["fluidRegion"][0][2];
	to_x = m_config["fluidRegion"][1][0];
	to_y = m_config["fluidRegion"][1][1];
	to_z = m_config["fluidRegion"][1][2];
	simpleSort(from_x, to_x);
	simpleSort(from_y, to_y);
	simpleSort(from_z, to_z);
}

void SimConfig::setRandomSeed(int seed) {
	m_config["randomSeed"] = seed;
}

int SimConfig::getRandomSeed() {
	return m_config["randomSeed"];
}