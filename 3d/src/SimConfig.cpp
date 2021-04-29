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
		setGridResolution(40, 40, 40);
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
	if (!m_config.contains("maxParticlesDisplay"))
		setMaxParticlesDisplay(4242);
}

void SimConfig::setExportMeshes(bool v) {
	m_config["exportMeshes"] = v;
}

bool SimConfig::getExportMeshes() {
	return m_config["exportMeshes"];
}

void SimConfig::setSystemSize(double x, double y, double z) {
	m_config["systemSize"] = {x, y, z};
}

void SimConfig::getSystemSize(double &x, double &y, double &z) {
	x = m_config["systemSize"][0];
	y = m_config["systemSize"][1];
	z = m_config["systemSize"][2];
}

void SimConfig::setGridResolution(int x, int y, int z) {
	m_config["gridResolution"] = {x, y, z};
}

void SimConfig::getGridResolution(int &x, int &y, int &z) {
	x = m_config["gridResolution"][0];
	y = m_config["gridResolution"][1];
	z = m_config["gridResolution"][2];
}

void SimConfig::setJitterParticles(bool jitter) {
	m_config["jitterParticles"] = jitter;
}

bool SimConfig::getJitterParticles() {
	return m_config["jitterParticles"];
}

void SimConfig::setTimeStep(double dt) {
	m_config["timeStep"] = dt;
}

double SimConfig::getTimeStep() {
	return m_config["timeStep"];
}

void SimConfig::setAlpha(double alpha) {
	m_config["alpha"] = alpha;
}

double SimConfig::getAlpha() {
	return m_config["alpha"];
}

void SimConfig::setDensity(double density) {
	m_config["density"] = density;
}

double SimConfig::getDensity() {
	return m_config["density"];
}

void SimConfig::setGravity(double gravity) {
	m_config["gravity"] = gravity;
}

double SimConfig::getGravity() {
	return m_config["gravity"];
}

void SimConfig::setDisplayGrid(bool display) {
	m_config["displayGrid"] = display;
}

bool SimConfig::getDisplayGrid() {
	return m_config["displayGrid"];
}

void SimConfig::setMaxParticlesDisplay(int maxParticles) {
	m_config["maxParticlesDisplay"] = maxParticles;
}

int SimConfig::getMaxParticlesDisplay() {
	return m_config["maxParticlesDisplay"];
}